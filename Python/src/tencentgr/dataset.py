"""Local TencentGR-10M subset loader. Persist seq/item/user/indexer/mm/samples."""

from __future__ import annotations

import glob
import json
import pickle
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.compute as pc
import pyarrow.parquet as pq
import torch
from torch.utils.data import Dataset

from .config import DEFAULT_CFG, mm_emb_dirname, resolve_root

USER_SCALAR_COLS = ["103", "104", "105", "109"]
USER_LIST_COLS = ["106", "107", "108", "110"]


def _as_int(x: Any, default: int = 0) -> int:
    try:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return default
        return int(x)
    except (TypeError, ValueError):
        return default


def _as_float(x: Any, default: float = 0.0) -> float:
    try:
        if x is None or (isinstance(x, float) and np.isnan(x)):
            return default
        return float(x)
    except (TypeError, ValueError):
        return default


def user_row_to_vec(row: Any, dim: int) -> np.ndarray:
    """Turn encrypted user columns into a fixed vector (scalars + list summaries)."""
    get = row.get if isinstance(row, dict) else lambda c, default=None: row[c] if c in row else default
    vals: List[float] = []
    for c in USER_SCALAR_COLS:
        vals.append(_as_float(get(c)))
    for c in USER_LIST_COLS:
        v = get(c)
        if v is None or (isinstance(v, float) and np.isnan(v)):
            vals.extend([0.0, 0.0, 0.0, 0.0])
            continue
        arr = np.asarray(v).reshape(-1)
        if arr.dtype == object:
            arr = np.asarray([_as_float(x) for x in arr.tolist()], dtype=np.float32)
        else:
            arr = arr.astype(np.float32, copy=False)
        vals.append(float(arr.size))
        vals.append(float(arr.mean()) if arr.size else 0.0)
        vals.append(float(arr[0]) if arr.size else 0.0)
        vals.append(float(arr[-1]) if arr.size else 0.0)
    vec = np.asarray(vals, dtype=np.float32)
    out = np.zeros(dim, dtype=np.float32)
    n = min(dim, vec.size)
    out[:n] = vec[:n]
    return out


def _expand_globs(root: Path, patterns: Sequence[str]) -> List[str]:
    files: List[str] = []
    for pat in patterns:
        p = pat
        if not Path(p).is_absolute():
            p = str(root / p)
        files.extend(sorted(glob.glob(p)))
    return files


def _read_parquet_limited(files: Sequence[str], max_rows: Optional[int] = None) -> pd.DataFrame:
    tables = []
    n = 0
    for path in files:
        pf = pq.ParquetFile(path)
        for rg in range(pf.num_row_groups):
            table = pf.read_row_group(rg)
            if max_rows is not None:
                remain = max_rows - n
                if remain <= 0:
                    break
                if table.num_rows > remain:
                    table = table.slice(0, remain)
            tables.append(table)
            n += table.num_rows
        if max_rows is not None and n >= max_rows:
            break
    if not tables:
        return pd.DataFrame()
    return pa.concat_tables(tables).to_pandas()


def _load_indexer_maps(path: Path) -> Dict[str, dict]:
    with open(path, "rb") as f:
        raw = pickle.load(f)
    if not isinstance(raw, dict):
        raise TypeError(f"indexer.pkl is {type(raw)}, expected dict with keys u/i")
    maps = {
        "u": raw.get("u", {}) or raw.get("oid_to_rid_user", {}),
        "i": raw.get("i", {}) or raw.get("oid_to_rid", {}),
    }
    del raw
    return maps


def _cid_to_oid(cid: Any) -> Optional[int]:
    if cid is None:
        return None
    if isinstance(cid, (int, np.integer)):
        return int(cid)
    text = str(cid).strip()
    if not text:
        return None
    try:
        return int(text)
    except ValueError:
        return None


def _events_from_row(seq_val: Any) -> List[dict]:
    if seq_val is None:
        return []
    if isinstance(seq_val, np.ndarray):
        seq_val = seq_val.tolist()
    if isinstance(seq_val, dict):
        return [seq_val]
    if not isinstance(seq_val, list):
        return []
    out = []
    for e in seq_val:
        if isinstance(e, dict):
            out.append(e)
    return out


class TencentGRDataset(Dataset):
    """Load a local TencentGR subset and return padded next-item examples.

    Each item:
      user_features (D_user,)
      history_items (T,)
      history_actions (T,)  # 0/1/2
      history_embs (T, emb_dim)
      history_mask (T,)
      target_item ()
      target_action ()
      target_emb (emb_dim,)
    """

    def __init__(self, cfg: Optional[Dict[str, Any]] = None, split: str = "train"):
        self.cfg = dict(DEFAULT_CFG)
        if cfg:
            self.cfg.update(cfg)
        self.split = split
        self.root = resolve_root(self.cfg)
        self.emb_dim = int(self.cfg["emb_dim"])
        self.max_seq_len = int(self.cfg["max_seq_len"])
        self.user_feat_dim = int(self.cfg["user_feat_dim"])
        self.pad_id = int(self.cfg.get("pad_id", 0))
        self.seq_df = pd.DataFrame()
        self.item_feat = pd.DataFrame()
        self.user_feat = pd.DataFrame()
        self.indexer: Dict[str, dict] = {"u": {}, "i": {}}
        self.mm_embeddings: Dict[int, np.ndarray] = {}
        self.samples: List[dict] = []
        self._user_vec: Dict[int, np.ndarray] = {}
        self._emb_mat = np.zeros((1, self.emb_dim), dtype=np.float32)
        self._rid_to_row: Dict[int, int] = {}

    @classmethod
    def from_raw(cls, cfg: Optional[Dict[str, Any]] = None, split: str = "train") -> "TencentGRDataset":
        ds = cls(cfg, split=split)
        ds._load_raw()
        return ds

    @classmethod
    def from_cache(cls, cache_dir: str | Path, cfg: Optional[Dict[str, Any]] = None, split: str = "train") -> "TencentGRDataset":
        ds = cls(cfg, split=split)
        ds.load_cache(cache_dir)
        ds._apply_split()
        return ds

    def _load_raw(self) -> None:
        max_users = self.cfg.get("max_users")
        seq_files = _expand_globs(self.root, self.cfg["seq_shards"])
        if not seq_files:
            raise FileNotFoundError(f"no seq shards under {self.root} matching {self.cfg['seq_shards']}")
        print(f"[tencentgr] seq files={len(seq_files)} max_users={max_users}")
        self.seq_df = _read_parquet_limited(seq_files, max_rows=max_users)
        print(f"[tencentgr] seq_df={self.seq_df.shape}")

        self.samples = self._build_samples(self.seq_df)
        self._apply_split(rebuild_from_all=True)
        keep_users = {int(s["user_id"]) for s in self.samples}
        keep_items = set()
        for s in self.samples:
            keep_items.update(int(i) for i in s["history_items"])
            keep_items.add(int(s["target_item"]))
        print(f"[tencentgr] split={self.split} samples={len(self.samples)} users={len(keep_users)} items={len(keep_items)}")

        item_files = _expand_globs(self.root, ["item_feat/*.parquet"])
        self.item_feat = self._load_table_filtered(item_files, "item_id", keep_items)
        user_files = _expand_globs(self.root, ["user_feat/*.parquet"])
        self.user_feat = self._load_table_filtered(user_files, "user_id", keep_users)
        print(f"[tencentgr] item_feat={self.item_feat.shape} user_feat={self.user_feat.shape}")

        indexer_path = self.root / "indexer.pkl"
        full_maps = _load_indexer_maps(indexer_path)
        item_map = full_maps.get("i", {})
        oid_to_rid = {}
        rid_to_oid = {}
        for oid, rid in item_map.items():
            rid_i = _as_int(rid, default=-1)
            if rid_i in keep_items:
                oid_i = _cid_to_oid(oid)
                if oid_i is None:
                    continue
                oid_to_rid[oid_i] = rid_i
                rid_to_oid[rid_i] = oid_i
        del full_maps, item_map
        self.indexer = {"u": {}, "i": oid_to_rid, "rid_to_oid": rid_to_oid}
        print(f"[tencentgr] indexer subset oids={len(oid_to_rid)}")

        self._index_user_vectors()
        self.mm_embeddings = self._load_mm_embeddings(oid_to_rid)
        self._rebuild_emb_matrix()
        print(f"[tencentgr] mm_embeddings={len(self.mm_embeddings)}")

    def _load_table_filtered(self, files: Sequence[str], id_col: str, keep_ids: Iterable[int]) -> pd.DataFrame:
        keep = set(int(x) for x in keep_ids)
        if not files or not keep:
            return pd.DataFrame()
        keep_arr = pa.array(list(keep), type=pa.int64())
        tables = []
        for path in files:
            table = pq.read_table(path)
            if id_col not in table.column_names:
                continue
            col = table[id_col]
            if not pa.types.is_int64(col.type):
                col = pc.cast(col, pa.int64())
                table = table.set_column(table.schema.get_field_index(id_col), id_col, col)
            filt = table.filter(pc.is_in(table[id_col], value_set=keep_arr))
            if filt.num_rows:
                tables.append(filt)
        if not tables:
            return pd.DataFrame()
        return pa.concat_tables(tables).to_pandas()

    def _index_user_vectors(self) -> None:
        self._user_vec = {}
        if self.user_feat.empty or "user_id" not in self.user_feat.columns:
            return
        for row in self.user_feat.to_dict(orient="records"):
            uid = _as_int(row.get("user_id"), default=-1)
            if uid < 0:
                continue
            self._user_vec[uid] = user_row_to_vec(row, self.user_feat_dim)

    def _load_mm_embeddings(self, oid_to_rid: Dict[int, int]) -> Dict[int, np.ndarray]:
        dims: Dict[str, int] = dict(self.cfg.get("mm_emb_dims") or {})
        rid_parts: Dict[int, Dict[str, np.ndarray]] = {}
        if not oid_to_rid:
            print("[tencentgr] no OID/RID overlap; mm embeddings empty")
            return {}
        needed_cids = pa.array([str(oid) for oid in oid_to_rid.keys()], type=pa.string())
        for cfg_name in self.cfg["mm_emb_configs"]:
            dim = int(dims.get(cfg_name, 0))
            folder = self.root / "mm_emb" / mm_emb_dirname(cfg_name)
            files = sorted(glob.glob(str(folder / "*.parquet")))
            print(f"[tencentgr] mm {cfg_name} files={len(files)} dir={folder}", flush=True)
            n_hit = 0
            for fi, path in enumerate(files, start=1):
                table = pq.read_table(path, columns=["anonymous_cid", "emb"])
                filt = table.filter(pc.is_in(table["anonymous_cid"], value_set=needed_cids))
                if filt.num_rows == 0:
                    if fi == 1 or fi % 10 == 0 or fi == len(files):
                        print(f"[tencentgr]   {cfg_name} {fi}/{len(files)} hits=0", flush=True)
                    continue
                pdf = filt.to_pandas()
                oids = pd.to_numeric(pdf["anonymous_cid"], errors="coerce")
                keep = oids.notna()
                if not bool(keep.any()):
                    continue
                oids_i = oids.loc[keep].astype(np.int64).to_numpy()
                embs = pdf.loc[keep, "emb"].to_numpy()
                for oid, emb in zip(oids_i, embs):
                    rid = oid_to_rid.get(int(oid))
                    if rid is None:
                        continue
                    vec = np.asarray(emb, dtype=np.float32).reshape(-1)
                    if dim and vec.size != dim:
                        fixed = np.zeros(dim, dtype=np.float32)
                        n = min(dim, vec.size)
                        fixed[:n] = vec[:n]
                        vec = fixed
                    rid_parts.setdefault(int(rid), {})[cfg_name] = vec
                    n_hit += 1
                if fi == 1 or fi % 10 == 0 or fi == len(files):
                    print(f"[tencentgr]   {cfg_name} {fi}/{len(files)} file_rows={filt.num_rows} total_hits={n_hit}", flush=True)
            print(f"[tencentgr] mm {cfg_name} matched_rows={n_hit}", flush=True)
        return self._concat_mm_parts(rid_parts, dims)

    def _concat_mm_parts(self, rid_parts: Dict[int, Dict[str, np.ndarray]], dims: Dict[str, int]) -> Dict[int, np.ndarray]:
        cfg_names = list(self.cfg["mm_emb_configs"])
        out: Dict[int, np.ndarray] = {}
        for rid, parts in rid_parts.items():
            chunks = []
            for name in cfg_names:
                dim = int(dims.get(name, 0))
                if name in parts:
                    chunks.append(parts[name])
                else:
                    chunks.append(np.zeros(dim or 0, dtype=np.float32))
            vec = np.concatenate(chunks, axis=-1).astype(np.float32, copy=False) if chunks else np.zeros(self.emb_dim, dtype=np.float32)
            if vec.size != self.emb_dim:
                fixed = np.zeros(self.emb_dim, dtype=np.float32)
                n = min(self.emb_dim, vec.size)
                fixed[:n] = vec[:n]
                vec = fixed
            out[int(rid)] = vec
        return out

    def _rebuild_emb_matrix(self) -> None:
        if not self.mm_embeddings:
            self._emb_mat = np.zeros((1, self.emb_dim), dtype=np.float32)
            self._rid_to_row = {}
            return
        rids = np.fromiter(self.mm_embeddings.keys(), dtype=np.int64)
        mat = np.zeros((len(rids) + 1, self.emb_dim), dtype=np.float32)
        rid_to_row = {}
        for i, rid in enumerate(rids, start=1):
            mat[i] = self.mm_embeddings[int(rid)]
            rid_to_row[int(rid)] = i
        self._emb_mat = mat
        self._rid_to_row = rid_to_row

    def _emb_of(self, rid: int) -> np.ndarray:
        row = self._rid_to_row.get(int(rid), 0)
        return self._emb_mat[row]

    def _build_samples(self, seq_df: pd.DataFrame) -> List[dict]:
        samples: List[dict] = []
        max_t = self.max_seq_len
        if seq_df.empty:
            return samples
        for row in seq_df.itertuples(index=False):
            user_id = _as_int(getattr(row, "user_id", 0))
            events = _events_from_row(getattr(row, "seq", None))
            if len(events) < 2:
                continue
            item_ids = [_as_int(e.get("item_id")) for e in events]
            actions = [_as_int(e.get("action_type")) for e in events]
            timestamps = [_as_int(e.get("timestamp")) for e in events]
            if len(item_ids) < 2:
                continue
            samples.append(
                {
                    "user_id": user_id,
                    "history_items": item_ids[:-1][-max_t:],
                    "history_actions": actions[:-1][-max_t:],
                    "target_item": item_ids[-1],
                    "target_action": actions[-1],
                    "last_timestamp": timestamps[-1],
                }
            )
        return samples

    def _is_train_user(self, user_id: int) -> bool:
        ratio = float(self.cfg.get("split_ratio", 0.8))
        return ((int(user_id) * 1_000_003) % 100) < int(round(ratio * 100))

    def _apply_split(self, rebuild_from_all: bool = False) -> None:
        if rebuild_from_all:
            all_samples = self.samples
        else:
            all_samples = self.samples
        if self.split == "all":
            return
        want_train = self.split != "test"
        self.samples = [s for s in all_samples if self._is_train_user(s["user_id"]) == want_train]

    def dump_cache(self, cache_dir: str | Path) -> Path:
        cache = Path(cache_dir)
        cache.mkdir(parents=True, exist_ok=True)
        seq_path = cache / "seq_df.parquet"
        item_path = cache / "item_feat.parquet"
        user_path = cache / "user_feat.parquet"
        self.seq_df.to_parquet(seq_path, index=False)
        self.item_feat.to_parquet(item_path, index=False)
        self.user_feat.to_parquet(user_path, index=False)
        with open(cache / "indexer.pkl", "wb") as f:
            pickle.dump(self.indexer, f, protocol=pickle.HIGHEST_PROTOCOL)
        rids = np.fromiter(self.mm_embeddings.keys(), dtype=np.int64) if self.mm_embeddings else np.zeros((0,), dtype=np.int64)
        if len(rids):
            embs = np.stack([self.mm_embeddings[int(r)] for r in rids], axis=0)
        else:
            embs = np.zeros((0, self.emb_dim), dtype=np.float32)
        np.savez_compressed(cache / "mm_embeddings.npz", rids=rids, embs=embs)
        with open(cache / "samples.pkl", "wb") as f:
            pickle.dump(self.samples, f, protocol=pickle.HIGHEST_PROTOCOL)
        meta = {
            "cfg": {k: (str(v) if isinstance(v, Path) else v) for k, v in self.cfg.items()},
            "split": self.split,
            "n_seq": int(len(self.seq_df)),
            "n_item_feat": int(len(self.item_feat)),
            "n_user_feat": int(len(self.user_feat)),
            "n_mm": int(len(self.mm_embeddings)),
            "n_samples": int(len(self.samples)),
            "emb_dim": self.emb_dim,
            "user_feat_dim": self.user_feat_dim,
        }
        (cache / "meta.json").write_text(json.dumps(meta, indent=2, default=str))
        print(f"[tencentgr] wrote cache -> {cache}")
        return cache

    def load_cache(self, cache_dir: str | Path) -> None:
        cache = Path(cache_dir)
        self.seq_df = pd.read_parquet(cache / "seq_df.parquet")
        self.item_feat = pd.read_parquet(cache / "item_feat.parquet")
        self.user_feat = pd.read_parquet(cache / "user_feat.parquet")
        with open(cache / "indexer.pkl", "rb") as f:
            self.indexer = pickle.load(f)
        blob = np.load(cache / "mm_embeddings.npz", allow_pickle=False)
        rids = blob["rids"]
        embs = blob["embs"]
        self.mm_embeddings = {int(r): embs[i].astype(np.float32, copy=False) for i, r in enumerate(rids)}
        with open(cache / "samples.pkl", "rb") as f:
            self.samples = pickle.load(f)
        meta_path = cache / "meta.json"
        if meta_path.exists():
            meta = json.loads(meta_path.read_text())
            self.cfg.update(meta.get("cfg") or {})
            self.emb_dim = int(meta.get("emb_dim", self.emb_dim))
            self.user_feat_dim = int(meta.get("user_feat_dim", self.user_feat_dim))
        self._index_user_vectors()
        self._rebuild_emb_matrix()
        print(
            f"[tencentgr] loaded cache {cache} "
            f"seq={len(self.seq_df)} item={len(self.item_feat)} user={len(self.user_feat)} "
            f"mm={len(self.mm_embeddings)} samples={len(self.samples)}"
        )

    def __len__(self) -> int:
        return len(self.samples)

    def __getitem__(self, idx: int) -> Dict[str, torch.Tensor]:
        s = self.samples[idx]
        user_vec = self._user_vec.get(int(s["user_id"]))
        if user_vec is None:
            user_vec = np.zeros(self.user_feat_dim, dtype=np.float32)
        hist_items = [int(i) for i in s["history_items"]]
        hist_actions = [int(a) for a in s["history_actions"]]
        t = min(len(hist_items), self.max_seq_len)
        hist_items = hist_items[-t:]
        hist_actions = hist_actions[-t:]
        items_pad = np.full(self.max_seq_len, self.pad_id, dtype=np.int64)
        acts_pad = np.zeros(self.max_seq_len, dtype=np.int64)
        mask = np.zeros(self.max_seq_len, dtype=np.float32)
        embs = np.zeros((self.max_seq_len, self.emb_dim), dtype=np.float32)
        if t:
            items_pad[-t:] = np.asarray(hist_items, dtype=np.int64)
            acts_pad[-t:] = np.asarray(hist_actions, dtype=np.int64)
            mask[-t:] = 1.0
            for i, iid in enumerate(hist_items):
                embs[-t + i] = self._emb_of(iid)
        target_item = int(s["target_item"])
        return {
            "user_features": torch.tensor(user_vec, dtype=torch.float32),
            "history_items": torch.tensor(items_pad, dtype=torch.long),
            "history_actions": torch.tensor(acts_pad, dtype=torch.long),
            "history_embs": torch.tensor(embs, dtype=torch.float32),
            "history_mask": torch.tensor(mask, dtype=torch.float32),
            "target_item": torch.tensor(target_item, dtype=torch.long),
            "target_action": torch.tensor(int(s["target_action"]), dtype=torch.long),
            "target_emb": torch.tensor(self._emb_of(target_item), dtype=torch.float32),
            "last_timestamp": torch.tensor(int(s.get("last_timestamp", 0)), dtype=torch.long),
            "user_id": torch.tensor(int(s["user_id"]), dtype=torch.long),
        }


def pooled_history_emb(ds: TencentGRDataset, sample: dict) -> np.ndarray:
    items = sample.get("history_items") or []
    if not items:
        return np.zeros(ds.emb_dim, dtype=np.float32)
    mat = np.stack([ds._emb_of(int(i)) for i in items], axis=0)
    return mat.mean(axis=0).astype(np.float32)


def design_matrix_from_dataset(
    ds: TencentGRDataset,
    *,
    max_n: Optional[int] = None,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """X = user ⊕ mean(history_emb) ⊕ target_emb; Y = click; W = late vs early timestamp."""
    n = len(ds.samples) if max_n is None else min(len(ds.samples), int(max_n))
    d_user = ds.user_feat_dim
    d = d_user + 2 * ds.emb_dim
    X = np.zeros((n, d), dtype=np.float32)
    Y = np.zeros(n, dtype=np.float32)
    ts = np.zeros(n, dtype=np.int64)
    user_ids = np.zeros(n, dtype=np.int64)
    for i in range(n):
        s = ds.samples[i]
        uid = int(s["user_id"])
        user_ids[i] = uid
        uvec = ds._user_vec.get(uid, np.zeros(d_user, dtype=np.float32))
        hist = pooled_history_emb(ds, s)
        tgt = ds._emb_of(int(s["target_item"]))
        X[i, :d_user] = uvec
        X[i, d_user : d_user + ds.emb_dim] = hist
        X[i, d_user + ds.emb_dim :] = tgt
        Y[i] = 1.0 if int(s["target_action"]) >= 1 else 0.0
        ts[i] = int(s.get("last_timestamp", 0))
    med = np.median(ts) if n else 0
    W = (ts >= med).astype(np.int64)
    return X, Y, W, user_ids
