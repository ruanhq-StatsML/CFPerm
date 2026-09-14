#!/usr/bin/env python3
"""Hop FSDS prototype: user — order — item — merchant.

原数据没有店。item 随机挂到 merchant；每笔转化发一个 order_number。
特征从 user/item 事件聚到店。归因不是图谱 GNN，是跳：
  item → merchant（上架）
  cnv  → order → item → merchant

FSDS = PO-risk + LOGO(跳) + LOCO(列)。φ=(Y-μ)(W-e)。R 不当检验。

  python3 scripts/tencent_gr/hop_fsds_proto.py
  python3 scripts/tencent_gr/hop_fsds_proto.py --names faker
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from po_fs_logo import po_risk_fit, rf_domain  # noqa: E402
from merchant_name import generate_catalog_names  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
EV = ROOT / "results/tencent_gr_fs150/tables/ev.parquet"
POST = ROOT / "results/tencent_gr_fs150/tables/post.parquet"
USER = ROOT / "results/tencent_gr_fs150/tables/user.parquet"
OUT = ROOT / "results/tencent_gr_fs150/hop"
SEED = 0
N_MERCH = 400
EXP, CLK, CNV = 0, 1, 2
HL7 = 7 * 86400.0
LN2 = math.log(2.0)
MIN_EXP = 3
SPLIT = 0.5

USER_HOP = [
    "life_ctr",
    "life_cnv_share",
    "dec_hl7d_dec_clk",
    "sess_bounce_rate",
    "post_clk_same_sess_rate",
]
ORDER_X = [
    "log1p_price",
    "empty_any",
    "sess_clk_before",
    "n_prior_cnv",
    "n_clk_before_1d",
]
MERCH_X = [
    "m_n_exp",
    "m_n_clk",
    "m_n_cnv",
    "m_n_order",
    "m_n_user",
    "m_n_item",
    "m_ctr",
    "m_cvr",
    "m_cnv_share",
    "m_dec_clk",
    "m_dec_cnv",
    "m_star",
    "m_open_days",
    "m_user_life_ctr",
    "m_user_life_cnv_share",
    "m_user_dec_clk",
    "m_user_bounce",
    "m_user_same_sess",
]


def hop_of(name: str) -> str:
    if name.startswith("m_") or name in MERCH_X:
        return "merchant"
    if name.startswith("u_") or name in USER_HOP:
        return "user"
    return "order"


def merch_family(name: str) -> str:
    if name in {"m_n_exp", "m_n_clk", "m_n_cnv", "m_n_order", "m_n_user", "m_n_item"}:
        return "volume"
    if name in {"m_ctr", "m_cvr", "m_cnv_share"}:
        return "structure"
    if name.startswith("m_dec"):
        return "recency"
    if name.startswith("m_user"):
        return "user_hop"
    return "shop_draw"


def shop_names(n: int, seed: int, *, backend: str = "faker") -> list[str]:
    if backend == "vllm":
        try:
            return _names_vllm(n, seed)
        except Exception as exc:
            print(f"vLLM names fallback ({exc})", flush=True)
            backend = "faker"
    return generate_catalog_names(n, seed, backend=backend)


def _names_vllm(n: int, seed: int) -> list[str]:
    from vllm import LLM, SamplingParams

    llm = LLM(model="facebook/opt-125m", max_model_len=64, seed=seed)
    sp = SamplingParams(temperature=1.1, max_tokens=8, seed=seed)
    prompts = [f"English shop name {i}:" for i in range(n)]
    outs = llm.generate(prompts, sp)
    names = []
    for o in outs:
        t = o.outputs[0].text.strip().split("\n")[0][:40] or "Nameless Mart"
        names.append(t)
    return names


def catalog(item_ids: np.ndarray, n_merch: int, seed: int, *, backend: str = "faker") -> pd.DataFrame:
    """上架在先：item 先挂到店，后面的点击/转化才能走到 merchant。"""
    item_ids = np.unique(np.asarray(item_ids, dtype=np.int64))
    mid = (item_ids * (10**9 + 7) + int(seed)) % int(n_merch)
    names = shop_names(n_merch, seed, backend=backend)
    rng = np.random.default_rng(seed + 3)
    shop = pd.DataFrame(
        {
            "merchant_id": np.arange(n_merch),
            "merchant_name": names,
            "m_star": rng.uniform(2.5, 5.0, n_merch),
            "m_open_days": rng.integers(30, 2400, n_merch).astype(float),
        }
    )
    imap = pd.DataFrame({"item_id": item_ids, "merchant_id": mid})
    return imap.merge(shop, on="merchant_id", how="left")


def mint_orders(cnv: pd.DataFrame, seed: int) -> pd.DataFrame:
    """转化当下才发单号。order 指向当时的 item，item 已经有 merchant。"""
    n = len(cnv)
    rng = np.random.default_rng(seed + 9)
    nums = rng.choice(np.arange(100_000_000, 200_000_000), size=n, replace=False)
    out = cnv.copy()
    out["order_number"] = nums.astype(str)
    return out


def _rate(n, d) -> pd.Series:
    d = d.replace(0, np.nan)
    return (n / d).fillna(0.0)


def merchant_x(left: pd.DataFrame, users: pd.DataFrame) -> pd.DataFrame:
    """左窗事件 → 店一行。量/结构/衰减 + 进店用户的均值（user→merchant 跳）。"""
    g = left.groupby("merchant_id", sort=False)
    is_exp, is_clk, is_cnv = left.act.eq(EXP), left.act.eq(CLK), left.act.eq(CNV)
    t_left = left.groupby("merchant_id")["ts"].transform("max")
    age = (t_left - left["ts"]).clip(lower=0)
    w7 = np.exp(-(LN2 / HL7) * age.to_numpy())
    dec_clk = pd.Series(w7, index=left.index)[is_clk].groupby(left.loc[is_clk, "merchant_id"]).sum()
    dec_cnv = pd.Series(w7, index=left.index)[is_cnv].groupby(left.loc[is_cnv, "merchant_id"]).sum()
    n_exp = is_exp.groupby(left.merchant_id).sum()
    n_clk = is_clk.groupby(left.merchant_id).sum()
    n_cnv = is_cnv.groupby(left.merchant_id).sum()
    n_tot = g.size()
    rec = pd.DataFrame(
        {
            "m_n_exp": n_exp,
            "m_n_clk": n_clk,
            "m_n_cnv": n_cnv,
            "m_n_order": n_cnv,
            "m_n_user": left.groupby("merchant_id")["user_id"].nunique(),
            "m_n_item": left.groupby("merchant_id")["item_id"].nunique(),
            "m_ctr": _rate(n_clk, n_exp),
            "m_cvr": _rate(n_cnv, n_clk),
            "m_cnv_share": _rate(n_cnv, n_tot),
            "m_dec_clk": dec_clk,
            "m_dec_cnv": dec_cnv,
        }
    ).fillna(0.0)
    touch = left[["merchant_id", "user_id"]].drop_duplicates()
    keep = [c for c in USER_HOP if c in users.columns]
    if keep:
        uh = touch.merge(users[["user_id"] + keep], on="user_id", how="left")
        um = uh.groupby("merchant_id")[keep].mean()
        um = um.rename(
            columns={
                "life_ctr": "m_user_life_ctr",
                "life_cnv_share": "m_user_life_cnv_share",
                "dec_hl7d_dec_clk": "m_user_dec_clk",
                "sess_bounce_rate": "m_user_bounce",
                "post_clk_same_sess_rate": "m_user_same_sess",
            }
        )
        rec = rec.join(um, how="left")
    return rec.fillna(0.0).reset_index()


def merchant_y(right: pd.DataFrame) -> pd.DataFrame:
    y = right.groupby("merchant_id").agg(
        y_n_exp=("act", lambda a: int((a == EXP).sum())),
        y_n_clk=("act", lambda a: int((a == CLK).sum())),
        y_n_cnv=("act", lambda a: int((a == CNV).sum())),
    )
    y["y_ctr"] = np.where(y.y_n_exp >= MIN_EXP, y.y_n_clk / y.y_n_exp, np.nan)
    y["y_any_cnv"] = (y.y_n_cnv > 0).astype(float)
    return y.reset_index()


def split_merchant(ev: pd.DataFrame, frac: float = SPLIT):
    mm = ev.groupby("merchant_id")["ts"].agg(t0="min", t_end="max").reset_index()
    mm["t_cut"] = mm["t0"] + (mm["t_end"] - mm["t0"]) * frac
    e = ev.merge(mm, on="merchant_id")
    left = e.loc[e.ts <= e.t_cut].drop(columns=["t0", "t_end", "t_cut"])
    right = e.loc[e.ts > e.t_cut]
    return left, right, mm


def _xy(df: pd.DataFrame, cols: list[str], ycol: str, clock: str, need: str | None):
    sub = df if need is None else df.loc[df[need].notna()]
    names = [c for c in cols if c in sub.columns]
    X = sub[names].replace([np.inf, -np.inf], np.nan).fillna(0.0).to_numpy(np.float64)
    y = pd.to_numeric(sub[ycol], errors="coerce").fillna(0.0).to_numpy(np.float64)
    clk = sub[clock].to_numpy(np.float64)
    w = (clk > np.median(clk)).astype(int)
    return sub, names, X, y, w


def board(X, y, w, names, *, hop_fn, seed: int, do_loco: bool) -> dict:
    families = sorted({hop_fn(n) for n in names})
    groups = {g: [i for i, n in enumerate(names) if hop_fn(n) == g] for g in families}
    print(f"  n={len(y)} p={len(names)} pos={float(y.mean()):.3f} W1={w.mean():.3f}", flush=True)
    auc, v_rf = rf_domain(X, w, seed=seed)
    print(f"  RF-domain AUC={auc:.3f}", flush=True)
    full = po_risk_fit(X, y, w, seed=seed)
    print(f"  PO-risk={full['risk']:.6f}", flush=True)
    logo = {}
    for g, ix in groups.items():
        keep = [j for j in range(X.shape[1]) if j not in set(ix)]
        if not keep:
            logo[g] = {"R_minus": float("nan"), "delta": float("nan")}
            continue
        rm = po_risk_fit(X[:, keep], y, w, seed=seed + 11 + families.index(g))["risk"]
        logo[g] = {"R_minus": rm, "delta": float(full["risk"] - rm)}
        print(f"  LOGO {g:10s} Δ={logo[g]['delta']:+.6f}", flush=True)
    pos = {g: max(logo[g]["delta"] or 0.0, 0.0) for g in families}
    s = sum(pos.values()) + 1e-12
    share = {g: pos[g] / s for g in families}
    po_v, rf_v = full["vimp"], v_rf
    ranked = [
        {
            "rank": r + 1,
            "name": names[j],
            "hop": hop_fn(names[j]),
            "po_vimp": float(po_v[j]),
            "rf_vimp": float(rf_v[j]),
        }
        for r, j in enumerate(np.argsort(-po_v))
    ]
    loco = []
    if do_loco:
        for j, n in enumerate(names):
            keep = [k for k in range(X.shape[1]) if k != j]
            rm = po_risk_fit(X[:, keep], y, w, seed=seed + 80 + j)["risk"]
            d = float(full["risk"] - rm)
            loco.append(
                {"name": n, "hop": hop_fn(n), "impurity": float(po_v[j]), "loco_dR": d}
            )
            print(f"  LOCO {n:28s} ΔR={d:+.6e}", flush=True)
        loco.sort(key=lambda r: -r["loco_dR"])
    mass_po, mass_rf = {}, {}
    for g in families:
        mass_po[g] = float(sum(po_v[i] for i in groups[g]))
        mass_rf[g] = float(sum(rf_v[i] for i in groups[g]))
    return {
        "n": int(len(y)),
        "p": int(len(names)),
        "pos_rate": float(y.mean()),
        "w1_rate": float(w.mean()),
        "rf_domain_auc": float(auc),
        "po_risk": float(full["risk"]),
        "family_n": {g: len(groups[g]) for g in families},
        "po_mass": {k: v / (sum(mass_po.values()) + 1e-12) for k, v in mass_po.items()},
        "rf_mass": {k: v / (sum(mass_rf.values()) + 1e-12) for k, v in mass_rf.items()},
        "logo": logo,
        "logo_share": share,
        "ranked_po": ranked,
        "loco": loco,
    }


def _md_board(title: str, rec: dict) -> list[str]:
    lines = [
        f"## {title}",
        "",
        f"n={rec['n']} p={rec['p']} pos={rec['pos_rate']:.3f} W1={rec['w1_rate']:.3f}  "
        f"RF-domain **{rec['rf_domain_auc']:.3f}**  PO-risk **{rec['po_risk']:.6f}**",
        "",
        "| hop | n | RF-domain | PO-VIMP | LOGO Δ | LOGO share |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    fam = rec["logo_share"]
    for g in sorted(fam, key=lambda x: -fam[x]):
        d = rec["logo"][g]["delta"]
        ds = f"{d:+.5f}" if d == d else "nan"
        lines.append(
            f"| {g} | {rec['family_n'][g]} | {rec['rf_mass'][g]:.3f} | "
            f"{rec['po_mass'][g]:.3f} | {ds} | {fam[g]:.3f} |"
        )
    lines += ["", "| rank | feat | hop | PO-VIMP | RF-VIMP |", "|---|---|---|---:|---:|"]
    for r in rec["ranked_po"][:16]:
        lines.append(
            f"| {r['rank']} | `{r['name']}` | {r['hop']} | {r['po_vimp']:.4f} | {r['rf_vimp']:.4f} |"
        )
    if rec["loco"]:
        lines += ["", "LOCO ΔR（>0 才是拿掉后 R 下降）：", "", "| feat | hop | impurity | LOCO ΔR |", "|---|---|---:|---:|"]
        for r in rec["loco"][:12]:
            lines.append(
                f"| `{r['name']}` | {r['hop']} | {r['impurity']:.4f} | {r['loco_dR']:+.6e} |"
            )
    lines.append("")
    return lines


def demo() -> None:
    t0 = 1_700_000_000
    rows = []
    for u in (1, 2, 3):
        for k in range(6):
            iid = 10 + (k % 4)
            rows.append((u, iid, EXP, t0 + u * 1000 + k * 400))
            if k % 2 == 0:
                rows.append((u, iid, CLK, t0 + u * 1000 + k * 400 + 20))
            if k == 4:
                rows.append((u, iid, CNV, t0 + u * 1000 + k * 400 + 80))
        rows.append((u, 99, EXP, t0 + 20 * 86400))
    ev = pd.DataFrame(rows, columns=["user_id", "item_id", "act", "ts"])
    cat = catalog(ev.item_id.to_numpy(), 3, 0, backend="faker")
    ev = ev.merge(cat[["item_id", "merchant_id"]], on="item_id")
    cnv = ev.loc[ev.act.eq(CNV), ["user_id", "item_id", "ts"]].rename(columns={"ts": "cnv_ts"})
    orders = mint_orders(cnv, 0)
    assert orders.order_number.nunique() == len(orders)
    assert set(ev.merchant_id.unique()) <= {0, 1, 2}
    print("demo orders", orders[["user_id", "item_id", "order_number"]].to_string(index=False))
    print("demo shops", cat[["item_id", "merchant_id", "merchant_name"]].drop_duplicates().head(8).to_string(index=False))


def main() -> None:
    import argparse

    ap = argparse.ArgumentParser()
    ap.add_argument("--n-merch", type=int, default=N_MERCH)
    ap.add_argument("--names", default="faker", choices=["faker", "template", "messy", "vllm"])
    ap.add_argument("--no-loco", action="store_true")
    args = ap.parse_args()
    demo()
    if not EV.exists():
        print("no ev.parquet, demo only")
        return

    import po_fs_logo

    po_fs_logo.TREES = 20
    po_fs_logo.DEPTH = 5

    ev = pd.read_parquet(EV)
    post = pd.read_parquet(POST)
    users = pd.read_parquet(USER)
    cat = catalog(ev.item_id.to_numpy(), args.n_merch, SEED, backend=args.names)
    ev = ev.merge(cat[["item_id", "merchant_id"]], on="item_id", how="left")
    cnv = ev.loc[ev.act.eq(CNV), ["user_id", "item_id", "merchant_id", "ts", "price"]].rename(
        columns={"ts": "cnv_ts"}
    )
    orders = mint_orders(cnv, SEED)
    orders = orders.merge(
        cat[["item_id", "merchant_name"]].drop_duplicates(),
        on="item_id",
        how="left",
    )
    keys = ["user_id", "item_id", "cnv_ts"]
    orders = orders.merge(post, on=keys, how="left", suffixes=("", "_p"))

    left, right, mm = split_merchant(ev)
    mx = merchant_x(left, users)
    shop = (
        cat[["merchant_id", "merchant_name", "m_star", "m_open_days"]]
        .drop_duplicates("merchant_id")
        .merge(mx, on="merchant_id", how="left")
        .merge(merchant_y(right), on="merchant_id", how="left")
        .merge(mm, on="merchant_id", how="left")
    )
    for c in ["y_n_exp", "y_n_clk", "y_n_cnv", "y_any_cnv"]:
        shop[c] = shop[c].fillna(0.0)

    OUT.mkdir(parents=True, exist_ok=True)
    cat.to_parquet(OUT / "item_merchant.parquet", index=False)
    orders.to_parquet(OUT / "orders.parquet", index=False)
    shop.to_parquet(OUT / "merchant.parquet", index=False)
    print(
        f"items={len(cat)} shops={shop.merchant_id.nunique()} orders={len(orders)} "
        f"names e.g. {shop.merchant_name.head(3).tolist()}",
        flush=True,
    )

    print("=== merchant grain (shop timeline 50% split) ===", flush=True)
    mcols = [c for c in MERCH_X if c in shop.columns]
    _, nm, Xm, ym, wm = _xy(shop, mcols, "y_ctr", "t_end", "y_ctr")
    rec_m = board(Xm, ym, wm, nm, hop_fn=merch_family, seed=SEED, do_loco=not args.no_loco)

    print("=== order grain (user / merchant / order hops) ===", flush=True)
    o = orders.merge(shop[["merchant_id"] + mcols], on="merchant_id", how="left")
    ukeep = [c for c in USER_HOP if c in users.columns]
    o = o.merge(users[["user_id"] + ukeep], on="user_id", how="left")
    xcols = ukeep + mcols + [c for c in ORDER_X if c in o.columns]
    o["_y"] = pd.to_numeric(o["y_post_clk_1d"], errors="coerce")
    _, no, Xo, yo, wo = _xy(o, xcols, "_y", "t_end", "_y")
    rec_o = board(Xo, yo, wo, no, hop_fn=hop_of, seed=SEED + 1, do_loco=not args.no_loco)

    payload = {
        "hops": "user — order — item — merchant",
        "n_merch": args.n_merch,
        "merchant": rec_m,
        "order": rec_o,
        "sample_shops": shop[["merchant_id", "merchant_name", "m_n_item", "m_n_cnv"]]
        .head(8)
        .to_dict(orient="records"),
    }
    (OUT / "HOP_FSDS.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    md = [
        "# Hop FSDS prototype（user / order / item / merchant）",
        "",
        "不是图谱 GNN。跳：`item → merchant`（上架），`cnv → order → item → merchant`。",
        "店名默认 `faker.company()`（标签，不进 X）。`--names template|messy|vllm`。",
        "PO-risk φ=(Y−μ)(W−e)；LOGO=整跳拿掉；LOCO=一列拿掉。R 不当检验。",
        "",
        "店粒：该店自己的时间 50% 切开，左 X 右 Y=`ctr`（右窗 n_clk/n_exp），W=`1{t_end>中位}`。",
        "订单粒：Y=满窗 `y_post_clk_1d`，X=用户历史 + 店聚合 + 当单 ctx。",
        "",
    ]
    md += _md_board("店粒", rec_m)
    md += _md_board("订单粒（三跳）", rec_o)
    md += [
        "样例店名： " + ", ".join(f"{r['merchant_name']}" for r in payload["sample_shops"][:4]),
        "",
        "`python3 scripts/tencent_gr/hop_fsds_proto.py`",
        "",
    ]
    (OUT / "HOP_FSDS.md").write_text("\n".join(md), encoding="utf-8")
    print("wrote", OUT / "HOP_FSDS.md")


if __name__ == "__main__":
    main()
