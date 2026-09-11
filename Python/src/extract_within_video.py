#!/usr/bin/env python3
"""Sliding-window multimodal embeddings for MSR-VTT + Layer-1 FSDS attribution.

Each video is cut into N_WINDOWS temporal windows. Every window is encoded as
    video (ViViT/ViT, 768) ⊕ audio (CLAP, 512) ⊕ text (GPT-2, 768) = 2048-d
and the N_WINDOWS vectors are concatenated, giving
    (N_videos * N_windows, 2048)
which is the matrix fed to the modality-attribution RF / FSDS step.

Example:
    python Python/src/extract_within_video.py \\
        --n-videos 8 --n-windows 100 --video-encoder vivit
"""
from __future__ import annotations

import argparse
import json
import random
import subprocess
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import torchaudio
from PIL import Image
from sklearn.ensemble import RandomForestClassifier


VIDEO_DIM = 768
AUDIO_DIM = 512
TEXT_DIM = 768
CONCAT_DIM = VIDEO_DIM + AUDIO_DIM + TEXT_DIM


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser()
    p.add_argument("--zip-path", default="/workspace/data/msrvtt/msrvtt_drive.bin")
    p.add_argument("--video-root", default="/workspace/data/msrvtt/TrainValVideo")
    p.add_argument("--audio-root", default="/workspace/data/msrvtt/audio")
    p.add_argument("--caption-json", default="/workspace/data/hf_ann/msrvtt_train_7k.json")
    p.add_argument("--output-dir", default="/workspace/data/msrvtt/features_window")
    p.add_argument("--repo-out", default="/workspace/experiments/msrvtt")
    p.add_argument("--n-windows", type=int, default=100)
    p.add_argument(
        "--n-videos",
        type=int,
        default=8,
        help="Total videos to encode, keeping any per_video checkpoints and sampling the rest.",
    )
    p.add_argument("--t-frames", type=int, default=32)
    p.add_argument("--stride", type=int, default=2)
    p.add_argument("--frame-size", type=int, default=224)
    p.add_argument("--sample-rate", type=int, default=16000)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument(
        "--video-encoder",
        choices=("vivit", "vit_mean"),
        default="vivit",
        help="vivit = google/vivit-b-16x2-kinetics400; vit_mean = ViT-B/16 temporal mean (CPU fallback).",
    )
    p.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    p.add_argument("--video-batch-size", type=int, default=1)
    p.add_argument("--text-batch-size", type=int, default=16)
    p.add_argument("--n-bootstrap", type=int, default=10, help="Bootstrap replicates for FSDS mean/variance.")
    p.add_argument("--n-estimators", type=int, default=100)
    p.add_argument(
        "--fsds-only",
        action="store_true",
        help="Skip encoding; rerun Layer-1 FSDS from saved video/audio/text npy.",
    )
    return p.parse_args()


def load_captions(path: str) -> dict[str, list[str]]:
    df = pd.read_json(path)
    vid_col = "video_id" if "video_id" in df.columns else "image_id"
    cap_col = "caption" if "caption" in df.columns else "captions"
    out: dict[str, list[str]] = {}
    for _, row in df.iterrows():
        vid = str(row[vid_col]).replace(".mp4", "")
        cap = row[cap_col]
        if isinstance(cap, list):
            caps = [str(c) for c in cap if str(c).strip()]
        else:
            caps = [str(cap)] if str(cap).strip() else []
        out[vid] = caps
    return out


def ensure_video(video_id: str, video_root: Path, zip_path: Path) -> Path:
    dest = video_root / f"{video_id}.mp4"
    if dest.exists() and dest.stat().st_size > 1000:
        return dest
    video_root.mkdir(parents=True, exist_ok=True)
    member = f"TrainValVideo/{video_id}.mp4"
    if not zip_path.exists():
        raise FileNotFoundError(f"missing video {dest} and zip {zip_path}")
    with zipfile.ZipFile(zip_path) as zf:
        if member not in zf.namelist():
            raise FileNotFoundError(member)
        zf.extract(member, video_root.parent)
    return dest


def load_wav_mono(path: Path, sr_expect: int = 16000) -> torch.Tensor:
    import wave

    with wave.open(str(path), "r") as wf:
        sr = wf.getframerate()
        nch = wf.getnchannels()
        sw = wf.getsampwidth()
        raw = wf.readframes(wf.getnframes())
    if sw == 2:
        arr = np.frombuffer(raw, dtype=np.int16).astype(np.float32) / 32768.0
    elif sw == 1:
        arr = np.frombuffer(raw, dtype=np.uint8).astype(np.float32) / 128.0 - 1.0
    else:
        arr = np.frombuffer(raw, dtype=np.int32).astype(np.float32) / 2147483648.0
    if nch > 1:
        arr = arr.reshape(-1, nch).mean(axis=1)
    wav = torch.from_numpy(np.ascontiguousarray(arr))
    if sr != sr_expect and wav.numel() > 0:
        wav = torchaudio.functional.resample(wav, sr, sr_expect)
    return wav


def _write_silence_wav(path: Path, sr: int, n_seconds: float = 10.0) -> None:
    import wave

    n = max(sr, int(sr * n_seconds))
    with wave.open(str(path), "w") as wf:
        wf.setnchannels(1)
        wf.setsampwidth(2)
        wf.setframerate(sr)
        wf.writeframes(b"\x00\x00" * n)


def _video_duration_sec(video_path: Path, default: float = 10.0) -> float:
    try:
        import av

        c = av.open(str(video_path))
        if c.duration and c.duration > 0:
            return float(c.duration) / 1e6
    except Exception:
        pass
    return default


def ensure_audio(video_id: str, video_path: Path, audio_root: Path, sr: int) -> Path:
    audio_root.mkdir(parents=True, exist_ok=True)
    out = audio_root / f"{video_id}.wav"
    if out.exists() and out.stat().st_size > 1000:
        return out
    cmd = [
        "ffmpeg", "-i", str(video_path),
        "-vn", "-acodec", "pcm_s16le",
        "-ar", str(sr), "-ac", "1",
        str(out), "-y",
    ]
    r = subprocess.run(cmd, capture_output=True)
    if r.returncode != 0 or not out.exists() or out.stat().st_size < 1000:
        dur = _video_duration_sec(video_path)
        _write_silence_wav(out, sr, dur)
        print(f"  [audio] no audio stream in {video_id}; wrote {dur:.1f}s silence")
    return out


def compute_windows(n_frames: int, t: int, stride: int, n_windows: int) -> list[np.ndarray]:
    span = t * stride
    max_start = max(0, n_frames - span)
    starts = np.linspace(0, max_start, n_windows).astype(int)
    windows = []
    for s in starts:
        idx = s + np.arange(t) * stride
        idx = np.clip(idx, 0, max(0, n_frames - 1))
        windows.append(idx.astype(int))
    return windows


def load_frames_pil(video_path: Path) -> tuple[list[Image.Image], float]:
    try:
        import av
    except ImportError as exc:
        raise ImportError("PyAV (av) is required to decode videos") from exc

    container = av.open(str(video_path))
    stream = container.streams.video[0]
    fps = float(stream.average_rate) if stream.average_rate else 30.0
    frames: list[Image.Image] = []
    for frame in container.decode(video=0):
        frames.append(frame.to_image().convert("RGB"))
    if not frames:
        frames = [Image.new("RGB", (224, 224), 0)]
        fps = 30.0
    return frames, fps


def slice_audio(wav: torch.Tensor, sr: int, t_start: float, t_end: float) -> tuple[torch.Tensor, torch.Tensor]:
    start = max(0, int(t_start * sr))
    end = min(wav.numel(), int(t_end * sr))
    if end <= start:
        return torch.zeros(sr), torch.tensor(0.0)
    seg = wav[start:end]
    mask = torch.tensor(1.0 if float(seg.abs().max()) > 1e-4 else 0.0)
    return seg, mask


class ViViTEncoder(torch.nn.Module):
    def __init__(self, device: str):
        super().__init__()
        from transformers import VivitImageProcessor, VivitModel

        self.processor = VivitImageProcessor.from_pretrained("google/vivit-b-16x2-kinetics400")
        self.model = VivitModel.from_pretrained("google/vivit-b-16x2-kinetics400").eval()
        for p in self.model.parameters():
            p.requires_grad = False
        self.out_dim = int(self.model.config.hidden_size)
        self.device = device
        self.model.to(device)

    @torch.inference_mode()
    def encode_windows(self, frames: list[Image.Image], windows: list[np.ndarray], batch_size: int = 1) -> torch.Tensor:
        feats = []
        for i in range(0, len(windows), batch_size):
            batch_imgs = []
            for idx in windows[i : i + batch_size]:
                batch_imgs.append([frames[int(j)] for j in idx])
            inputs = self.processor(batch_imgs, return_tensors="pt")
            pixel = inputs["pixel_values"].to(self.device)
            out = self.model(pixel_values=pixel)
            pooled = out.last_hidden_state[:, 0]
            feats.append(pooled.cpu().float())
        return torch.cat(feats, dim=0)


class ViTMeanEncoder(torch.nn.Module):
    """CPU-friendly 768-d video encoder: ViT-B/16 CLS, mean over T frames."""

    def __init__(self, device: str):
        super().__init__()
        from transformers import ViTImageProcessor, ViTModel

        self.processor = ViTImageProcessor.from_pretrained("google/vit-base-patch16-224")
        self.model = ViTModel.from_pretrained("google/vit-base-patch16-224").eval()
        for p in self.model.parameters():
            p.requires_grad = False
        self.out_dim = int(self.model.config.hidden_size)
        self.device = device
        self.model.to(device)

    @torch.inference_mode()
    def encode_windows(self, frames: list[Image.Image], windows: list[np.ndarray], batch_size: int = 8) -> torch.Tensor:
        feats = []
        for idx in windows:
            imgs = [frames[int(j)] for j in idx]
            inputs = self.processor(imgs, return_tensors="pt")
            pixel = inputs["pixel_values"].to(self.device)
            out = self.model(pixel_values=pixel)
            pooled = out.last_hidden_state[:, 0].mean(0)
            feats.append(pooled.cpu().float())
        return torch.stack(feats, dim=0)


class CLAPEncoder(torch.nn.Module):
    def __init__(self, device: str, native_sr: int = 16000):
        super().__init__()
        from transformers import ClapModel, ClapProcessor

        self.processor = ClapProcessor.from_pretrained("laion/clap-htsat-unfused")
        self.model = ClapModel.from_pretrained("laion/clap-htsat-unfused").eval()
        for p in self.model.parameters():
            p.requires_grad = False
        self.out_dim = int(self.model.config.projection_dim)
        self.device = device
        self.native_sr = native_sr
        self.model.to(device)
        self.target_sr = 48000

    @torch.inference_mode()
    def encode_segments(self, segments: list[torch.Tensor]) -> torch.Tensor:
        feats = []
        for seg in segments:
            wav = seg
            if wav.dim() > 1:
                wav = wav.mean(0)
            wav = wav.cpu().float()
            if self.native_sr != self.target_sr:
                wav = torchaudio.functional.resample(wav, self.native_sr, self.target_sr)
            try:
                arr = wav.numpy()
                inputs = self.processor(
                    audio=arr,
                    sampling_rate=self.target_sr,
                    return_tensors="pt",
                    padding=True,
                )
                inputs = {k: v.to(self.device) if torch.is_tensor(v) else v for k, v in inputs.items()}
                out = self.model.get_audio_features(**inputs)
                if not isinstance(out, torch.Tensor):
                    out = (
                        out.pooler_output
                        if getattr(out, "pooler_output", None) is not None
                        else out.last_hidden_state.mean(1)
                    )
                vec = out if out.dim() == 1 else out.squeeze(0)
                feats.append(vec.cpu().float().reshape(-1)[:AUDIO_DIM])
            except Exception:
                feats.append(torch.zeros(AUDIO_DIM))
        stacked = torch.stack(feats, dim=0)
        if stacked.shape[-1] != AUDIO_DIM:
            stacked = F.pad(stacked, (0, max(0, AUDIO_DIM - stacked.shape[-1])))[:, :AUDIO_DIM]
        return stacked


class GPT2Encoder(torch.nn.Module):
    def __init__(self, device: str):
        super().__init__()
        from transformers import GPT2Model, GPT2Tokenizer

        self.tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
        self.model = GPT2Model.from_pretrained("gpt2").eval()
        for p in self.model.parameters():
            p.requires_grad = False
        self.out_dim = int(self.model.config.hidden_size)
        self.device = device
        self.model.to(device)

    @torch.inference_mode()
    def encode_texts(self, texts: list[str], batch_size: int = 16) -> torch.Tensor:
        feats = []
        clean = [t if t.strip() else " " for t in texts]
        for i in range(0, len(clean), batch_size):
            batch = clean[i : i + batch_size]
            toks = self.tokenizer(
                batch,
                padding=True,
                truncation=True,
                max_length=32,
                return_tensors="pt",
            ).to(self.device)
            out = self.model(**toks)
            hidden = out.last_hidden_state
            lengths = toks["attention_mask"].sum(1) - 1
            idx = torch.arange(hidden.shape[0], device=hidden.device)
            feats.append(hidden[idx, lengths].cpu().float())
        return torch.cat(feats, dim=0)


def modality_shares(vimp: np.ndarray) -> dict[str, float]:
    v = float(vimp[:VIDEO_DIM].sum())
    a = float(vimp[VIDEO_DIM : VIDEO_DIM + AUDIO_DIM].sum())
    t = float(vimp[VIDEO_DIM + AUDIO_DIM :].sum())
    total = v + a + t + 1e-12
    return {"video_share": v / total, "audio_share": a / total, "text_share": t / total}


def fit_rf_vimp(X: np.ndarray, y: np.ndarray, seed: int = 42, n_estimators: int = 100) -> np.ndarray:
    rf = RandomForestClassifier(
        n_estimators=n_estimators,
        max_features="sqrt",
        random_state=seed,
        n_jobs=-1,
    )
    rf.fit(X, y)
    return rf.feature_importances_


def _summarize_share_rows(rows: list[dict]) -> dict[str, float]:
    keys = ("video_share", "audio_share", "text_share")
    out: dict[str, float] = {"n_bootstrap": float(len(rows))}
    for k in keys:
        arr = np.asarray([r[k] for r in rows], dtype=np.float64)
        mean = float(arr.mean()) if arr.size else 0.0
        var = float(arr.var(ddof=1)) if arr.size > 1 else 0.0
        out[k] = mean
        out[f"{k}_mean"] = mean
        out[f"{k}_var"] = var
        out[f"{k}_std"] = float(np.sqrt(var))
    return out


def bootstrap_shares(
    X: np.ndarray,
    y: np.ndarray,
    n_bootstrap: int = 10,
    seed: int = 42,
    n_estimators: int = 100,
) -> dict[str, float]:
    """Row-bootstrap RF VIMP n_bootstrap times; return mean/variance of modality shares."""
    rng = np.random.RandomState(seed)
    n = int(len(y))
    rows: list[dict] = []
    b = 0
    attempts = 0
    while len(rows) < n_bootstrap and attempts < n_bootstrap * 5:
        attempts += 1
        idx = rng.randint(0, n, size=n)
        yb = y[idx]
        if len(np.unique(yb)) < 2:
            continue
        vimp = fit_rf_vimp(X[idx], yb, seed=seed + b, n_estimators=n_estimators)
        rows.append(modality_shares(vimp))
        b += 1
    if not rows:
        vimp = fit_rf_vimp(X, y, seed=seed, n_estimators=n_estimators)
        rows.append(modality_shares(vimp))
    return _summarize_share_rows(rows)


def _fmt_share(rec: dict, name: str) -> str:
    return (
        f"{name}={rec[f'{name}_share_mean']:.3f} "
        f"(var={rec[f'{name}_share_var']:.4g})"
    )


def run_fsds(
    video_mat: np.ndarray,
    audio_mat: np.ndarray,
    text_mat: np.ndarray,
    video_ids: list[str],
    out_dir: Path,
    n_bootstrap: int = 10,
    n_estimators: int = 100,
    seed: int = 42,
) -> dict:
    """Layer-1 multimodal attribution on concatenated 2048-d windows."""
    n_videos, n_windows, _ = video_mat.shape
    concat_3d = np.concatenate([video_mat, audio_mat, text_mat], axis=-1)

    within = []
    for vi, vid in enumerate(video_ids):
        X = concat_3d[vi]
        y = np.array([0] * (n_windows // 2) + [1] * (n_windows - n_windows // 2))
        if len(np.unique(y)) < 2:
            continue
        rec = {"video_id": vid, **bootstrap_shares(X, y, n_bootstrap, seed + vi, n_estimators)}
        within.append(rec)
        print(
            f"  [within] {vid}: {_fmt_share(rec, 'video')} "
            f"{_fmt_share(rec, 'audio')} {_fmt_share(rec, 'text')}"
        )

    early = concat_3d[:, : n_windows // 2].reshape(-1, CONCAT_DIM)
    late = concat_3d[:, n_windows // 2 :].reshape(-1, CONCAT_DIM)
    Xp = np.concatenate([early, late], axis=0)
    yp = np.array([0] * len(early) + [1] * len(late))
    pooled = bootstrap_shares(Xp, yp, n_bootstrap, seed + 101, n_estimators)
    print(
        f"  [pooled early vs late] {_fmt_share(pooled, 'video')} "
        f"{_fmt_share(pooled, 'audio')} {_fmt_share(pooled, 'text')}"
    )

    mixture = None
    if n_videos >= 2:
        mid = n_videos // 2
        a = concat_3d[:mid].reshape(-1, CONCAT_DIM)
        b = concat_3d[mid:].reshape(-1, CONCAT_DIM)
        Xm = np.concatenate([a, b], axis=0)
        ym = np.array([0] * len(a) + [1] * len(b))
        mixture = bootstrap_shares(Xm, ym, n_bootstrap, seed + 202, n_estimators)
        print(
            f"  [mixture videos 0..{mid-1} vs {mid}..{n_videos-1}] "
            f"{_fmt_share(mixture, 'video')} {_fmt_share(mixture, 'audio')} {_fmt_share(mixture, 'text')}"
        )

    avg = None
    if within:
        clip_means = [
            {
                "video_share": r["video_share_mean"],
                "audio_share": r["audio_share_mean"],
                "text_share": r["text_share_mean"],
            }
            for r in within
        ]
        avg = _summarize_share_rows(clip_means)
        # Plot/report bootstrap uncertainty as the mean of per-clip bootstrap variances,
        # not the across-clip spread of the means.
        for key in ("video_share", "audio_share", "text_share"):
            boot_vars = np.asarray([r[f"{key}_var"] for r in within], dtype=np.float64)
            avg[f"{key}_across_clip_var"] = avg[f"{key}_var"]
            avg[f"{key}_across_clip_std"] = avg[f"{key}_std"]
            avg[f"{key}_var"] = float(boot_vars.mean()) if boot_vars.size else 0.0
            avg[f"{key}_std"] = float(np.sqrt(avg[f"{key}_var"]))
        avg["n_bootstrap"] = float(n_bootstrap)

    payload = {
        "n_bootstrap": n_bootstrap,
        "n_estimators": n_estimators,
        "within_video": within,
        "within_video_average": avg,
        "pooled_temporal": pooled,
        "mixture_shift": mixture,
        "dims": {"video": VIDEO_DIM, "audio": AUDIO_DIM, "text": TEXT_DIM, "concat": CONCAT_DIM},
    }
    (out_dir / "within_video_fsds.json").write_text(json.dumps(payload, indent=2))
    return payload


def plot_shares(payload: dict, out_png: Path) -> None:
    import matplotlib.pyplot as plt

    rows = []
    if payload.get("within_video_average"):
        rows.append(("Within-video\n(avg early vs late)", payload["within_video_average"]))
    if payload.get("pooled_temporal"):
        rows.append(("Pooled temporal\n(all early vs late)", payload["pooled_temporal"]))
    if payload.get("mixture_shift"):
        rows.append(("User mixture\n(video group A vs B)", payload["mixture_shift"]))
    if not rows:
        return

    n_boot = int(payload.get("n_bootstrap", 10))
    labels = [r[0] for r in rows]

    def _mean_std(rec: dict, key: str) -> tuple[float, float]:
        mean = rec.get(f"{key}_mean", rec.get(key, 0.0))
        std = rec.get(f"{key}_std", float(np.sqrt(rec.get(f"{key}_var", 0.0))))
        return float(mean), float(std)

    video_m, video_s = zip(*[_mean_std(r[1], "video_share") for r in rows])
    audio_m, audio_s = zip(*[_mean_std(r[1], "audio_share") for r in rows])
    text_m, text_s = zip(*[_mean_std(r[1], "text_share") for r in rows])
    x = np.arange(len(labels))
    w = 0.25
    err = dict(capsize=3.0, ecolor="#1A2332", linewidth=0.8)

    fig, ax = plt.subplots(figsize=(9.2, 4.8), dpi=160)
    ax.bar(x - w, video_m, w, yerr=video_s, label="Video (768)", color="#E07A3D", **err)
    ax.bar(x, audio_m, w, yerr=audio_s, label="Audio (512)", color="#2C4A6E", **err)
    ax.bar(x + w, text_m, w, yerr=text_s, label="Text (768)", color="#2F6B4F", **err)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylim(0, 1.08)
    ax.set_ylabel("Modality share of RF VIMP")
    ax.set_title(f"MSR-VTT Layer-1 multimodal attribution  (n_bootstrap={n_boot}, mean ± sd)")
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(out_png, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def main() -> None:
    args = parse_args()
    torch.set_grad_enabled(False)
    if args.device == "cpu":
        torch.set_num_threads(max(1, min(4, torch.get_num_threads())))

    random.seed(args.seed)
    np.random.seed(args.seed)
    torch.manual_seed(args.seed)

    video_root = Path(args.video_root)
    audio_root = Path(args.audio_root)
    zip_path = Path(args.zip_path)
    out_dir = Path(args.output_dir)
    repo_out = Path(args.repo_out)
    out_dir.mkdir(parents=True, exist_ok=True)
    repo_out.mkdir(parents=True, exist_ok=True)

    if args.fsds_only:
        video_mat = np.load(out_dir / "video_feat_3d.npy")
        audio_mat = np.load(out_dir / "audio_feat_3d.npy")
        text_mat = np.load(out_dir / "text_feat_3d.npy")
        meta = json.loads((out_dir / "meta.json").read_text())
        picked = list(meta["video_ids"])
        print(f"[fsds-only] {video_mat.shape} videos={picked} n_bootstrap={args.n_bootstrap}")
        print("\nLayer-1 FSDS modality attribution")
        payload = run_fsds(
            video_mat, audio_mat, text_mat, picked, out_dir,
            n_bootstrap=args.n_bootstrap, n_estimators=args.n_estimators, seed=args.seed,
        )
        plot_png = repo_out / "msrvtt_modality_attribution.png"
        plot_shares(payload, plot_png)
        (repo_out / "within_video_fsds.json").write_text(json.dumps(payload, indent=2))
        meta["n_bootstrap"] = args.n_bootstrap
        meta["n_estimators"] = args.n_estimators
        (repo_out / "meta.json").write_text(json.dumps(meta, indent=2))
        (out_dir / "meta.json").write_text(json.dumps(meta, indent=2))
        print(f"wrote {plot_png}")
        return
    ckpt_dir = out_dir / "per_video"
    ckpt_dir.mkdir(parents=True, exist_ok=True)

    captions = load_captions(args.caption_json)
    available = sorted(
        {
            p.stem
            for p in video_root.glob("*.mp4")
        }
        | ({n.split("/")[-1].replace(".mp4", "") for n in zipfile.ZipFile(zip_path).namelist() if n.endswith(".mp4")} if zip_path.exists() else set())
    )
    available = [v for v in available if v in captions]
    if not available:
        raise RuntimeError("no videos with captions found")
    existing = [
        p.stem
        for p in sorted(ckpt_dir.glob("*.npz"))
        if p.stem in captions
    ]
    rng = random.Random(args.seed)
    remaining = [v for v in available if v not in set(existing)]
    need = max(0, args.n_videos - len(existing))
    extra = rng.sample(remaining, min(need, len(remaining))) if need and remaining else []
    picked = existing + extra
    picked = sorted(
        picked,
        key=lambda x: int(x.replace("video", "")) if x.replace("video", "").isdigit() else x,
    )
    print(
        f"[pick] n={len(picked)} resume={len(existing)} new={len(extra)} {picked}"
    )

    print(f"[encoders] video={args.video_encoder} device={args.device}")
    if args.video_encoder == "vivit":
        video_enc = ViViTEncoder(args.device)
    else:
        video_enc = ViTMeanEncoder(args.device)
    audio_enc = CLAPEncoder(args.device, native_sr=args.sample_rate)
    text_enc = GPT2Encoder(args.device)
    print(f"[dims] video={video_enc.out_dim} audio={audio_enc.out_dim} text={text_enc.out_dim}")

    all_video, all_audio, all_text, all_mask = [], [], [], []
    meta_videos = []

    for vi, vid in enumerate(picked):
        ckpt = ckpt_dir / f"{vid}.npz"
        print(f"\n[{vi+1}/{len(picked)}] {vid}")
        if ckpt.exists():
            z = np.load(ckpt, allow_pickle=True)
            vf = torch.from_numpy(z["video"])
            af = torch.from_numpy(z["audio"])
            tf = torch.from_numpy(z["text"])
            am = torch.from_numpy(z["mask"])
            raw_meta = z["meta"]
            if isinstance(raw_meta, np.ndarray):
                raw_meta = raw_meta.item()
            rec = json.loads(raw_meta) if isinstance(raw_meta, str) else dict(raw_meta)
            audio_dead = bool(torch.as_tensor(af).abs().max() < 1e-8)
            if audio_dead:
                print(f"  resume video/text {tuple(vf.shape)}; re-encoding audio")
                vpath = ensure_video(vid, video_root, zip_path)
                apath = ensure_audio(vid, vpath, audio_root, args.sample_rate)
                frames, fps = load_frames_pil(vpath)
                n_frames = len(frames)
                windows = compute_windows(n_frames, args.t_frames, args.stride, args.n_windows)
                wav = load_wav_mono(apath, args.sample_rate)
                segs, masks = [], []
                half_span = (args.t_frames * args.stride) / (2.0 * fps)
                for idx in windows:
                    center_t = float(idx[args.t_frames // 2]) / fps
                    t_start = max(0.0, center_t - half_span)
                    t_end = t_start + 2 * half_span
                    seg, mask = slice_audio(wav, args.sample_rate, t_start, t_end)
                    segs.append(seg)
                    masks.append(mask)
                af = audio_enc.encode_segments(segs)
                am = torch.stack(masks)
                rec["n_frames"] = n_frames
                rec["fps"] = fps
                np.savez_compressed(
                    ckpt,
                    video=vf.numpy(),
                    audio=af.numpy(),
                    text=tf.numpy(),
                    mask=am.numpy(),
                    meta=json.dumps(rec),
                )
                print(f"  audio_feat {tuple(af.shape)} mask_mean={float(am.mean()):.2f} norm={float(af.norm(dim=1).mean()):.3f}")
            else:
                print(f"  resume {tuple(vf.shape)}")
            all_video.append(vf)
            all_audio.append(af)
            all_text.append(tf)
            all_mask.append(am)
            meta_videos.append(rec)
            continue

        vpath = ensure_video(vid, video_root, zip_path)
        apath = ensure_audio(vid, vpath, audio_root, args.sample_rate)
        frames, fps = load_frames_pil(vpath)
        n_frames = len(frames)
        windows = compute_windows(n_frames, args.t_frames, args.stride, args.n_windows)
        print(f"  frames={n_frames} fps={fps:.2f} windows={len(windows)} span={args.t_frames * args.stride / fps:.2f}s")

        vf = video_enc.encode_windows(frames, windows, batch_size=args.video_batch_size)
        print(f"  video_feat {tuple(vf.shape)}")

        wav = load_wav_mono(apath, args.sample_rate)
        segs, masks = [], []
        half_span = (args.t_frames * args.stride) / (2.0 * fps)
        for idx in windows:
            center_t = float(idx[args.t_frames // 2]) / fps
            t_start = max(0.0, center_t - half_span)
            t_end = t_start + 2 * half_span
            seg, mask = slice_audio(wav, args.sample_rate, t_start, t_end)
            segs.append(seg)
            masks.append(mask)
        af = audio_enc.encode_segments(segs)
        am = torch.stack(masks)
        print(f"  audio_feat {tuple(af.shape)} mask_mean={float(am.mean()):.2f}")

        caps = captions.get(vid, [""]) or [""]
        texts = [caps[wi % len(caps)] for wi in range(args.n_windows)]
        tf = text_enc.encode_texts(texts, batch_size=args.text_batch_size)
        print(f"  text_feat {tuple(tf.shape)}")

        rec = {
            "video_id": vid,
            "n_windows": args.n_windows,
            "n_frames": n_frames,
            "fps": fps,
            "window_span_frames": args.t_frames * args.stride,
            "window_span_sec": args.t_frames * args.stride / fps,
        }
        np.savez_compressed(
            ckpt,
            video=vf.numpy(),
            audio=af.numpy(),
            text=tf.numpy(),
            mask=am.numpy(),
            meta=json.dumps(rec),
        )
        all_video.append(vf)
        all_audio.append(af)
        all_text.append(tf)
        all_mask.append(am)
        meta_videos.append(rec)

    video_mat = torch.stack(all_video, dim=0).numpy()
    audio_mat = torch.stack(all_audio, dim=0).numpy()
    text_mat = torch.stack(all_text, dim=0).numpy()
    mask_mat = torch.stack(all_mask, dim=0).numpy()
    n_videos, n_windows = video_mat.shape[0], video_mat.shape[1]
    video_flat = video_mat.reshape(n_videos * n_windows, -1)
    audio_flat = audio_mat.reshape(n_videos * n_windows, -1)
    text_flat = text_mat.reshape(n_videos * n_windows, -1)
    concat_flat = np.concatenate([video_flat, audio_flat, text_flat], axis=1)
    labels_flat = np.repeat(np.arange(n_videos), n_windows)
    window_flat = np.tile(np.arange(n_windows), n_videos)

    np.save(out_dir / "video_feat_3d.npy", video_mat)
    np.save(out_dir / "audio_feat_3d.npy", audio_mat)
    np.save(out_dir / "text_feat_3d.npy", text_mat)
    np.save(out_dir / "audio_mask_3d.npy", mask_mat)
    np.save(out_dir / "video_feat.npy", video_flat)
    np.save(out_dir / "audio_feat.npy", audio_flat)
    np.save(out_dir / "text_feat.npy", text_flat)
    np.save(out_dir / "concat_feat.npy", concat_flat)
    np.save(out_dir / "video_labels.npy", labels_flat)
    np.save(out_dir / "window_index.npy", window_flat)

    meta = {
        "n_videos": n_videos,
        "n_windows": n_windows,
        "video_ids": picked,
        "T": args.t_frames,
        "stride": args.stride,
        "video_encoder": args.video_encoder,
        "device": args.device,
        "dims": {"video": VIDEO_DIM, "audio": AUDIO_DIM, "text": TEXT_DIM, "concat": CONCAT_DIM},
        "shapes": {
            "video_feat_3d": list(video_mat.shape),
            "audio_feat_3d": list(audio_mat.shape),
            "text_feat_3d": list(text_mat.shape),
            "concat_feat": list(concat_flat.shape),
        },
        "videos": meta_videos,
    }
    (out_dir / "meta.json").write_text(json.dumps(meta, indent=2))
    print("\n" + "=" * 60)
    print(f"saved to {out_dir}")
    print(f"  video_feat_3d {video_mat.shape}")
    print(f"  audio_feat_3d {audio_mat.shape}")
    print(f"  text_feat_3d  {text_mat.shape}")
    print(f"  concat_feat   {concat_flat.shape}   (N_videos*N_windows, 2048)")
    print("=" * 60)

    print("\nLayer-1 FSDS modality attribution")
    payload = run_fsds(
        video_mat, audio_mat, text_mat, picked, out_dir,
        n_bootstrap=args.n_bootstrap, n_estimators=args.n_estimators, seed=args.seed,
    )
    plot_png = repo_out / "msrvtt_modality_attribution.png"
    plot_shares(payload, plot_png)
    (repo_out / "within_video_fsds.json").write_text(json.dumps(payload, indent=2))
    (repo_out / "meta.json").write_text(json.dumps(meta, indent=2))
    # compact features for the repo (float16)
    np.savez_compressed(
        repo_out / "concat_feat_f16.npz",
        concat=concat_flat.astype(np.float16),
        labels=labels_flat,
        windows=window_flat,
        video_ids=np.array(picked),
    )
    print(f"wrote {plot_png}")


if __name__ == "__main__":
    main()
