#!/usr/bin/env python3
"""MSR-VTT within-video sliding-window embeddings.

Each video → N_WINDOWS windows; each window concat (ViViT 768 + CLAP 512 + GPT2 768) = 2048.
Also runs within-video early/late FSDS modality mass.

  python3 scripts/extract_msrvtt_within_video.py
"""
from __future__ import annotations

import json
import random
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
import torchaudio
from PIL import Image
from sklearn.ensemble import RandomForestClassifier
from torchvision import transforms
from torchvision.transforms import InterpolationMode

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "msrvtt"
ANNOTATION_JSON = RAW / "captions" / "msrvtt_train_7k.json"
ANNOTATION_CSV = RAW / "captions" / "train_7k.csv"
VIDEO_ROOT = RAW / "TrainValVideo"
AUDIO_ROOT = RAW / "audio"
OUTPUT_DIR = ROOT / "data" / "msrvtt" / "features_window"
RESULTS = ROOT / "results" / "msrvtt_window"

T = 32
STRIDE = 2
FRAME_SIZE = 224
SAMPLE_RATE = 16000
TOTAL_FRAMES = 300
N_WINDOWS = 100  # per video
FPS = 30.0
SEED = 42
N_VIDEOS = 5  # start with 5; raise when GPU available
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)
AUDIO_ROOT.mkdir(parents=True, exist_ok=True)

random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)


def ensure_audio(video_id: str) -> Path:
    wav = AUDIO_ROOT / f"{video_id}.wav"
    if wav.exists() and wav.stat().st_size > 1000:
        return wav
    mp4 = VIDEO_ROOT / f"{video_id}.mp4"
    if not mp4.exists():
        return wav
    cmd = [
        "ffmpeg",
        "-i",
        str(mp4),
        "-vn",
        "-acodec",
        "pcm_s16le",
        "-ar",
        str(SAMPLE_RATE),
        "-ac",
        "1",
        str(wav),
        "-y",
    ]
    subprocess.run(cmd, capture_output=True, check=False)
    return wav


def compute_windows(
    total_frames=TOTAL_FRAMES, T=T, stride=STRIDE, n_windows=N_WINDOWS
):
    span = T * stride
    max_start = max(0, total_frames - span)
    starts = np.linspace(0, max_start, n_windows).astype(int)
    windows = []
    for s in starts:
        idx = s + np.arange(T) * stride
        idx = np.clip(idx, 0, total_frames - 1)
        windows.append(idx)
    return windows


frame_tf = transforms.Compose(
    [
        transforms.Resize(
            int(FRAME_SIZE * 256 / 224), interpolation=InterpolationMode.BICUBIC
        ),
        transforms.CenterCrop(FRAME_SIZE),
        transforms.ToTensor(),
        transforms.Normalize([0.485, 0.456, 0.406], [0.229, 0.224, 0.225]),
    ]
)


def load_video_window(video_id, indices):
    path = VIDEO_ROOT / f"{video_id}.mp4"
    try:
        import decord

        vr = decord.VideoReader(str(path))
        n = len(vr)
        idx = [int(min(max(i, 0), n - 1)) for i in indices.tolist()]
        frames = vr.get_batch(idx).asnumpy()
        imgs = [Image.fromarray(f) for f in frames]
    except Exception as e:
        print(f"  [warn] video read fail {video_id}: {e}", flush=True)
        imgs = [Image.new("RGB", (FRAME_SIZE, FRAME_SIZE), 0) for _ in indices]
    return torch.stack([frame_tf(im) for im in imgs])


def load_audio_window(video_id, t_start, t_end, sr=SAMPLE_RATE):
    path = AUDIO_ROOT / f"{video_id}.wav"
    if not path.exists():
        return torch.zeros(sr), torch.tensor(0.0)
    try:
        wav, orig_sr = torchaudio.load(str(path))
        if orig_sr != sr:
            wav = torchaudio.functional.resample(wav, orig_sr, sr)
        if wav.shape[0] > 1:
            wav = wav.mean(0, keepdim=True)
        wav = wav.squeeze(0)
        start = max(0, int(t_start * sr))
        end = min(wav.shape[0], int(t_end * sr))
        seg = wav[start:end]
        if seg.numel() == 0:
            return torch.zeros(sr), torch.tensor(0.0)
        mask = torch.tensor(1.0 if seg.abs().max() > 1e-4 else 0.0)
        return seg, mask
    except Exception:
        return torch.zeros(sr), torch.tensor(0.0)


class ViViTEncoder(torch.nn.Module):
    def __init__(self):
        super().__init__()
        from transformers import VivitModel

        self.model = VivitModel.from_pretrained(
            "google/vivit-b-16x2-kinetics400"
        ).eval()
        for p in self.model.parameters():
            p.requires_grad = False
        self.out_dim = self.model.config.hidden_size

    @torch.no_grad()
    def forward(self, video):
        out = self.model(pixel_values=video)
        return out.last_hidden_state[:, 0]


class CLAPEncoder(torch.nn.Module):
    def __init__(self):
        super().__init__()
        from transformers import ClapModel

        self.model = ClapModel.from_pretrained("laion/clap-htsat-unfused").eval()
        for p in self.model.parameters():
            p.requires_grad = False
        self.out_dim = self.model.config.projection_dim
        cfg = self.model.config.audio_config
        self.spec_size = cfg.spec_size
        self.mel_bins = cfg.num_mel_bins
        self.target_samples = 48000 * 10
        self.hop_length = self.target_samples // self.spec_size
        self.mel = torchaudio.transforms.MelSpectrogram(
            sample_rate=48000,
            n_fft=2048,
            hop_length=self.hop_length,
            n_mels=self.mel_bins,
            f_min=50,
            f_max=14000,
            power=2.0,
        )
        self.amp = torchaudio.transforms.AmplitudeToDB(top_db=80)

    @torch.no_grad()
    def forward(self, audio, sample_rate=16000):
        if sample_rate != 48000:
            audio = torchaudio.functional.resample(audio, sample_rate, 48000)
        if audio.dim() == 1:
            audio = audio.unsqueeze(0)
        if audio.shape[-1] > self.target_samples:
            audio = audio[..., : self.target_samples]
        elif audio.shape[-1] < self.target_samples:
            audio = F.pad(audio, (0, self.target_samples - audio.shape[-1]))
        mel = self.amp(self.mel(audio))
        mel = mel.transpose(1, 2).unsqueeze(1)
        if mel.shape[2] > self.spec_size:
            mel = mel[:, :, : self.spec_size, :]
        elif mel.shape[2] < self.spec_size:
            mel = F.pad(mel, (0, 0, 0, self.spec_size - mel.shape[2]))
        mel = (mel - mel.mean(-1, keepdim=True)) / (mel.std(-1, keepdim=True) + 1e-6)
        out = self.model.get_audio_features(
            input_features=mel,
            is_longer=torch.zeros(mel.shape[0], dtype=torch.bool, device=mel.device),
        )
        if isinstance(out, torch.Tensor):
            return out
        if hasattr(out, "pooler_output") and out.pooler_output is not None:
            return out.pooler_output
        return out.last_hidden_state.mean(1)


class GPT2Encoder(torch.nn.Module):
    def __init__(self):
        super().__init__()
        from transformers import GPT2Model, GPT2Tokenizer

        self.tokenizer = GPT2Tokenizer.from_pretrained("gpt2")
        if self.tokenizer.pad_token is None:
            self.tokenizer.pad_token = self.tokenizer.eos_token
        self.model = GPT2Model.from_pretrained("gpt2").eval()
        for p in self.model.parameters():
            p.requires_grad = False
        self.out_dim = self.model.config.hidden_size

    @torch.no_grad()
    def forward(self, texts):
        toks = self.tokenizer(
            texts,
            padding=True,
            truncation=True,
            max_length=32,
            return_tensors="pt",
        ).to(self.model.device)
        out = self.model(**toks)
        hidden = out.last_hidden_state
        lengths = toks["attention_mask"].sum(1) - 1
        idx = torch.arange(hidden.shape[0], device=hidden.device)
        return hidden[idx, lengths]


def main():
    print(f"DEVICE={DEVICE}", flush=True)
    if ANNOTATION_JSON.exists():
        rows = json.loads(ANNOTATION_JSON.read_text())
        video_captions = {}
        for r in rows:
            vid = str(r["video_id"]).replace(".mp4", "")
            cap = r.get("caption", [])
            if isinstance(cap, str):
                try:
                    import ast

                    cap = ast.literal_eval(cap)
                except Exception:
                    cap = [cap]
            video_captions[vid] = list(cap) if isinstance(cap, list) else [str(cap)]
    else:
        df = pd.read_csv(ANNOTATION_CSV)
        vid_col = "video_id" if "video_id" in df.columns else "image_id"
        cap_col = "caption" if "caption" in df.columns else "captions"
        video_captions = {}
        for vid, g in df.groupby(vid_col):
            caps = []
            for c in g[cap_col].tolist():
                if isinstance(c, str) and c.startswith("["):
                    import ast

                    try:
                        caps.extend(ast.literal_eval(c))
                    except Exception:
                        caps.append(c)
                else:
                    caps.append(str(c))
            video_captions[str(vid).replace(".mp4", "")] = caps

    # only videos that exist on disk
    available = [
        v
        for v in sorted(video_captions.keys())
        if (VIDEO_ROOT / f"{v}.mp4").exists()
    ]
    print(f"available videos with captions: {len(available)}", flush=True)
    n_take = min(N_VIDEOS, len(available))
    picked = random.sample(available, n_take)
    print(f"[pick] {picked}", flush=True)

    print("extract audio…", flush=True)
    for vid in picked:
        ensure_audio(vid)

    print("\n[load encoders]", flush=True)
    video_enc = ViViTEncoder().to(DEVICE)
    audio_enc = CLAPEncoder().to(DEVICE)
    text_enc = GPT2Encoder().to(DEVICE)
    print(
        f"[dims] video={video_enc.out_dim}, audio={audio_enc.out_dim}, "
        f"text={text_enc.out_dim}",
        flush=True,
    )

    windows = compute_windows()
    print(
        f"[window] N_WINDOWS={N_WINDOWS}, span={T * STRIDE} frames "
        f"= {T * STRIDE / FPS:.2f}s",
        flush=True,
    )

    all_video_feats, all_audio_feats, all_text_feats = [], [], []
    all_audio_masks, all_meta = [], []

    for vi, vid in enumerate(picked):
        print(f"\n[{vi+1}/{n_take}] {vid}", flush=True)
        caps = video_captions.get(vid, [])
        video_windows, audio_windows, text_windows, audio_mask_windows = [], [], [], []

        for wi, indices in enumerate(windows):
            video = load_video_window(vid, indices)
            video_t = video.unsqueeze(0).to(DEVICE)
            vf = video_enc(video_t).cpu().float().squeeze(0)
            video_windows.append(vf)

            center_frame = int(indices[T // 2])
            center_t = center_frame / FPS
            half_span = (T * STRIDE) / (2 * FPS)
            t_start = max(0.0, center_t - half_span)
            t_end = t_start + 2 * half_span
            wav, amask = load_audio_window(vid, t_start, t_end)
            af = audio_enc(wav.unsqueeze(0).to(DEVICE)).cpu().float().squeeze(0)
            audio_windows.append(af)
            audio_mask_windows.append(amask)

            text = caps[wi % len(caps)] if caps else ""
            tf = text_enc([text]).cpu().float().squeeze(0)
            text_windows.append(tf)

            if wi == 0 or (wi + 1) % 20 == 0:
                print(
                    f"  window {wi}: v={tuple(vf.shape)} a={tuple(af.shape)} "
                    f"t={tuple(tf.shape)} mask={float(amask)}",
                    flush=True,
                )

        all_video_feats.append(torch.stack(video_windows, 0))
        all_audio_feats.append(torch.stack(audio_windows, 0))
        all_text_feats.append(torch.stack(text_windows, 0))
        all_audio_masks.append(torch.stack(audio_mask_windows))
        all_meta.append(
            {
                "video_id": vid,
                "n_windows": N_WINDOWS,
                "window_span_frames": T * STRIDE,
                "window_span_sec": T * STRIDE / FPS,
            }
        )

    video_mat = torch.stack(all_video_feats, 0).numpy()
    audio_mat = torch.stack(all_audio_feats, 0).numpy()
    text_mat = torch.stack(all_text_feats, 0).numpy()
    audio_mask_mat = torch.stack(all_audio_masks, 0).numpy()

    N = video_mat.shape[0] * video_mat.shape[1]
    video_flat = video_mat.reshape(N, -1)
    audio_flat = audio_mat.reshape(N, -1)
    text_flat = text_mat.reshape(N, -1)
    labels_flat = np.repeat(np.arange(n_take), N_WINDOWS)

    np.save(OUTPUT_DIR / "video_feat_3d.npy", video_mat)
    np.save(OUTPUT_DIR / "audio_feat_3d.npy", audio_mat)
    np.save(OUTPUT_DIR / "text_feat_3d.npy", text_mat)
    np.save(OUTPUT_DIR / "audio_mask_3d.npy", audio_mask_mat)
    np.save(OUTPUT_DIR / "video_feat.npy", video_flat)
    np.save(OUTPUT_DIR / "audio_feat.npy", audio_flat)
    np.save(OUTPUT_DIR / "text_feat.npy", text_flat)
    np.save(OUTPUT_DIR / "video_labels.npy", labels_flat)

    # concat 2048
    concat = np.concatenate([video_flat, audio_flat, text_flat], 1)
    np.save(OUTPUT_DIR / "concat_feat.npy", concat)

    meta = {
        "n_videos": n_take,
        "n_windows": N_WINDOWS,
        "video_ids": picked,
        "window_span_frames": T * STRIDE,
        "window_span_sec": T * STRIDE / FPS,
        "T": T,
        "stride": STRIDE,
        "fps": FPS,
        "dims": {"video": 768, "audio": 512, "text": 768, "concat": 2048},
        "device": DEVICE,
        "per_video_meta": all_meta,
    }
    (OUTPUT_DIR / "meta.json").write_text(json.dumps(meta, indent=2))
    print(f"\nsaved features → {OUTPUT_DIR}", flush=True)
    print(
        f"  3d video/audio/text: {video_mat.shape} / {audio_mat.shape} / {text_mat.shape}",
        flush=True,
    )
    print(f"  flat concat: {concat.shape}", flush=True)

    # within-video FSDS
    print("\nWithin-video FSDS (early vs late windows)", flush=True)
    results = []
    for vi, vid in enumerate(picked):
        vf, af, tf = video_mat[vi], audio_mat[vi], text_mat[vi]
        Nw = vf.shape[0]
        if Nw < 4:
            continue
        labels = np.array([0] * (Nw // 2) + [1] * (Nw - Nw // 2))
        X = np.concatenate([vf, af, tf], 1)
        rf = RandomForestClassifier(
            n_estimators=100, max_features="sqrt", random_state=42, n_jobs=-1
        )
        rf.fit(X, labels)
        vimp = rf.feature_importances_
        v = float(vimp[:768].sum())
        a = float(vimp[768 : 768 + 512].sum())
        t = float(vimp[768 + 512 :].sum())
        total = v + a + t + 1e-12
        row = {
            "video_id": vid,
            "video_share": v / total,
            "audio_share": a / total,
            "text_share": t / total,
        }
        results.append(row)
        print(
            f"{vid}: video={row['video_share']:.3f}, audio={row['audio_share']:.3f}, "
            f"text={row['text_share']:.3f}",
            flush=True,
        )

    summary = {
        "n_videos": len(results),
        "n_windows": N_WINDOWS,
        "avg_video": float(np.mean([r["video_share"] for r in results])) if results else None,
        "avg_audio": float(np.mean([r["audio_share"] for r in results])) if results else None,
        "avg_text": float(np.mean([r["text_share"] for r in results])) if results else None,
        "per_video": results,
    }
    (OUTPUT_DIR / "within_video_fsds.json").write_text(json.dumps(summary, indent=2))
    (RESULTS / "within_video_fsds.json").write_text(json.dumps(summary, indent=2))

    md = [
        "# MSR-VTT within-video window embeddings\n\n",
        f"- videos: {n_take}, windows/video: {N_WINDOWS}\n",
        f"- dims: video 768 + audio 512 + text 768 = 2048\n",
        f"- device: {DEVICE}\n\n",
        "## Within-video FSDS mass (early vs late windows)\n\n",
        "| video | video | audio | text |\n|---|---:|---:|---:|\n",
    ]
    for r in results:
        md.append(
            f"| {r['video_id']} | {r['video_share']:.3f} | {r['audio_share']:.3f} | "
            f"{r['text_share']:.3f} |\n"
        )
    if results:
        md.append(
            f"| **avg** | {summary['avg_video']:.3f} | {summary['avg_audio']:.3f} | "
            f"{summary['avg_text']:.3f} |\n"
        )
    (RESULTS / "README.md").write_text("".join(md))
    (OUTPUT_DIR / "README.md").write_text("".join(md))
    print("done", summary, flush=True)


if __name__ == "__main__":
    main()
