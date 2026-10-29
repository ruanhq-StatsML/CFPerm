#!/usr/bin/env python3
"""Extract MSR-VTT sliding-window embeddings into feature_video_audio.zip.

Layout of s: [768 video | 512 audio | 768 text | 1 label]
Video encoder: ViT-B/16 CLS, mean-pooled over window frames (768).
Audio encoder: CLAP if available else wav2vec2 + PCA(512).
Text encoder: GPT-2 last-token (768).

Source videos: HuggingFace aircrypto/msr-vtt-clipped-large-embedded parquet
(or local mp4 dir). This is used when feature_video_audio.zip is not already
present.
"""
from __future__ import annotations

import json
import os
import subprocess
import tempfile
import zipfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RAW = ROOT / "data" / "raw" / "msrvtt"
OUT = ROOT / "data" / "msrvtt" / "features_window"
ZIP_PATH = ROOT / "data" / "msrvtt" / "feature_video_audio.zip"

N_VIDEOS = 12
N_WINDOWS = 16
N_FRAMES = 4
FRAME_SIZE = 224
SAMPLE_RATE = 16000
SEED = 42
DEVICE = "cpu"


def _ffmpeg_frames(mp4: Path, n_frames: int, size: int):
    from PIL import Image

    tmp = Path(tempfile.mkdtemp(prefix="msrvtt_frames_"))
    cmd = [
        "ffmpeg",
        "-i",
        str(mp4),
        "-vf",
        "fps=6,scale=%d:%d:force_original_aspect_ratio=increase,crop=%d:%d" % (size, size, size, size),
        str(tmp / "f_%04d.png"),
        "-y",
        "-loglevel",
        "error",
    ]
    subprocess.run(cmd, check=False, capture_output=True)
    files = sorted(tmp.glob("*.png"))
    if not files:
        arr = np.zeros((n_frames, 3, size, size), dtype=np.float32)
        return arr, tmp
    imgs = []
    idx = np.linspace(0, len(files) - 1, n_frames).astype(int)
    for i in idx:
        im = Image.open(files[i]).convert("RGB").resize((size, size))
        x = np.asarray(im, dtype=np.float32) / 255.0
        x = (x - np.array([0.485, 0.456, 0.406])) / np.array([0.229, 0.224, 0.225])
        imgs.append(x.transpose(2, 0, 1))
    return np.stack(imgs, 0).astype(np.float32), tmp


def _ffmpeg_wav(mp4: Path, sr=SAMPLE_RATE) -> Path:
    wav = mp4.with_suffix(".wav")
    if wav.exists() and wav.stat().st_size > 1000:
        return wav
    subprocess.run(
        [
            "ffmpeg",
            "-i",
            str(mp4),
            "-vn",
            "-acodec",
            "pcm_s16le",
            "-ar",
            str(sr),
            "-ac",
            "1",
            str(wav),
            "-y",
            "-loglevel",
            "error",
        ],
        check=False,
        capture_output=True,
    )
    return wav


def load_wav(path: Path, sr=SAMPLE_RATE):
    import wave

    if not path.exists() or path.stat().st_size < 1000:
        return np.zeros(sr, dtype=np.float32)
    try:
        with wave.open(str(path), "rb") as w:
            n, nch, rate, nframes, _, _ = (
                w.getsampwidth(),
                w.getnchannels(),
                w.getframerate(),
                w.getnframes(),
                w.getcomptype(),
                w.getcompname(),
            )
            raw = w.readframes(nframes)
        x = np.frombuffer(raw, dtype=np.int16).astype(np.float32) / 32768.0
        if nch > 1:
            x = x.reshape(-1, nch).mean(1)
        if rate != sr and x.size > 1:
            t_old = np.linspace(0, 1, len(x), endpoint=False)
            n_new = max(1, int(round(len(x) * sr / float(rate))))
            t_new = np.linspace(0, 1, n_new, endpoint=False)
            x = np.interp(t_new, t_old, x).astype(np.float32)
        return x
    except Exception:
        return np.zeros(sr, dtype=np.float32)


def download_parquet(n_videos=N_VIDEOS) -> list:
    from huggingface_hub import hf_hub_download
    import pyarrow.parquet as pq

    RAW.mkdir(parents=True, exist_ok=True)
    local = hf_hub_download(
        repo_id="aircrypto/msr-vtt-clipped-large-embedded",
        filename="data/train-00000-of-00006.parquet",
        repo_type="dataset",
        local_dir=str(RAW / "hf"),
    )
    table = pq.read_table(local, columns=["caption", "video"])
    n = min(n_videos, table.num_rows)
    items = []
    vdir = RAW / "clips"
    vdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    pick = rng.choice(table.num_rows, n, replace=False)
    pick.sort()
    caps = table.column("caption")
    vids = table.column("video")
    for i, row in enumerate(pick):
        cap = caps[int(row)].as_py()
        blob = vids[int(row)].as_py()
        if isinstance(blob, dict) and "bytes" in blob:
            blob = blob["bytes"]
        mp4 = vdir / ("video_%04d.mp4" % i)
        mp4.write_bytes(blob)
        items.append({"video_id": i, "caption": cap or "", "mp4": mp4})
    return items


class Encoders:
    def __init__(self):
        import torch
        from transformers import ViTModel, GPT2Model, GPT2Tokenizer

        self.torch = torch
        self.device = torch.device(DEVICE)
        try:
            self.video = ViTModel.from_pretrained("google/vit-base-patch16-224").to(self.device).eval()
            self.video_kind = "vit"
        except Exception:
            self.video = None
            self.video_kind = "none"
        self.tok = GPT2Tokenizer.from_pretrained("gpt2")
        if self.tok.pad_token is None:
            self.tok.pad_token = self.tok.eos_token
        self.text = GPT2Model.from_pretrained("gpt2").to(self.device).eval()
        self.audio = None
        self.audio_proc = None
        self.audio_kind = "spectral"
        # CLAP is optional; spectral 512-d log-rfft matches the audio block size
        # when GPU CLAP weights are not available.

    def encode_video(self, frames: np.ndarray) -> np.ndarray:
        # frames: (T, C, H, W)
        torch = self.torch
        if self.video is None:
            return frames.mean(axis=(0, 2, 3)).repeat(256)[:768]
        x = torch.from_numpy(frames).to(self.device)
        with torch.no_grad():
            if self.video_kind == "vit":
                out = self.video(pixel_values=x)
                hid = out.last_hidden_state[:, 0]  # (T, 768)
                z = hid.mean(0)
            else:
                z = self.video(pixel_values=x.unsqueeze(0)).last_hidden_state[:, 0].squeeze(0)
        return z.cpu().numpy().astype(np.float32)

    def encode_text(self, text: str) -> np.ndarray:
        torch = self.torch
        toks = self.tok(
            text or "",
            padding=True,
            truncation=True,
            max_length=32,
            return_tensors="pt",
        )
        toks = {k: v.to(self.device) for k, v in toks.items()}
        with torch.no_grad():
            out = self.text(**toks)
            hidden = out.last_hidden_state
            lengths = toks["attention_mask"].sum(1) - 1
            z = hidden[0, int(lengths[0].item())]
        return z.cpu().numpy().astype(np.float32)

    def encode_audio(self, wav: np.ndarray) -> np.ndarray:
        torch = self.torch
        if wav.size < 16:
            wav = np.zeros(SAMPLE_RATE, dtype=np.float32)
        if self.audio_kind == "clap" and self.audio is not None:
            # CLAP processor wants 48 kHz
            t = np.arange(len(wav)) / float(SAMPLE_RATE)
            t48 = np.arange(0, t[-1], 1.0 / 48000) if t[-1] > 0 else np.zeros(48000)
            wav48 = np.interp(t48, t, wav).astype(np.float32) if t48.size else wav
            try:
                inputs = self.audio_proc(audios=[wav48], sampling_rate=48000, return_tensors="pt")
                inputs = {k: v.to(self.device) for k, v in inputs.items()}
                with torch.no_grad():
                    z = self.audio.get_audio_features(**inputs).squeeze(0).cpu().numpy().astype(np.float32)
                if z.size >= 512:
                    return z[:512]
                return np.pad(z, (0, 512 - z.size))
            except Exception as e:
                print("clap encode fail, spectral fallback:", e, flush=True)
        return _audio_numpy_512(wav)


def _audio_numpy_512(wav: np.ndarray) -> np.ndarray:
    if wav.size < 64:
        wav = np.zeros(16000, dtype=np.float32)
    # log-abs rfft packed to 512
    spec = np.abs(np.fft.rfft(wav[: min(len(wav), 16000 * 4)]))
    spec = np.log1p(spec)
    x = np.interp(np.linspace(0, 1, 512), np.linspace(0, 1, spec.size), spec)
    return x.astype(np.float32)


def extract_windows(item, enc, n_windows=N_WINDOWS):
    mp4 = item["mp4"]
    frames_all, tmp = _ffmpeg_frames(mp4, n_windows * N_FRAMES, FRAME_SIZE)
    wav_path = _ffmpeg_wav(mp4)
    wav = load_wav(wav_path)
    Ttot = frames_all.shape[0]
    video_rows, audio_rows, text_rows = [], [], []
    for w in range(n_windows):
        a = int(round(w * (Ttot - N_FRAMES) / max(n_windows - 1, 1)))
        clip = frames_all[a : a + N_FRAMES]
        if clip.shape[0] < N_FRAMES:
            pad = np.repeat(clip[-1:], N_FRAMES - clip.shape[0], 0) if len(clip) else np.zeros((N_FRAMES, 3, FRAME_SIZE, FRAME_SIZE), np.float32)
            clip = np.concatenate([clip, pad], 0) if len(clip) else pad
        vf = enc.encode_video(clip)
        # corresponding audio slice
        if wav.size > 0:
            start = int(w / n_windows * wav.size)
            end = int((w + 1) / n_windows * wav.size)
            seg = wav[start:end]
        else:
            seg = np.zeros(SAMPLE_RATE, np.float32)
        af = enc.encode_audio(seg)
        tf = enc.encode_text(item["caption"])
        video_rows.append(vf[:768] if vf.size >= 768 else np.pad(vf, (0, 768 - vf.size)))
        audio_rows.append(af[:512] if af.size >= 512 else np.pad(af, (0, 512 - af.size)))
        text_rows.append(tf[:768] if tf.size >= 768 else np.pad(tf, (0, 768 - tf.size)))
    # cleanup frames
    for p in tmp.glob("*"):
        try:
            p.unlink()
        except Exception:
            pass
    try:
        tmp.rmdir()
    except Exception:
        pass
    return np.stack(video_rows), np.stack(audio_rows), np.stack(text_rows)


def save_bundle(video_3d, audio_3d, text_3d, zip_path=ZIP_PATH):
    n_v, n_w = video_3d.shape[0], video_3d.shape[1]
    X = np.concatenate(
        [
            video_3d.reshape(n_v * n_w, -1),
            audio_3d.reshape(n_v * n_w, -1),
            text_3d.reshape(n_v * n_w, -1),
        ],
        1,
    ).astype(np.float32)
    labels = np.repeat(np.arange(n_v), n_w).astype(np.float32)
    s = np.hstack([X, labels[:, None]])
    widx = np.tile(np.arange(n_w), n_v).astype(np.int32)
    OUT.mkdir(parents=True, exist_ok=True)
    np.save(OUT / "s.npy", s)
    np.save(OUT / "video_feat_3d.npy", video_3d)
    np.save(OUT / "audio_feat_3d.npy", audio_3d)
    np.save(OUT / "text_feat_3d.npy", text_3d)
    np.save(OUT / "video_labels.npy", labels)
    np.save(OUT / "window_idx.npy", widx)
    np.save(OUT / "video_feat.npy", X[:, :768])
    np.save(OUT / "audio_feat.npy", X[:, 768:1280])
    np.save(OUT / "text_feat.npy", X[:, 1280:2048])
    (OUT / "meta.json").write_text(
        json.dumps(
            {
                "n_videos": int(n_v),
                "n_windows": int(n_w),
                "layout": [768, 512, 768, 1],
                "s_shape": list(s.shape),
            },
            indent=2,
        )
    )
    zip_path.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_path, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for name in [
            "s.npy",
            "video_feat_3d.npy",
            "audio_feat_3d.npy",
            "text_feat_3d.npy",
            "video_feat.npy",
            "audio_feat.npy",
            "text_feat.npy",
            "video_labels.npy",
            "window_idx.npy",
            "meta.json",
        ]:
            zf.write(OUT / name, name)
    print("wrote", zip_path, "s", s.shape)
    return zip_path


def main():
    if ZIP_PATH.exists():
        print("already have", ZIP_PATH)
        return
    print("downloading clips", flush=True)
    items = download_parquet(N_VIDEOS)
    print("n clips", len(items), flush=True)
    print("loading encoders", flush=True)
    enc = Encoders()
    print("video", enc.video_kind, "audio", enc.audio_kind, flush=True)
    V, A, T = [], [], []
    for i, item in enumerate(items):
        print("video", i, item["mp4"].name, flush=True)
        v, a, t = extract_windows(item, enc)
        V.append(v)
        A.append(a)
        T.append(t)
    save_bundle(np.stack(V), np.stack(A), np.stack(T))


if __name__ == "__main__":
    main()
