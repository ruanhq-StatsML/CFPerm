#!/usr/bin/env python3
"""Extract 16kHz mono wav audio from MSR-VTT TrainVal videos via ffmpeg."""
from __future__ import annotations

import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
VIDEO_DIR = ROOT / "data" / "raw" / "msrvtt" / "TrainValVideo"
AUDIO_DIR = ROOT / "data" / "raw" / "msrvtt" / "audio"
WORKERS = 8

AUDIO_DIR.mkdir(parents=True, exist_ok=True)


def one(f: Path) -> str:
    out = AUDIO_DIR / f"{f.stem}.wav"
    if out.exists() and out.stat().st_size > 1000:
        return "skip"
    cmd = [
        "ffmpeg",
        "-i",
        str(f),
        "-vn",
        "-acodec",
        "pcm_s16le",
        "-ar",
        "16000",
        "-ac",
        "1",
        str(out),
        "-y",
    ]
    r = subprocess.run(cmd, capture_output=True)
    return "ok" if r.returncode == 0 else "fail"


def main():
    videos = sorted(VIDEO_DIR.glob("*.mp4"))
    print(f"found {len(videos)} videos → {AUDIO_DIR}", flush=True)
    ok = skip = fail = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        futs = [ex.submit(one, f) for f in videos]
        for i, fu in enumerate(as_completed(futs), 1):
            s = fu.result()
            ok += s == "ok"
            skip += s == "skip"
            fail += s == "fail"
            if i % 100 == 0 or i == len(videos):
                print(
                    f"processed {i}/{len(videos)} ok={ok} skip={skip} fail={fail}",
                    flush=True,
                )
    print("done", flush=True)


if __name__ == "__main__":
    main()
