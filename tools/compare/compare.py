#!/usr/bin/env python3
"""Compare two folders of renders from openpiano-render.

    python3 tools/compare/compare.py renders/upstream renders/new --out report

For every scenario found in both folders this writes a figure with the two spectrograms and
their level over time, and a report.md with peak level, decay time and pitch error per file.

Needs numpy, scipy and matplotlib (pip install -r tools/compare/requirements.txt).
"""

import argparse
import math
import re
from pathlib import Path

import numpy as np
from scipy.io import wavfile
from scipy.signal import spectrogram

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NOTE_NAMES = {"C": 0, "D": 2, "E": 4, "F": 5, "G": 7, "A": 9, "B": 11}


def note_to_hz(name):
    """'C4' -> 261.63 Hz (equal temperament, A4 = 440 Hz)."""
    m = re.fullmatch(r"([A-G])(s?)(\d)", name)
    if not m:
        return None
    midi = 12 * (int(m.group(3)) + 1) + NOTE_NAMES[m.group(1)] + (1 if m.group(2) else 0)
    return 440.0 * 2 ** ((midi - 69) / 12)


def load(path):
    rate, data = wavfile.read(path)
    data = data.astype(np.float64)
    if data.ndim > 1:
        data = data[:, 0]
    return rate, data


def db(x):
    return 20 * np.log10(np.maximum(np.abs(x), 1e-12))


def envelope_db(x, rate, window_s=0.02):
    hop = max(1, int(rate * window_s))
    frames = len(x) // hop
    rms = np.sqrt(np.mean(x[: frames * hop].reshape(frames, hop) ** 2, axis=1))
    return np.arange(frames) * window_s, db(rms)


def decay_time(t, env):
    """Seconds for the level to fall 30 dB from its peak (T30), or None."""
    if len(env) == 0 or env.max() < -100:
        return None
    peak = int(np.argmax(env))
    below = np.where(env[peak:] < env[peak] - 30)[0]
    return None if len(below) == 0 else t[peak + below[0]] - t[peak]


def fundamental(x, rate, expected_hz):
    """Strongest spectral peak within four semitones of the expected pitch."""
    seg = x[int(0.05 * rate): int(1.05 * rate)]
    if len(seg) < rate // 4 or np.max(np.abs(seg)) == 0:
        return None
    n = 1 << int(math.ceil(math.log2(len(seg) * 8)))
    spec = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), n))
    freqs = np.fft.rfftfreq(n, 1 / rate)
    lo, hi = expected_hz * 2 ** (-4 / 12), expected_hz * 2 ** (4 / 12)
    band = (freqs >= lo) & (freqs <= hi)
    if not band.any() or spec[band].max() == 0:
        return None
    return freqs[band][np.argmax(spec[band])]


def analyse(path):
    rate, x = load(path)
    t, env = envelope_db(x, rate)
    info = {"rate": rate, "x": x, "t": t, "env": env, "peak_db": db(np.max(np.abs(x))) if len(x) else -240}
    info["t30"] = decay_time(t, env)
    m = re.match(r"single_([A-G]s?\d)_", path.stem)
    expected = note_to_hz(m.group(1)) if m else None
    f = fundamental(x, rate, expected) if expected else None
    info["cents"] = 1200 * math.log2(f / expected) if f and expected else None
    return info


def fmt(v, unit="", digits=1):
    return "–" if v is None or (isinstance(v, float) and not math.isfinite(v)) else f"{v:.{digits}f}{unit}"


def plot(name, a, b, label_a, label_b, out_path):
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), constrained_layout=True)
    for ax, info, label in ((axes[0], a, label_a), (axes[1], b, label_b)):
        f, t, s = spectrogram(info["x"], info["rate"], nperseg=4096, noverlap=3072)
        ax.pcolormesh(t, f, 10 * np.log10(s + 1e-20), shading="auto", vmin=-160, vmax=-40, cmap="magma")
        ax.set_ylim(0, 5000)
        ax.set_ylabel("Hz")
        ax.set_title(f"{label}: {name}")
    axes[2].plot(a["t"], a["env"], label=label_a)
    axes[2].plot(b["t"], b["env"], label=label_b, alpha=0.8)
    axes[2].set_ylim(-140, 20)
    axes[2].set_xlabel("seconds")
    axes[2].set_ylabel("RMS dBFS")
    axes[2].legend()
    axes[2].grid(alpha=0.3)
    fig.savefig(out_path, dpi=90)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a", type=Path, help="first folder (e.g. renders/upstream)")
    ap.add_argument("b", type=Path, help="second folder (e.g. renders/new)")
    ap.add_argument("--out", type=Path, default=Path("report"))
    args = ap.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    la, lb = args.a.name, args.b.name
    names = sorted({p.stem for p in args.a.glob("*.wav")} & {p.stem for p in args.b.glob("*.wav")})
    if not names:
        raise SystemExit("No scenarios in common between the two folders.")

    rows = [f"# Render comparison: {la} vs {lb}", "",
            "Peak is the loudest sample. T30 is how long the sound takes to fall 30 dB from its peak. "
            "Pitch is how far the strongest peak within four semitones of the note is from equal temperament, "
            "in cents. Low notes with a weak fundamental can report a neighbouring partial instead.", "",
            f"| scenario | peak {la} | peak {lb} | T30 {la} | T30 {lb} | pitch {la} | pitch {lb} | identical |",
            "|---|---|---|---|---|---|---|---|"]

    for name in names:
        a = analyse(args.a / f"{name}.wav")
        b = analyse(args.b / f"{name}.wav")
        same = len(a["x"]) == len(b["x"]) and np.array_equal(a["x"], b["x"])
        rows.append(f"| [{name}]({name}.png) | {fmt(a['peak_db'], ' dB')} | {fmt(b['peak_db'], ' dB')} "
                    f"| {fmt(a['t30'], ' s', 2)} | {fmt(b['t30'], ' s', 2)} "
                    f"| {fmt(a['cents'], ' c')} | {fmt(b['cents'], ' c')} | {'yes' if same else 'no'} |")
        plot(name, a, b, la, lb, args.out / f"{name}.png")
        print(f"{name}: done")

    (args.out / "report.md").write_text("\n".join(rows) + "\n")
    print(f"Wrote {args.out / 'report.md'}")


if __name__ == "__main__":
    main()
