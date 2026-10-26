#!/usr/bin/env python3
"""Plot the harmonic level over time, real against model, from analyze.py's metrics.json.

    python3 tools/calibrate/plot.py report/calibrate/metrics.json --vel 92 --out report/calibrate/decay_v92.png
    python3 tools/calibrate/plot.py before.json after.json --vel 92 --out compare.png   (real, before, after)
"""

import argparse
import json

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NAMES = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("metrics", nargs="+")
    ap.add_argument("--vel", type=int, default=92)
    ap.add_argument("--out", required=True)
    ap.add_argument("--seconds", type=float, default=8)
    args = ap.parse_args()
    sets = [json.load(open(p)) for p in args.metrics]
    rows = [r for r in sets[0] if r["vel"] == args.vel]
    cols = 4
    nrows = (len(rows) + cols - 1) // cols
    fig, axs = plt.subplots(nrows, cols, figsize=(4 * cols, 2.6 * nrows), squeeze=False)
    for i, r in enumerate(rows):
        ax = axs[i // cols][i % cols]
        curves = [("real", r["real"])]
        for j, s in enumerate(sets):
            m = next(q for q in s if q["midi"] == r["midi"] and q["vel"] == r["vel"])
            curves.append(("model" if len(sets) == 1 else ["before", "after"][min(j, 1)], m["model"]))
        for label, d in curves:
            env = np.array(d.get("env", []))
            if env.size:
                ax.plot(np.arange(env.size) * 0.05, env - env.max(), label=label, lw=1)
        ax.set_title(f"{NAMES[r['midi'] % 12]}{r['midi'] // 12 - 1}  v{r['vel']}", fontsize=9)
        ax.set_xlim(0, args.seconds)
        ax.set_ylim(-80, 3)
        ax.grid(True, alpha=0.3)
        if i == 0:
            ax.legend(fontsize=7)
    for k in range(len(rows), nrows * cols):
        axs[k // cols][k % cols].axis("off")
    fig.suptitle("Harmonic level (first 20 partials) relative to its peak [dB] over time [s]")
    fig.tight_layout()
    fig.savefig(args.out, dpi=70)


if __name__ == "__main__":
    main()
