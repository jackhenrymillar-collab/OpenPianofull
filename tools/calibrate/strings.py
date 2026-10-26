#!/usr/bin/env python3
"""Measure string behaviour on single notes: inharmonicity and phantom partials.

    python3 tools/calibrate/strings.py --salamander PATH --model renders/notes [--out report.md]

Per note it measures:
  B        inharmonicity, fitted to the partial frequencies: f_n = n f0 sqrt(1 + B n^2)
  phantom  level of the phantom partials at 2 f_n (twice a partial's frequency), relative to
           the normal partial f_2n just above them, averaged over the partials where the two
           are far enough apart to tell apart [dB], at soft and loud velocities. Phantoms come
           from the string stretching (longitudinal motion), so they grow with velocity.
  sumtone  the same for sum tones f_n + f_(n+1), relative to f_(2n+1)

Needs numpy, scipy and soundfile. Run it from tools/calibrate (it imports analyze.py).
"""

import argparse
import math
from pathlib import Path

import numpy as np

import analyze as an


def spectrum(x, rate, t0, t1, pad=4):
    seg = x[int(t0 * rate):int(t1 * rate)]
    n = len(seg)
    w = np.hanning(n)
    N = 1 << int(math.ceil(math.log2(n * pad)))
    spec = np.abs(np.fft.rfft(seg * w, N)) ** 2
    return spec, rate / N


def onset_trim(x):
    env = np.abs(x)
    onset = int(np.argmax(env > env.max() * 10 ** (-30 / 20)))
    return x[onset:]


def peak_near(spec, df, f, tol):
    lo, hi = int((f - tol) / df), int((f + tol) / df) + 1
    lo = max(lo, 1)
    if hi >= len(spec) or hi <= lo:
        return None, None
    i = lo + int(np.argmax(spec[lo:hi]))
    # parabolic refinement on log magnitude
    if 0 < i < len(spec) - 1:
        a, b, c = np.log(spec[i - 1:i + 2] + 1e-30)
        d = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0
        return (i + d) * df, spec[i]
    return i * df, spec[i]


def fit_b(x, rate, midi):
    """Grid search for the B and pitch that put the most energy on the predicted partials."""
    f0 = 440.0 * 2 ** ((midi - 69) / 12)
    t1 = min(len(x) / rate, 1.0 if f0 < 300 else 0.6)
    spec, df = spectrum(x, rate, 0.02, t1)
    L = 10 * np.log10(spec + 1e-30)
    fmax = min(8000.0, 0.45 * rate)
    nmax = int(min(40, fmax / f0))
    if nmax < 5:
        return None, nmax
    from scipy.ndimage import maximum_filter1d
    Lmax = maximum_filter1d(L, 5)
    best = None
    for B in np.geomspace(2e-5, 2e-2, 160):
        for fs in np.linspace(0.985, 1.025, 81):
            n = np.arange(1, nmax + 1)
            f = n * f0 * fs * np.sqrt(1 + B * n * n)
            f = f[f < fmax]
            if len(f) < 5:
                continue
            i = (f / df).astype(int)
            score = Lmax[i].mean()
            if best is None or score > best[0]:
                best = (score, B, fs)
    B, fs = best[1], best[2]
    # Refine on the actual peak positions: least squares on (f_n / n)^2 = C (1 + B n^2)
    for _ in range(3):
        pts = []
        for n in range(1, nmax + 1):
            fe = n * f0 * fs * math.sqrt(1 + B * n * n)
            if fe > fmax:
                break
            f, p = peak_near(spec, df, fe, max(4 * df, 0.004 * fe))
            if f and 10 * math.log10(p + 1e-30) > Lmax.max() - 60:
                pts.append((n, f))
        if len(pts) < 5:
            break
        n = np.array([q[0] for q in pts], float)
        f = np.array([q[1] for q in pts])
        a = np.polyfit(n * n, (f / n) ** 2, 1)
        B = max(1e-7, a[0] / a[1])
        fs = math.sqrt(a[1]) / f0
    return B, fs


def phantoms(x, rate, midi, B, fs=1.0):
    """Mean level of phantoms 2 f_n relative to f_2n, and of f_n + f_(n+1) relative to f_(2n+1).
    Partials are located on the actual spectrum (real strings stray a little from the B law);
    the phantom sits at the sum of the two measured partial frequencies."""
    f0 = 440.0 * 2 ** ((midi - 69) / 12) * fs
    spec, df = spectrum(x, rate, 0.02, min(len(x) / rate, 1.0))
    found = {}

    def part(n):
        if n not in found:
            fe = n * f0 * math.sqrt(1 + B * n * n)
            f, p = peak_near(spec, df, fe, max(3 * df, 0.003 * fe))
            found[n] = (f or fe, p or 0.0)
        return found[n]

    out = {"double": [], "sum": []}
    for n in range(2, 40):
        for kind, a, b, c in (("double", n, n, 2 * n), ("sum", n, n + 1, 2 * n + 1)):
            fn, pn = part(c)
            if fn > 6000 or pn <= 0:
                continue
            fp = part(a)[0] + part(b)[0]
            gap = fn - fp
            if gap < 8.0:
                continue
            ip, i_n = int(round(fp / df)), int(round(fn / df))
            pp = spec[ip - 2:ip + 3].max()
            mid = (ip + i_n) // 2
            floor = spec[mid - 1:mid + 2].max()
            # a phantom hidden below the floor between the two counts as the floor
            out[kind].append(10 * math.log10(max(pp, floor) / pn))
    return {k: (float(np.median(v)) if len(v) >= 3 else None) for k, v in out.items()}


def inter_partial(x, rate, midi, B, fs):
    """Energy between the normal partials (phantoms, longitudinal modes, noise) relative to
    the energy on them, 300 Hz to 5 kHz, over the first 300 ms [dB]."""
    f0 = 440.0 * 2 ** ((midi - 69) / 12) * fs
    spec, df = spectrum(x, rate, 0.01, min(len(x) / rate, 0.3))
    f = np.arange(len(spec)) * df
    band = (f > max(300.0, 2.5 * f0)) & (f < 5000.0)
    on = np.zeros(len(spec), bool)
    n = np.arange(1, 400)
    fp = n * f0 * np.sqrt(1 + B * n * n)
    width = max(4 * df, 6.0)
    for q in fp[fp < 5200]:
        on |= np.abs(f - q) < width
    return 10 * math.log10(spec[band & ~on].sum() / max(1e-30, spec[band & on].sum()))


def measure(path, midi):
    x, rate = an.load(path)
    x = onset_trim(x)
    B, fs = fit_b(x, rate, midi)
    ph = phantoms(x, rate, midi, B, fs) if B else {"double": None, "sum": None}
    ph["inter"] = inter_partial(x, rate, midi, B, fs) if B else None
    return {"B": B, **ph}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--salamander", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--notes", default=",".join(str(n) for n in range(21, 109, 3)))
    ap.add_argument("--out")
    args = ap.parse_args()
    notes = [int(n) for n in args.notes.split(",")]
    lines = ["| note | B real | B model | phantom 2f soft real/model | loud real/model | sum tone loud real/model | between partials soft real/model | loud real/model |",
             "|---|---|---|---|---|---|---|---|"]
    f = lambda v, d=1: "-" if v is None else f"{v:.{d}f}"
    for midi in notes:
        r = {}
        for vel in (40, 124):
            sp = Path(args.salamander) / "samples" / f"{an.name_of(midi)}v{an.VELOCITY_LAYER[vel]}.flac"
            mp = Path(args.model) / f"n{midi}_v{vel}.wav"
            if sp.exists() and mp.exists():
                r[vel] = (measure(sp, midi), measure(mp, midi))
        if 124 not in r:
            continue
        rl, ml = r[124]
        rs, ms = r.get(40, ({}, {}))
        bf = lambda v: "-" if v is None else f"{v * 1e4:.2f}e-4"
        lines.append(f"| {an.name_of(midi)} | {bf(rl['B'])} | {bf(ml['B'])} | {f(rs.get('double'))} / {f(ms.get('double'))} "
                     f"| {f(rl['double'])} / {f(ml['double'])} | {f(rl['sum'])} / {f(ml['sum'])} "
                     f"| {f(rs.get('inter'))} / {f(ms.get('inter'))} | {f(rl['inter'])} / {f(ml['inter'])} |")
    text = "\n".join(lines)
    print(text)
    if args.out:
        Path(args.out).write_text(text + "\n")


if __name__ == "__main__":
    main()
