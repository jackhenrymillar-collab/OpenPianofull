#!/usr/bin/env python3
"""Measure single piano notes: decay rates and brightness, note by note.

Compares recorded single-note samples (the Salamander Grand Piano, FLAC, from
https://freepats.zenvoid.org/Piano/acoustic-grand-piano.html) with model renders from
openpiano-notes, at matching notes and MIDI velocities.

    build/tools/openpiano-notes --out renders/notes --seconds 8
    python3 tools/calibrate/analyze.py --salamander PATH/TO/SalamanderGrandPiano-SFZ+FLAC-V3+20200602 \
        --model renders/notes --out report/calibrate

Per note and velocity it measures, after the onset:
  early      decay of the harmonic level (first 20 partials, summed) from 0.1 to 1.0 s [dB/s]
  late       the same from 2 to 6 s (shorter for high notes) [dB/s]
  f1e, f1l   the same two decay rates for the first partial alone
  centroid   spectral centroid over 20-200 ms, in multiples of the fundamental
  hi         level above 2 kHz relative to the whole, over 20-200 ms [dB]
  tilt       level of partials 2-8 relative to the first, over 20-200 ms [dB]
  level      loudness over the first 200 ms [dB]; compared relative to C4 at velocity 92

Needs numpy, scipy and soundfile (pip install numpy scipy soundfile).
"""

import argparse
import json
import math
import re
from pathlib import Path

import numpy as np
import soundfile as sf
from scipy.signal import stft

NOTE_NAMES = {"C": 0, "D": 2, "E": 4, "F": 5, "G": 7, "A": 9, "B": 11}
# Salamander velocity layers -> the MIDI velocity in the middle of each layer's range
LAYER_VELOCITY = {1: 13, 2: 30, 3: 35, 4: 40, 5: 45, 6: 48, 7: 53, 8: 60, 9: 68, 10: 76,
                  11: 84, 12: 92, 13: 100, 14: 108, 15: 116, 16: 124}
VELOCITY_LAYER = {v: k for k, v in LAYER_VELOCITY.items()}


def midi_of(name):
    m = re.fullmatch(r"([A-G])(#?)(\d)", name)
    return 12 * (int(m.group(3)) + 1) + NOTE_NAMES[m.group(1)] + (1 if m.group(2) else 0)


def name_of(midi):
    names = ["C", "C#", "D", "D#", "E", "F", "F#", "G", "G#", "A", "A#", "B"]
    return f"{names[midi % 12]}{midi // 12 - 1}"


def load(path):
    data, rate = sf.read(str(path), dtype="float64", always_2d=True)
    return data.mean(axis=1), rate


def level_db(x):
    return 10 * np.log10(np.maximum(x, 1e-20))


def fit_slope(t, y, floor):
    """Least-squares slope in dB/s of points above the noise floor; None if too few."""
    ok = y > floor
    if ok.sum() < 5:
        return None
    a = np.polyfit(t[ok], y[ok], 1)
    return -a[0]


def analyze(x, rate, midi):
    f0 = 440.0 * 2 ** ((midi - 69) / 12)
    # Onset: the first moment within 30 dB of the peak (samples have a little lead-in)
    env = np.abs(x)
    peak = env.max()
    onset = int(np.argmax(env > peak * 10 ** (-30 / 20)))
    x = x[onset:]
    dur = len(x) / rate

    # Overall level, 10 ms hop, 40 ms RMS window
    hop, win = int(0.01 * rate), int(0.04 * rate)
    frames = max(1, (len(x) - win) // hop)
    rms = np.array([np.mean(x[i * hop:i * hop + win] ** 2) for i in range(frames)])
    t = np.arange(frames) * hop / rate + win / (2 * rate)
    L = level_db(rms)
    floor = np.percentile(L[-max(5, frames // 20):], 50) + 10 if dur > 3 else L.max() - 80
    floor = max(floor, L.max() - 75)

    def window_slope(t0, t1):
        sel = (t >= t0) & (t <= min(t1, dur))
        if sel.sum() < 10:
            return None
        return fit_slope(t[sel], L[sel], floor)

    # Late window: 2 to 6 s, brought earlier for high notes whose recordings are short
    late_end = min(6.0, dur - 0.3)
    late_start = min(2.0, max(1.0, late_end - 2.0))
    out = {"bb_early": window_slope(0.1, 1.0), "bb_late": window_slope(late_start, late_end)}

    # Spectrum: STFT with a window long enough to separate partials in the bass
    nfft = int(2 ** np.ceil(np.log2(max(2048, 6 * rate / f0))))
    nfft = min(nfft, 16384)
    fr, ft, Z = stft(x, fs=rate, nperseg=nfft, noverlap=nfft - int(0.01 * rate), boundary=None, padded=False)
    P = np.abs(Z) ** 2
    df = fr[1] - fr[0]

    # Partials: find each near its expected place, refining the inharmonicity as we go
    early_sel = (ft >= 0.02 + nfft / (2 * rate)) & (ft <= 0.2 + nfft / (2 * rate))
    if early_sel.sum() == 0:
        early_sel = np.zeros_like(ft, dtype=bool); early_sel[:3] = True
    spec_early = P[:, early_sel].mean(axis=1)
    found = []
    B = 0.0
    f1 = f0
    for n in range(1, 41):
        guess = n * f1 * math.sqrt((1 + B * n * n) / (1 + B))
        if guess > min(rate / 2 - 500, 16000):
            break
        lo, hi = guess * (1 - 0.03 / max(1, n ** 0.3)) - 2 * df, guess * (1 + 0.03 / max(1, n ** 0.3)) + 2 * df
        lo = max(lo, guess - 0.4 * f1); hi = min(hi, guess + 0.4 * f1)
        band = (fr >= lo) & (fr <= hi)
        if not band.any():
            continue
        k = np.argmax(np.where(band, spec_early, 0))
        fk = fr[k]
        if n == 1:
            f1 = fk
        found.append((n, fk, k))
        if len(found) >= 4:
            ns = np.array([q[0] for q in found], float)
            fs = np.array([q[1] for q in found], float)
            # (f_n / (n f1))^2 = (1 + B n^2) / (1 + B) ~ 1 + B (n^2 - 1)
            yv = (fs / (ns * f1)) ** 2 - 1
            xv = ns ** 2 - 1
            if (xv > 0).any():
                B = max(0.0, float(np.sum(xv * yv) / np.sum(xv * xv)))

    # Refine B on a long, zero-padded window (0.05-1.5 s) with interpolated peaks, fitting
    # f_n^2 / n^2 = F^2 (1 + B n^2) without trusting the (often weak) fundamental in the bass
    seg = x[int(0.05 * rate):int(min(1.5, dur) * rate)]
    if len(seg) > rate // 4 and len(found) >= 6:
        nz = 1 << 18
        S = np.abs(np.fft.rfft(seg * np.hanning(len(seg)), nz))
        fz = np.fft.rfftfreq(nz, 1 / rate)
        dz = fz[1]
        ns_, fs_ = [], []
        for n, fk, _ in found:
            if n > 25:
                break
            lo_, hi_ = fk - 2.5 * df, fk + 2.5 * df
            band = np.where((fz >= lo_) & (fz <= hi_))[0]
            if len(band) < 3:
                continue
            j = band[np.argmax(S[band])]
            if 0 < j < len(S) - 1:
                a_, b_, c_ = np.log(S[j - 1] + 1e-12), np.log(S[j] + 1e-12), np.log(S[j + 1] + 1e-12)
                d_ = 0.5 * (a_ - c_) / (a_ - 2 * b_ + c_) if (a_ - 2 * b_ + c_) != 0 else 0.0
                ns_.append(n); fs_.append(fz[j] + d_ * dz)
        if len(ns_) >= 6:
            ns_ = np.array(ns_, float); y_ = (np.array(fs_) / ns_) ** 2
            for _ in range(2):
                A_ = np.vstack([np.ones_like(ns_), ns_ ** 2]).T
                c0, c1 = np.linalg.lstsq(A_, y_, rcond=None)[0]
                r_ = y_ - (c0 + c1 * ns_ ** 2)
                keep = np.abs(r_) <= 2.5 * np.median(np.abs(r_)) + 1e-12
                if keep.sum() < 6 or keep.all():
                    break
                ns_, y_ = ns_[keep], y_[keep]
            if c0 > 0:
                B = max(0.0, float(c1 / c0))

    def partial_track(k):
        band = slice(max(0, k - 2), k + 3)
        return P[band, :].max(axis=0)

    tE = ft - nfft / (2 * rate)  # time of each frame's start, from the onset
    out["f1"] = f1
    out["B"] = B
    if found:
        # Harmonic level: the partials summed in narrow bands. Unlike the broadband level
        # this stays well clear of the recordings' hiss, so late decay can be measured too.
        tracks = np.array([partial_track(k) for _, _, k in found[:20]])
        harm = level_db(tracks.sum(axis=0))
        out["env"] = [float(v) for v in harm[::5]]  # every 50 ms, for plots
        sel_e = (tE >= 0.1) & (tE <= 1.0)
        sel_l = (tE >= late_start) & (tE <= late_end)
        out["early"] = fit_slope(tE[sel_e], harm[sel_e], harm.max() - 70) if sel_e.sum() > 5 else None
        out["late"] = fit_slope(tE[sel_l], harm[sel_l], harm.max() - 70) if sel_l.sum() > 5 else None
        tr = level_db(tracks[0])
        out["f1e"] = fit_slope(tE[sel_e], tr[sel_e], tr.max() - 70) if sel_e.sum() > 5 else None
        out["f1l"] = fit_slope(tE[sel_l], tr[sel_l], tr.max() - 70) if sel_l.sum() > 5 else None
        lv = {n: level_db(spec_early[max(0, k - 2):k + 3].max()) for n, _, k in found}
        p1 = lv.get(1)
        upper = [lv[n] for n in range(2, 9) if n in lv]
        if p1 is not None and upper:
            out["tilt"] = float(level_db(np.sum(10 ** (np.array(upper) / 10))) - p1)
        # decay of each partial over the first second (for the loss curve)
        out["partials"] = []
        for n, fk, k in found[:24]:
            tr = level_db(partial_track(k))
            s = fit_slope(tE[sel_e], tr[sel_e], tr.max() - 60) if sel_e.sum() > 5 else None
            out["partials"].append((n, float(fk), None if s is None else float(s), float(lv[n])))

    # Brightness: ignore the knock and rumble below the fundamental
    spec_early = np.where(fr >= 0.7 * f1, spec_early, 0)
    tot = spec_early.sum()
    out["centroid"] = float((fr * spec_early).sum() / tot / f1)
    out["hi"] = float(level_db(spec_early[fr >= 2000].sum() / tot)) if f1 < 1000 else None
    out["level"] = float(level_db(np.mean(x[:int(0.2 * rate)] ** 2)))
    return out


def collect(salamander, model, notes, velocities):
    rows = []
    for midi in notes:
        for vel in velocities:
            layer = VELOCITY_LAYER[vel]
            sp = Path(salamander) / "samples" / f"{name_of(midi)}v{layer}.flac"
            mp = Path(model) / f"n{midi}_v{vel}.wav"
            if not sp.exists() or not mp.exists():
                continue
            xs, rs = load(sp)
            xm, rm = load(mp)
            rows.append({"midi": midi, "vel": vel, "real": analyze(xs, rs, midi), "model": analyze(xm, rm, midi)})
    return rows


def fmt(v, digits=1):
    return "  -  " if v is None else f"{v:.{digits}f}"


def report(rows, title):
    ref = {}
    for r in rows:
        if r["midi"] == 60 and r["vel"] == 92:
            ref = {"real": r["real"]["level"], "model": r["model"]["level"]}
    lines = [f"# {title}", "",
             "Each cell is real / model. Decay in dB/s (higher = dies faster). Centroid in multiples of f1.",
             "Level is relative to C4 at velocity 92 on the same instrument.", "",
             "| note | vel | early | late | f1 early | f1 late | centroid | >2k dB | tilt dB | level dB |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        a, b = r["real"], r["model"]
        la = a["level"] - ref.get("real", 0)
        lb = b["level"] - ref.get("model", 0)
        lines.append(f"| {name_of(r['midi'])} | {r['vel']} | {fmt(a['early'])} / {fmt(b['early'])} | {fmt(a['late'])} / {fmt(b['late'])} "
                     f"| {fmt(a.get('f1e'))} / {fmt(b.get('f1e'))} | {fmt(a.get('f1l'))} / {fmt(b.get('f1l'))} "
                     f"| {fmt(a['centroid'])} / {fmt(b['centroid'])} | {fmt(a.get('hi'))} / {fmt(b.get('hi'))} "
                     f"| {fmt(a.get('tilt'))} / {fmt(b.get('tilt'))} | {fmt(la)} / {fmt(lb)} |")
    return "\n".join(lines)


def summary(rows):
    """Mean absolute and signed (model - real) differences over all rows."""
    keys = ["early", "late", "f1e", "f1l", "centroid", "hi", "tilt"]
    out = {}
    for k in keys:
        d = [r["model"].get(k) - r["real"].get(k) for r in rows
             if r["model"].get(k) is not None and r["real"].get(k) is not None]
        if k == "centroid":  # in percent of an octave... as 100 * log2(model / real)
            d = [100 * math.log2(r["model"][k] / r["real"][k]) for r in rows]
        out[k] = (float(np.mean(np.abs(d))), float(np.mean(d))) if d else (None, None)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--salamander", required=True)
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", default="report/calibrate")
    ap.add_argument("--notes", default="21,24,27,30,33,36,39,42,45,48,51,54,57,60,63,66,69,72,75,78,81,84,87,90,93,96,99,102,105,108")
    ap.add_argument("--vel", default="40,60,92,124")
    ap.add_argument("--title", default="Salamander vs model")
    args = ap.parse_args()
    notes = [int(v) for v in args.notes.split(",")]
    vels = [int(v) for v in args.vel.split(",")]
    rows = collect(args.salamander, args.model, notes, vels)
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    (out / "metrics.json").write_text(json.dumps(rows, indent=1))
    (out / "report.md").write_text(report(rows, args.title))
    s = summary(rows)
    print(f"{len(rows)} note/velocity pairs")
    for k, (mae, bias) in s.items():
        if mae is not None:
            print(f"  {k:9s} mean |diff| {mae:6.2f}   mean (model - real) {bias:+6.2f}")


if __name__ == "__main__":
    main()
