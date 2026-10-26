#!/usr/bin/env python3
"""Soundboard modes for the piano model, from a finite-difference plate model.

    python3 tools/soundboard/build_board.py [--out engine/src/piano/BoardModes.inc] [--plot board.png]

The board is a spruce plate with the outline of a grand's soundboard (straight spine on the
bass side, curved bentside, narrow tail), clamped all round where it is glued into the rim.
It is thickest in the middle (about 9 mm) and thinner towards the edges (about 6 mm).
The grain runs across the board (x); ribs glued across the grain run along the strings (y),
spaced about 13 cm apart, so above about 1 kHz the vibration is trapped between them.
The treble and bass bridges are hardwood beams glued on top: they stiffen and load the board
along their curves.

It solves K v = w^2 M v for the modes below fTop and writes, for each mode, its frequency,
its shape at the bridge zones (where each group of keys presses), its far-field radiation
towards two listening points and the net volume it pushes. Shapes are scaled so that the
mode's kinetic energy is 1 (unit modal mass), which ModalPiano then rescales to the measured
bridge admittance. The crown and string downbearing (which raise the low modes) and the
damping are applied in ModalPiano, as they change with the piano's age.

Needs numpy and scipy.
"""
import argparse, math
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy import ndimage

DX = 0.02                      # grid step [m]
A, B = 1.42, 1.62              # bounding box: across (treble 0 .. bass A), along the strings (keyboard 0 .. tail B)
E_L, E_T, G_LT, NU = 11e9, 0.75e9, 0.7e9, 0.4   # spruce: along / across the grain, shear, Poisson (LT)
RHO = 420.0
H_MID, H_EDGE = 0.009, 0.006
RIB_SPACING, RIB_W, RIB_T = 0.13, 0.022, 0.025
RIB_EI = 2.0 * E_L * RIB_W * RIB_T ** 3 / 12.0     # x2: glued on the plate, the neutral axis moves out
RIB_M = RHO * RIB_W * RIB_T
TREBLE_EI, TREBLE_M = 14e9 * 0.030 * 0.030 ** 3 / 12.0, 750.0 * 0.030 * 0.030   # beech bridge
BASS_EI, BASS_M = 14e9 * 0.035 * 0.035 ** 3 / 12.0, 750.0 * 0.035 * 0.035
F_TOP = 2000.0
N_BASS, N_TREBLE = 3, 13       # bridge zones (bass bridge keys A0..E2, treble bridge F2..C8)


def bentside(y):
    """Treble-side edge of the board: straight near the keyboard, then curving in to the tail."""
    t = np.clip((y - 0.25) / (B - 0.25), 0.0, 1.0)
    return 1.10 * t ** 1.8


def inside(x, y):
    return (y >= 0) & (y <= B) & (x <= A) & (x >= bentside(y))


def treble_bridge(t):
    """t = 0 at the top note (C8, near the keyboard on the treble side) .. 1 at F2."""
    x = 0.15 + 0.70 * t
    y = 0.10 + 0.78 * t ** 1.15 + 0.04 * np.sin(np.pi * t)
    return x, y


def bass_bridge(t):
    """t = 0 at E2 .. 1 at A0, near the tail on the bass side."""
    return 1.02 + 0.23 * t, 1.10 + 0.38 * t


def curve_points(fn, n):
    t = np.linspace(0, 1, 400)
    x, y = fn(t)
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    si = np.linspace(0, s[-1], n)
    return np.interp(si, s, x), np.interp(si, s, y), s[-1]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="engine/src/piano/BoardModes.inc")
    ap.add_argument("--plot")
    args = ap.parse_args()

    nx, ny = int(round(A / DX)) + 1, int(round(B / DX)) + 1
    X, Y = np.meshgrid(np.arange(nx) * DX, np.arange(ny) * DX, indexing="ij")
    board = inside(X, Y)
    # Clamped edge: the outermost ring of board nodes is held (w = 0), and so is everything outside
    free = ndimage.binary_erosion(board, iterations=1)
    idx = -np.ones((nx, ny), int)
    idx[free] = np.arange(free.sum())
    nf = int(free.sum())
    P = sp.csr_matrix((np.ones(nf), (np.flatnonzero(free.ravel()), np.arange(nf))), shape=(nx * ny, nf))

    # Thickness: 9 mm in the middle, 6 mm at the edge
    dist = ndimage.distance_transform_edt(board) * DX
    h = H_EDGE + (H_MID - H_EDGE) * np.clip(dist / dist.max(), 0, 1) ** 0.7
    Dx = E_L * h ** 3 / 12.0
    Dy = E_T * h ** 3 / 12.0
    D1 = NU * E_T * h ** 3 / 12.0
    Dk = G_LT * h ** 3 / 12.0
    mass = RHO * h * DX * DX
    # Ribs along y, every RIB_SPACING across the board
    for xr in np.arange(0.10, A - 0.05, RIB_SPACING):
        i = int(round(xr / DX))
        col = board[i, :]
        Dy[i, col] += RIB_EI / DX
        mass[i, col] += RIB_M * DX

    def lin(i, j):
        return i * ny + j

    # Curvature operators on the full grid, evaluated at every board node
    rows, cols, vals = [], [], []
    def op(stencil):
        r, c, v = [], [], []
        k = 0
        pts = np.argwhere(board)
        for (i, j) in pts:
            for (di, dj, w) in stencil:
                ii, jj = i + di, j + dj
                if 0 <= ii < nx and 0 <= jj < ny:
                    r.append(k); c.append(lin(ii, jj)); v.append(w / DX ** 2)
            k += 1
        return sp.csr_matrix((v, (r, c)), shape=(len(pts), nx * ny)), pts
    Cxx, pts = op([(-1, 0, 1), (0, 0, -2), (1, 0, 1)])
    Cyy, _ = op([(0, -1, 1), (0, 0, -2), (0, 1, 1)])
    Cxy, _ = op([(0, 0, 0.25 * 4), (1, 1, 0.25 * 4), (1, 0, -0.25 * 4), (0, 1, -0.25 * 4)])  # w_xy at the cell corner, x4/4
    Cxy = Cxy * 0.25
    area = DX * DX
    g = lambda F: sp.diags(F[pts[:, 0], pts[:, 1]] * area)
    K = (Cxx.T @ g(Dx) @ Cxx + Cyy.T @ g(Dy) @ Cyy + Cxx.T @ g(D1) @ Cyy + Cyy.T @ g(D1) @ Cxx
         + Cxy.T @ g(4 * Dk) @ Cxy)
    M = sp.diags(mass.ravel())

    # Bridges: beams along their curves, sampled every DX, bilinear onto the grid
    def interp_rows(xs, ys):
        r, c, v = [], [], []
        for k, (x, y) in enumerate(zip(xs, ys)):
            fx, fy = x / DX, y / DX
            i0, j0 = int(fx), int(fy)
            ax, ay = fx - i0, fy - j0
            for di, dj, w in ((0, 0, (1 - ax) * (1 - ay)), (1, 0, ax * (1 - ay)), (0, 1, (1 - ax) * ay), (1, 1, ax * ay)):
                r.append(k); c.append(lin(i0 + di, j0 + dj)); v.append(w)
        return sp.csr_matrix((v, (r, c)), shape=(len(xs), nx * ny))
    for fn, EI, mlin in ((treble_bridge, TREBLE_EI, TREBLE_M), (bass_bridge, BASS_EI, BASS_M)):
        _, _, L = curve_points(fn, 2)
        n = max(3, int(L / DX) + 1)
        xs, ys, L = curve_points(fn, n)
        ds = L / (n - 1)
        Bm = interp_rows(xs, ys)
        D2 = sp.diags([np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)], [0, 1, 2], shape=(n - 2, n)) / ds ** 2
        C = D2 @ Bm
        K = K + C.T @ C * (EI * ds)
        M = M + Bm.T @ sp.diags(np.full(n, mlin * ds)) @ Bm

    Kf = (P.T @ K @ P).tocsc()
    Mf = (P.T @ M @ P).tocsc()
    lam_top = (2 * math.pi * F_TOP) ** 2
    nev = 400
    vals, vecs = sla.eigsh(Kf, k=nev, M=Mf, sigma=0, which="LM")
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    keep = vals < lam_top
    vals, vecs = vals[keep], vecs[:, keep]
    f = np.sqrt(np.maximum(vals, 0)) / (2 * math.pi)
    print(f"{len(f)} modes below {F_TOP:.0f} Hz (asked {nev}); lowest {f[:6].round(1)}; "
          f"mean spacing {np.mean(np.diff(f)):.1f} Hz, below 1 kHz {np.mean(np.diff(f[f < 1000])):.1f} Hz")
    if keep.all():
        print("warning: raise nev, the top of the range is missing")
    # unit modal mass (eigsh returns M-orthonormal vectors already)
    full = P @ vecs                                         # (nx*ny, n)

    # Bridge zone points: the middle of each group of keys along each bridge
    zx, zy = [], []
    bx, by, _ = curve_points(bass_bridge, 2 * N_BASS + 1)
    for z in range(N_BASS):          # zone 0 = lowest keys (A0 end, t = 1)
        k = 2 * (N_BASS - 1 - z) + 1
        zx.append(bx[k]); zy.append(by[k])
    tx, ty, _ = curve_points(treble_bridge, 2 * N_TREBLE + 1)
    for z in range(N_TREBLE):        # from F2 (t = 1) up to C8 (t = 0)
        k = 2 * (N_TREBLE - 1 - z) + 1
        zx.append(tx[k]); zy.append(ty[k])
    Z = interp_rows(np.array(zx), np.array(zy))
    phi = (Z @ full)                                        # (zones, n)

    # Radiation: far-field (Rayleigh) integral towards two points above the keyboard
    dirL, dirR = (0.35, -0.5), (-0.35, -0.5)
    xs, ys = X.ravel(), Y.ravel()
    k = 2 * math.pi * f / 343.0
    radL = np.zeros(len(f)); radR = np.zeros(len(f)); vol = np.zeros(len(f))
    for i in range(len(f)):
        w = full[:, i]
        nz = w != 0
        pl = k[i] * (dirL[0] * xs[nz] + dirL[1] * ys[nz])
        pr = k[i] * (dirR[0] * xs[nz] + dirR[1] * ys[nz])
        lr, li = np.sum(w[nz] * np.cos(pl)), np.sum(w[nz] * np.sin(pl))
        rr, ri = np.sum(w[nz] * np.cos(pr)), np.sum(w[nz] * np.sin(pr))
        radL[i] = math.hypot(lr, li) * (-1 if lr < 0 else 1)
        radR[i] = math.hypot(rr, ri) * (-1 if rr < 0 else 1)
        vol[i] = np.sum(w[nz])

    def lit(v):
        t = "%.5g" % v
        if "." not in t and "e" not in t:
            t += ".0"
        return t + "f"

    with open(args.out, "w") as fo:
        fo.write("// Generated by tools/soundboard/build_board.py: do not edit by hand.\n")
        fo.write(f"// Soundboard modes below {F_TOP:.0f} Hz of a clamped, tapered, ribbed grand soundboard with bridges.\n")
        fo.write("// Per mode: frequency [Hz], shape at the 16 bridge zones (unit modal mass), radiation towards\n")
        fo.write("// the left and right listening points, and net volume (both in arbitrary units).\n")
        fo.write(f"static constexpr int kBoardModeCount = {len(f)};\n")
        fo.write(f"static constexpr int kBoardZones = {N_BASS + N_TREBLE};\n")
        fo.write("struct BoardModeData { float f; float phi[%d]; float radLeft, radRight, volume; };\n" % (N_BASS + N_TREBLE))
        fo.write("static const BoardModeData kBoardModes[kBoardModeCount] = {\n")
        for i in range(len(f)):
            fo.write("    { %.2ff, { %s }, %s, %s, %s },\n" % (
                f[i], ", ".join(lit(v) for v in phi[:, i]), lit(radL[i]), lit(radR[i]), lit(vol[i])))
        fo.write("};\n")
    print("wrote", args.out)

    if args.plot:
        import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
        fig, ax = plt.subplots(3, 4, figsize=(14, 12))
        sel = [0, 1, 2, 3, 5, 8, 12, 20, 40, 80, 150, len(f) - 1]
        for a_, i in zip(ax.ravel(), sel):
            a_.imshow(full[:, i].reshape(nx, ny).T, origin="lower", cmap="RdBu", extent=[0, A, 0, B])
            a_.plot(zx, zy, "k.", ms=3); a_.set_title(f"{f[i]:.0f} Hz")
        fig.savefig(args.plot, dpi=70)


if __name__ == "__main__":
    main()
