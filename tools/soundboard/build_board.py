#!/usr/bin/env python3
"""Soundboard modes for the piano model, from a finite-difference plate model.

    python3 tools/soundboard/build_board.py [--out engine/src/piano/BoardModes.inc] [--plot board.png]

The board is a spruce plate with the outline of a grand's soundboard (straight spine on the
bass side, curved bentside, narrow tail), glued into a rim that is itself modelled: the rim bends, carries the iron frame, is
held by back braces and stands on three legs, and the air in the case above the board is coupled in.
It is thickest in the middle (about 9 mm) and thinner towards the edges (about 6 mm).
The grain runs diagonally (35 degrees from the keyboard, roughly along the treble bridge);
ribs glued across the grain are spaced about 13 cm apart, so above about 1 kHz the vibration is trapped between them.
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

DX = 0.01                      # grid step [m]
A, B = 1.38, 1.47              # Yamaha C5 board's bounding box: across (treble 0 .. bass A), along the strings (belly rail 0 .. tail B)
E_L, E_T, G_LT, NU = 11e9, 0.75e9, 0.7e9, 0.4   # spruce: along / across the grain, shear, Poisson (LT)
GRAIN = math.radians(35.0)    # the grain runs diagonally, from the treble front towards the bass tail
RHO = 420.0
H_MID, H_EDGE = 0.009, 0.006
RIB_SPACING, RIB_W, RIB_T = 0.13, 0.022, 0.025
RIB_EI = 2.0 * E_L * RIB_W * RIB_T ** 3 / 12.0     # x2: glued on the plate, the neutral axis moves out
RIB_M = RHO * RIB_W * RIB_T
TREBLE_EI, TREBLE_M = 14e9 * 0.030 * 0.030 ** 3 / 12.0, 750.0 * 0.030 * 0.030   # beech bridge
BASS_EI, BASS_M = 14e9 * 0.035 * 0.035 ** 3 / 12.0, 750.0 * 0.035 * 0.035
F_TOP = 4000.0
# The body around the board (a Yamaha C5): the rim it is glued into, the braces under it,
# the iron frame bolted on the rim, the legs, and the air in the case above the board.
RIM_EI = 12e9 * 0.035 * 0.26 ** 3 / 12.0 * 3.0   # laminated maple, 26 cm deep; x3: the bolted iron frame ties it
RIM_M = 650.0 * 0.035 * 0.26                      # [kg/m]
FRAME_MASS, FRAME_BOLTS = 125.0, 10               # iron frame [kg], spread over its bolts on the rim
EDGE_K = 2.0e6                                    # glue joint holding the board's edge level with the rim [N/m per node]
BRACE_EI, BRACE_M = 12e9 * 0.09 * 0.11 ** 3 / 12.0, 650.0 * 0.09 * 0.11   # back braces (beech/spruce, 9x11 cm)
LEG_K = 3.0e6                                     # each leg on its castor [N/m]
CASE_H = 0.24 + 0.08                              # rim height above the board, plus the open top's end correction [m]
RHO_AIR, C_AIR = 1.2, 343.0
FLOOR_GAP = 0.62           # board underside to the floor [m] (the rim skirt hangs ~10 cm below the board, then it is open)
# Crown and downbearing. The ribs are planed to an arch (radius about 20 m), so the board is a
# shallow shell: bending it also stretches it, which stiffens the low modes. The strings press
# the bridges down (downbearing, about 1.7% of the total string tension, ~3 kN on a C5); that
# load flattens the arch where it pushes (less crown, less stiffening there) and squeezes the
# board along the ribs (compression, which softens bending). Both are solved here, not guessed.
CROWN_R = 20.0
DOWNBEARING = 3000.0       # [N], shared along both bridges (a quarter on the bass bridge)
N_BASS, N_TREBLE = 3, 13       # bridge zones (bass bridge keys A0..E2, treble bridge F2..C8)


def bentside(y):
    """Treble-side edge of the board: straight near the keyboard, then curving in to the tail."""
    t = np.clip((y - 0.25) / (B - 0.25), 0.0, 1.0)
    return 1.10 * t ** 1.8


def inside(x, y):
    return (y >= 0) & (y <= B) & (x <= A) & (x >= bentside(y))


def c5_lengths():
    """Speaking lengths per key (1..88), the same table the engine uses (ScaleDesign.cpp)."""
    import re, pathlib
    src = (pathlib.Path(__file__).resolve().parents[2] / "engine/src/piano/ScaleDesign.cpp").read_text()
    d = [float(v) for v in re.search(r"kSteinwayD\[88\] = \{(.*?)\};", src, re.S).group(1).replace("\n", " ").split(",") if v.strip()]
    def scale(k):
        if k >= 40: return 1.0
        if k > 20: return 0.80 + 0.20 * (k - 21) / 19.0
        return 0.66 + 0.06 * (k - 1) / 19.0
    return np.array([d[k - 1] * scale(k) for k in range(1, 89)])


LEN = c5_lengths()
CAPO_Y = 0.02     # the treble strings' speaking length starts just in front of the belly rail
FAN = 0.86        # the strings fan out a little towards the tail: along-the-board share of each length


def treble_bridge(t):
    """t = 0 at the top note (C8, near the belly rail on the treble side) .. 1 at F2.
    Each key's bridge pin sits one speaking length behind its capo/agraffe."""
    k = 88 - t * (88 - 21)
    L = np.interp(k, np.arange(1, 89), LEN)
    x = 0.12 + 0.80 * t
    return x, CAPO_Y + FAN * L


def bass_bridge(t):
    """t = 0 at E2 .. 1 at A0: the cross-strung bass runs at about 25 degrees over the tenor,
    towards the spine and the tail."""
    k = 20 - t * 19
    L = np.interp(k, np.arange(1, 89), LEN)
    ang = math.radians(25.0)
    return 1.00 + 0.20 * t, 0.10 + L * math.cos(ang)


def curve_points(fn, n):
    t = np.linspace(0, 1, 400)
    x, y = fn(t)
    s = np.concatenate([[0], np.cumsum(np.hypot(np.diff(x), np.diff(y)))])
    si = np.linspace(0, s[-1], n)
    return np.interp(si, s, x), np.interp(si, s, y), s[-1]


# --- The cast-iron frame (plate) ----------------------------------------------------------
# The frame holds the strings' far ends: the agraffes and capo bar at the front. It is a
# network of cast-iron bars (front bar over the pin block, bentside and spine flanges, the
# tail rail with the hitch pins, and three struts), bolted onto the rim. Its bending modes
# move the strings' far ends; the shapes at the 16 key zones along the front bar are what
# the engine needs.
IRON_E, IRON_RHO = 1.1e11, 7200.0
FRAME_A, FRAME_I = 0.0022, 2.6e-6          # bar section [m^2], bending [m^4] (8 cm deep ribbed bars)
FRONT_Y = -0.05                             # the front bar, just ahead of the belly rail
BOLT_K = 2.0e7                              # each bolt into the rim [N/m]


def frame_modes(zone_x, zone_is_bass, f_top):
    bars = [
        [(0.05, FRONT_Y), (1.35, FRONT_Y)],                                   # front bar (capo, agraffes)
        [(1.35, FRONT_Y), (1.35, 1.40)],                                      # spine flange
        [(0.98, 1.40), (1.35, 1.40)],                                         # tail rail (hitch pins)
        [(0.05, FRONT_Y)] + [(float(bentside(y)) + 0.04, float(y)) for y in np.linspace(0.05, 1.40, 12)],  # bentside flange
        [(0.35, FRONT_Y), (0.98, 1.38)], [(0.75, FRONT_Y), (1.15, 1.40)], [(1.05, FRONT_Y), (1.32, 1.00)],  # struts
    ]
    h = 0.04
    nodes, coords = [], []
    def node_at(x, y):
        for k, (a, b) in enumerate(coords):
            if abs(a - x) < 0.03 and abs(b - y) < 0.03:
                return k
        coords.append((x, y)); return len(coords) - 1
    chains = []
    for pts_ in bars:
        P_ = np.array(pts_)
        seg = np.hypot(*np.diff(P_, axis=0).T); s_ = np.concatenate([[0], np.cumsum(seg)])
        n_ = max(3, int(s_[-1] / h) + 1)
        si = np.linspace(0, s_[-1], n_)
        xs = np.interp(si, s_, P_[:, 0]); ys = np.interp(si, s_, P_[:, 1])
        chains.append(([node_at(x, y) for x, y in zip(xs, ys)], s_[-1] / (n_ - 1)))
    n = len(coords)
    K = np.zeros((n, n)); Mv = np.zeros(n)
    for ch, ds in chains:
        for a, b, c in zip(ch[:-2], ch[1:-1], ch[2:]):
            v = np.zeros(n); v[a] += 1; v[b] -= 2; v[c] += 1
            K += np.outer(v, v) * IRON_E * FRAME_I / ds ** 3
        for k in ch:
            Mv[k] += IRON_RHO * FRAME_A * ds
    C = np.array(coords)
    # Bolts: every ~30 cm along the bars that lie on the rim (front ends, spine, bentside, tail)
    on_rim = (np.abs(C[:, 0] - 1.35) < 0.02) | (C[:, 1] > 1.38) | (np.abs(C[:, 0] - (bentside(C[:, 1]) + 0.04)) < 0.03) | \
             ((np.abs(C[:, 1] - FRONT_Y) < 0.01) & ((C[:, 0] < 0.1) | (C[:, 0] > 1.3)))
    bolts = np.flatnonzero(on_rim)[::6]
    for b in bolts:
        K[b, b] += BOLT_K
    vals, vecs = np.linalg.eigh(K / np.sqrt(Mv[:, None] * Mv[None, :]))
    vecs = vecs / np.sqrt(Mv)[:, None]                       # unit modal mass
    keep = (vals > (2 * math.pi * 20) ** 2) & (vals < (2 * math.pi * f_top) ** 2)
    f = np.sqrt(vals[keep]) / (2 * math.pi); V = vecs[:, keep]
    front = [k for k in range(n) if abs(C[k, 1] - FRONT_Y) < 0.01]
    front = sorted(front, key=lambda k: C[k, 0])
    fx = C[front, 0]
    # Each zone's strings end on the front bar: treble zones straight ahead of their bridge
    # section, the cross-strung bass towards the spine end
    ax = np.where(zone_is_bass, 1.25 - 0.08 * np.arange(len(zone_x)), zone_x)
    phi = np.array([[np.interp(x, fx, V[front, i]) for x in ax] for i in range(V.shape[1])])
    rad = np.array([np.sum(V[:, i] * Mv) for i in range(V.shape[1])])
    rad = rad / (np.abs(rad).max() + 1e-12)
    print(f"frame: {len(coords)} nodes, {Mv.sum():.0f} kg, {len(bolts)} bolts; {len(f)} modes below {f_top:.0f} Hz, lowest {f[:6].round(0)}")
    return f, phi, rad


def write_frame(path, f, phi, rad):
    def lit(v):
        t = "%.5g" % v
        if "." not in t and "e" not in t:
            t += ".0"
        return t + "f"
    with open(path, "w") as fo:
        fo.write("// Generated by tools/soundboard/build_board.py: do not edit by hand.\n")
        fo.write("// Bending modes of the cast-iron frame (bars bolted on the rim). Per mode: frequency [Hz], shape at the\n")
        fo.write("// strings' far ends for the 16 key zones (unit modal mass, 1/sqrt(kg)), and a relative radiation weight.\n")
        fo.write(f"static constexpr int kFrameModeCount = {len(f)};\n")
        fo.write("struct FrameModeData { float f; float phi[%d]; float rad; };\n" % phi.shape[1])
        fo.write("static const FrameModeData kFrameModes[kFrameModeCount] = {\n")
        for i in range(len(f)):
            fo.write("    { %.2ff, { %s }, %s },\n" % (f[i], ", ".join(lit(v) for v in phi[i]), lit(rad[i])))
        fo.write("};\n")
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="engine/src/piano/BoardModes.inc")
    ap.add_argument("--frame-out", default="engine/src/piano/FrameModes.inc")
    ap.add_argument("--plot")
    ap.add_argument("--json", help="also write the board geometry and some mode shapes for the visualiser")
    args = ap.parse_args()

    nx, ny = int(round(A / DX)) + 1, int(round(B / DX)) + 1
    X, Y = np.meshgrid(np.arange(nx) * DX, np.arange(ny) * DX, indexing="ij")
    board = inside(X, Y)
    # The outermost ring of board nodes sits on the rim. The rim is not rigid: it is a deep
    # curved beam that bends, carries the iron frame and stands on three legs, so the ring
    # nodes move (with the rim's mass and stiffness added to them) instead of being held.
    interior = ndimage.binary_erosion(board, structure=np.ones((3, 3)), iterations=1)
    ring = board & ~interior
    free = board
    nf = int(free.sum())
    P = sp.csr_matrix((np.ones(nf), (np.flatnonzero(free.ravel()), np.arange(nf))), shape=(nx * ny, nf))

    # Thickness: 9 mm in the middle, 6 mm at the edge
    dist = ndimage.distance_transform_edt(board) * DX
    h = H_EDGE + (H_MID - H_EDGE) * np.clip(dist / dist.max(), 0, 1) ** 0.7
    D11 = E_L * h ** 3 / 12.0          # along the grain
    D22 = E_T * h ** 3 / 12.0          # across it
    D12 = NU * E_T * h ** 3 / 12.0
    D66 = G_LT * h ** 3 / 12.0
    mass = RHO * h * DX * DX

    def lin(i, j):
        return i * ny + j

    # Curvature operators on the full grid, evaluated at every board node
    rows, cols, vals = [], [], []
    def op(stencil):
        r, c, v = [], [], []
        k = 0
        pts = np.argwhere(interior)
        for (i, j) in pts:
            for (di, dj, w) in stencil:
                ii, jj = i + di, j + dj
                if 0 <= ii < nx and 0 <= jj < ny:
                    r.append(k); c.append(lin(ii, jj)); v.append(w / DX ** 2)
            k += 1
        return sp.csr_matrix((v, (r, c)), shape=(len(pts), nx * ny)), pts
    Cxx, pts = op([(-1, 0, 1), (0, 0, -2), (1, 0, 1)])
    Cyy, _ = op([(0, -1, 1), (0, 0, -2), (0, 1, 1)])
    Cxy, _ = op([(0, 0, 1), (1, 1, 1), (1, 0, -1), (0, 1, -1)])   # w_xy at the cell corner
    area = DX * DX
    g = lambda F: sp.diags(F[pts[:, 0], pts[:, 1]] * area)
    # Curvatures in the grain's axes: w_11 (along), w_22 (across) and 2 w_12
    c, s_ = math.cos(GRAIN), math.sin(GRAIN)
    C1 = c * c * Cxx + s_ * s_ * Cyy + 2 * c * s_ * Cxy
    C2 = s_ * s_ * Cxx + c * c * Cyy - 2 * c * s_ * Cxy
    C3 = -2 * c * s_ * Cxx + 2 * c * s_ * Cyy + 2 * (c * c - s_ * s_) * Cxy
    K = (C1.T @ g(D11) @ C1 + C2.T @ g(D22) @ C2 + C1.T @ g(D12) @ C2 + C2.T @ g(D12) @ C1
         + C3.T @ g(D66) @ C3)
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
    beams = [(treble_bridge, TREBLE_EI, TREBLE_M), (bass_bridge, BASS_EI, BASS_M)]
    # Ribs: straight beams across the grain, every RIB_SPACING along it, cut to the board
    gdir = np.array([c, s_]); rdir = np.array([-s_, c])
    centre = np.array([0.85, 0.75])
    for off in np.arange(-1.2, 1.2, RIB_SPACING):
        t = np.linspace(-1.5, 1.5, 1500)
        px = centre[0] + off * gdir[0] + t * rdir[0]
        py = centre[1] + off * gdir[1] + t * rdir[1]
        ok = inside(px, py) & (px > DX) & (py > DX) & (px < A - DX) & (py < B - DX)
        if ok.sum() < 10:
            continue
        ii = np.flatnonzero(ok)
        t0, t1 = t[ii[0]] + 0.03, t[ii[-1]] - 0.03   # ribs stop just short of the rim
        if t1 - t0 < 0.15:
            continue
        beams.append(((lambda u, o=off, a0=t0, a1=t1: (centre[0] + o * gdir[0] + (a0 + (a1 - a0) * u) * rdir[0],
                                                     centre[1] + o * gdir[1] + (a0 + (a1 - a0) * u) * rdir[1])),
                      RIB_EI, RIB_M))
    print(f"{len(beams) - 2} ribs")
    beam_lines = []
    for fn, EI, mlin in beams:
        u = np.linspace(0, 1, 40)
        bxs, bys = fn(u)
        beam_lines.append([[round(float(a_), 3), round(float(b_), 3)] for a_, b_ in zip(bxs, bys)])
    for fn, EI, mlin in beams:
        _, _, L = curve_points(fn, 2)
        n = max(3, int(L / DX) + 1)
        xs, ys, L = curve_points(fn, n)
        ds = L / (n - 1)
        Bm = interp_rows(xs, ys)
        D2 = sp.diags([np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)], [0, 1, 2], shape=(n - 2, n)) / ds ** 2
        C = D2 @ Bm
        K = K + C.T @ C * (EI * ds)
        M = M + Bm.T @ sp.diags(np.full(n, mlin * ds)) @ Bm

    # --- Crown and downbearing (see the constants above)
    cR, sR = -s_, c                                         # rib direction (across the grain)
    def grad_op(dx_, dy_):
        r, cc, v = [], [], []
        for k_, (i, j) in enumerate(pts):
            r += [k_, k_]; cc += [lin(i, j), lin(i + dx_, j + dy_)]; v += [-1.0 / DX, 1.0 / DX]
        return sp.csr_matrix((v, (r, cc)), shape=(len(pts), nx * ny))
    Gr = cR * grad_op(1, 0) + sR * grad_op(0, 1)            # slope along the ribs
    Crr = cR * cR * Cxx + sR * sR * Cyy + 2 * cR * sR * Cxy # curvature along the ribs
    Eh_r = E_T * h + E_L * RIB_W * RIB_T / RIB_SPACING      # membrane stiffness along the ribs [N/m]
    ipts = (pts[:, 0], pts[:, 1])
    # Static sag under the downbearing, board edge held: K_board w = F
    Fd = np.zeros(nx * ny)
    for fn, share in ((treble_bridge, 0.75), (bass_bridge, 0.25)):
        n_ = 60
        xs_, ys_, _ = curve_points(fn, n_)
        Fd -= interp_rows(xs_, ys_).T @ np.full(n_, DOWNBEARING * share / n_)
    found0 = sp.diags((Eh_r[ipts] / CROWN_R ** 2) * area)
    Ks = (C1.T @ g(D11) @ C1 + C2.T @ g(D22) @ C2 + C1.T @ g(D12) @ C2 + C2.T @ g(D12) @ C1 + C3.T @ g(D66) @ C3)
    Ks = Ks + sp.csr_matrix((found0.diagonal(), (np.array([lin(i, j) for i, j in pts]),) * 2), shape=(nx * ny, nx * ny))
    ii_ = np.flatnonzero(interior.ravel())
    Kss = (K - (C1.T @ g(D11) @ C1 + C2.T @ g(D22) @ C2 + C1.T @ g(D12) @ C2 + C2.T @ g(D12) @ C1 + C3.T @ g(D66) @ C3) + Ks).tocsr()[ii_][:, ii_]
    wsag = np.zeros(nx * ny); wsag[ii_] = sla.spsolve(Kss.tocsc(), Fd[ii_])
    kappa = 1.0 / CROWN_R - (Crr @ wsag)                    # the arch curves down (w_rr = -1/R); the sag bends the other way: flatter where pressed
    print(f"downbearing sag: max {(-wsag).max() * 1000:.2f} mm; crown radius left under the bridge: "
          f"{1 / np.clip(np.median(kappa[np.abs(Gr @ wsag) >= 0]), 1e-4, None):.0f} m (median), "
          f"{1 / np.clip(kappa.min(), 1e-4, None):.0f} m (flattest)")
    kappa = np.clip(kappa, 0.0, None)
    Nc = DOWNBEARING / (A * B * 0.75) * CROWN_R              # membrane compression along the ribs, p R [N/m]
    K = K + sp.csr_matrix(((Eh_r[ipts] * kappa ** 2) * area, (np.array([lin(i, j) for i, j in pts]),) * 2), shape=(nx * ny, nx * ny))
    K = K - Gr.T @ sp.diags(np.full(len(pts), Nc * area)) @ Gr
    print(f"crown stiffening ~{math.sqrt(np.median(Eh_r[ipts] * kappa ** 2) / (RHO * H_MID * 1.4)) / (2 * math.pi):.0f} Hz equivalent; compression {Nc / 1000:.0f} kN/m")

    # --- The rim: ring nodes ordered round the outline, a beam bending out of the board's plane
    rp = np.argwhere(ring)
    cx, cy = rp[:, 0].mean(), rp[:, 1].mean()
    rp = rp[np.argsort(np.arctan2(rp[:, 1] - cy, rp[:, 0] - cx))]
    nr = len(rp)
    rlin = np.array([lin(i, j) for i, j in rp])
    seg = np.hypot(*(np.roll(rp, -1, 0) - rp).T) * DX
    ds = float(seg.mean())
    R = sp.csr_matrix((np.ones(nr), (np.arange(nr), rlin)), shape=(nr, nx * ny))
    D2r = (sp.eye(nr, k=-1) - 2 * sp.eye(nr) + sp.eye(nr, k=1)).tolil()
    D2r[0, nr - 1] = 1; D2r[nr - 1, 0] = 1                 # closed loop
    Cr = D2r.tocsr() / ds ** 2 @ R
    K = K + Cr.T @ Cr * (RIM_EI * ds)
    rimMass = np.full(nr, RIM_M * ds)
    rimMass[np.linspace(0, nr, FRAME_BOLTS, endpoint=False).astype(int)] += FRAME_MASS / FRAME_BOLTS
    M = M + R.T @ sp.diags(rimMass) @ R
    # The glue joint: the board's edge leaves the rim level (a stiff hinge spring to its inner neighbour)
    er, ec, ev = [], [], []
    for (i, j) in np.argwhere(ring):
        nb = [(i + a, j + b) for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)) if 0 <= i + a < nx and 0 <= j + b < ny and interior[i + a, j + b]]
        for (ii, jj) in nb[:1]:
            r_ = len(er) // 2
            er += [r_, r_]; ec += [lin(i, j), lin(ii, jj)]; ev += [1.0, -1.0]
    Ej = sp.csr_matrix((ev, (er, ec)), shape=(len(er) // 2, nx * ny))
    K = K + Ej.T @ Ej * EDGE_K
    # Legs: two at the keyboard corners, one under the tail
    def nearest_ring(x, y):
        d = np.hypot(rp[:, 0] * DX - x, rp[:, 1] * DX - y)
        return int(np.argmin(d))
    legs = [nearest_ring(A, 0.0), nearest_ring(0.0, 0.0), nearest_ring(0.75, B)]
    for r_ in legs:
        K = K + sp.csr_matrix(([LEG_K], ([rlin[r_]], [rlin[r_]])), shape=(nx * ny, nx * ny))

    # --- Back braces: beams under the board (not touching it) from a hub towards the tail out
    # to the rim, each with its own nodes, the ends pinned onto the rim and the hub shared
    hub = np.array([0.70, 1.05])
    ends = [(A, 0.35), (A, 0.95), (0.95, B), (0.35, 1.05), (0.12, 0.45)]
    nb_nodes = 12
    nbr = len(ends) * nb_nodes + 1                          # + the shared hub node
    Kb = sp.lil_matrix((nbr, nbr)); Mb = np.zeros(nbr)
    tie = []                                                # (brace node, ring index)
    for e, (ex, ey) in enumerate(ends):
        L = float(np.hypot(ex - hub[0], ey - hub[1]))
        h_ = L / nb_nodes
        nodes = [nbr - 1] + [e * nb_nodes + k for k in range(nb_nodes)]   # hub .. rim end
        for k in range(1, nb_nodes):
            a, b, c_ = nodes[k - 1], nodes[k], nodes[k + 1]
            v = np.zeros(nbr); v[a] += 1; v[b] -= 2; v[c_] += 1
            Kb += sp.lil_matrix(np.outer(v, v) * BRACE_EI / h_ ** 3)
        for k in nodes:
            Mb[k] += BRACE_M * h_
        tie.append((nodes[-1], nearest_ring(ex, ey)))
    Kb = Kb.tocsr()
    nb_ = nbr
    # Assemble board+rim (free nodes) with the braces
    Kp = (P.T @ K @ P)
    Mp = (P.T @ M @ P)
    ring_free = np.flatnonzero(ring.ravel()[np.flatnonzero(free.ravel())])   # unused, kept for clarity
    Kf = sp.bmat([[Kp, None], [None, Kb]]).tolil()
    Mf = sp.bmat([[Mp, None], [None, sp.diags(Mb)]]).tolil()
    pos_of = -np.ones(nx * ny, int); pos_of[np.flatnonzero(free.ravel())] = np.arange(nf)
    for bn_, r_ in tie:
        a, b = pos_of[rlin[r_]], nf + bn_
        k_ = 50.0 * EDGE_K
        Kf[a, a] += k_; Kf[b, b] += k_; Kf[a, b] -= k_; Kf[b, a] -= k_
    Kf = Kf.tocsc(); Mf = Mf.tocsc()
    rimNodes = pos_of[rlin]
    lam_top = (2 * math.pi * F_TOP) ** 2
    nev = 1800
    vals, vecs = sla.eigsh(Kf, k=nev, M=Mf, sigma=0, which="LM")
    order = np.argsort(vals)
    vals, vecs = vals[order], vecs[:, order]
    keep = (vals < lam_top) & (vals > (2 * math.pi * 5.0) ** 2)
    vals, vecs = vals[keep], vecs[:, keep]
    f = np.sqrt(np.maximum(vals, 0)) / (2 * math.pi)
    print(f"{len(f)} modes below {F_TOP:.0f} Hz (asked {nev}); lowest {f[:6].round(1)}; "
          f"mean spacing {np.mean(np.diff(f)):.1f} Hz, below 1 kHz {np.mean(np.diff(f[f < 1000])):.1f} Hz")
    if keep.all():
        print("warning: raise nev, the top of the range is missing")
    # Energy shares of the structural modes: board, rim (with frame), braces
    rimE = np.array([np.sum(rimMass * vecs[rimNodes, i] ** 2) for i in range(vecs.shape[1])])
    brE = (Mb[:, None] * vecs[nf:, :] ** 2).sum(0)
    vecs = vecs[:nf]

    # --- The air in the case above the board, coupled to the structure. Pressure modes of the
    # outline (rigid rim walls) with a quarter-wave standing up to the open top; the board
    # pumps them and they push back on it (vibro-acoustic coupling, solved together).
    bi = np.flatnonzero(board.ravel())
    nbd = len(bi)
    pos = -np.ones(nx * ny, int); pos[bi] = np.arange(nbd)
    lr_, lc_, lv_ = [], [], []
    for k_, (i, j) in enumerate(np.argwhere(board)):
        for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + a, j + b
            if 0 <= ii < nx and 0 <= jj < ny and board[ii, jj]:
                q_ = pos[lin(ii, jj)]
                lr_ += [k_, k_]; lc_ += [k_, q_]; lv_ += [1.0, -1.0]
    Lap = sp.csr_matrix((lv_, (lr_, lc_)), shape=(nbd, nbd)) / DX ** 2      # Neumann Laplacian (positive)
    kz = math.pi / (2 * CASE_H)
    nair = 40
    av, aw = sla.eigsh(Lap.tocsc() + 1e-6 * sp.eye(nbd), k=nair, sigma=0, which="LM")
    aw = aw / np.sqrt((aw ** 2).sum(0) * DX * DX / (nbd * DX * DX))        # mean square 1 over the board
    fa = C_AIR * np.sqrt(np.maximum(av, 0) + kz ** 2) / (2 * math.pi)
    Abd = nbd * DX * DX
    Lam = Abd * CASE_H / 2.0                                                 # integral of psi^2 cos^2 over the volume
    full_s = P @ vecs
    Cba = (full_s[bi, :].T @ aw) * DX * DX                                   # (structural modes, air modes)
    ns = len(f)
    Ob = np.diag((2 * math.pi * f) ** 2); Oa = np.diag((2 * math.pi * fa) ** 2)
    Dg = RHO_AIR * C_AIR ** 2 / Lam
    # --- The air under the board: a layer between the board and the floor, open all round
    # the sides (pressure released just past the rim: p = 0 on the ring), with standing waves
    # across it and, from about 280 Hz, up and down between board and floor. The board pushes
    # it the opposite way to the air above.
    ii2 = np.flatnonzero(interior.ravel()[bi])                               # interior nodes, in board numbering
    Ld = Lap.tocsr()[ii2][:, ii2]
    nu_ = 40
    uv, uw = sla.eigsh(Ld.tocsc(), k=nu_, sigma=0, which="LM")
    psiU = np.zeros((nbd, nu_)); psiU[ii2] = uw
    psiU = psiU / np.sqrt((psiU ** 2).sum(0) / nbd)
    fu = C_AIR * np.sqrt(np.maximum(uv, 0)) / (2 * math.pi)
    LamU = Abd * FLOOR_GAP                                                   # the horizontal modes fill the whole layer evenly
    Cbu = (full_s[bi, :].T @ psiU) * DX * DX
    Ou = np.diag((2 * math.pi * fu) ** 2); Du = RHO_AIR * C_AIR ** 2 / LamU
    print(f"air under the board: lowest {fu[:4].round(0)} Hz")
    Z_ = np.zeros
    Kc = np.block([[Ob, Cba, -Cbu], [Z_((nair, ns)), Oa, Z_((nair, nu_))], [Z_((nu_, ns)), Z_((nu_, nair)), Ou]])
    Mc = np.block([[np.eye(ns), Z_((ns, nair)), Z_((ns, nu_))],
                   [-Dg * Cba.T, np.eye(nair), Z_((nair, nu_))],
                   [Du * Cbu.T, Z_((nu_, nair)), np.eye(nu_)]])
    import scipy.linalg as la
    lam, Vc = la.eig(Kc, Mc)
    lam, Vc = lam.real, Vc.real
    o = np.argsort(lam); lam, Vc = lam[o], Vc[:, o]
    keep = (lam > (2 * math.pi * 5.0) ** 2) & (lam < lam_top)
    lam, Vc = lam[keep], Vc[:, keep]
    Q = Vc[:ns]                                                              # structural part
    nq = np.sqrt((Q ** 2).sum(0)); nq[nq == 0] = 1
    # Share of the potential energy held by the air (structure: w^2 q^2 / 2; air: Lam P^2 / (2 rho c^2))
    eS = ((2 * math.pi * f)[:, None] ** 2 * Q ** 2).sum(0)
    eA = (Lam / (RHO_AIR * C_AIR ** 2)) * (Vc[ns:ns + nair] ** 2).sum(0) + (LamU / (RHO_AIR * C_AIR ** 2)) * (Vc[ns + nair:] ** 2).sum(0)
    airE = eA / (eA + eS)
    Q = Q / nq                                                               # unit structural modal mass
    f = np.sqrt(lam) / (2 * math.pi)
    rimShare = (Q ** 2 * rimE[:, None]).sum(0)
    braceShare = (Q ** 2 * brE[:, None]).sum(0)
    print(f"coupled body: {len(f)} modes; lowest {f[:6].round(1)}; air modes from {fa[0]:.0f} Hz; "
          f"rim-heavy modes {(rimShare > 0.3).sum()}, air-heavy {(airE > 0.5).sum()}")
    full = full_s @ Q                                       # (nx*ny, n)

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
    # Close microphones, as on studio recordings like the reference samples: a pair about
    # 35 cm above the strings, one over the bass and tenor, one over the treble, near the
    # hammers. Near-field Rayleigh integral (each board element a small piston) at each mic.
    mics = {"L": (1.05, 0.55, 0.35), "R": (0.40, 0.25, 0.35)}
    micL = np.zeros(len(f)); micR = np.zeros(len(f))
    bnz = board.ravel()
    for name, (mx, my, mz) in mics.items():
        r_ = np.sqrt((xs[bnz] - mx) ** 2 + (ys[bnz] - my) ** 2 + mz ** 2)
        for i in range(len(f)):
            ph_ = np.exp(-1j * k[i] * r_) / r_
            p_ = np.sum(full[bnz, i] * ph_) * k[i]           # x k: radiated pressure rises with frequency like the far field's
            val = abs(p_) * (-1 if p_.real < 0 else 1)
            (micL if name == "L" else micR)[i] = val

    def lit(v):
        t = "%.5g" % v
        if "." not in t and "e" not in t:
            t += ".0"
        return t + "f"

    # Radiation damping: power radiated from both faces, from the far field over the hemisphere
    th = (np.arange(8) + 0.5) / 8 * (math.pi / 2)
    ph = np.arange(24) / 24 * 2 * math.pi
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    dOmega = (math.pi / 2 / 8) * (2 * math.pi / 24) * np.sin(TH).ravel()
    ux = (np.sin(TH) * np.cos(PH)).ravel(); uy = (np.sin(TH) * np.sin(PH)).ravel()
    bn = board.ravel()
    bx, by_ = xs[bn], ys[bn]
    etaRad = np.zeros(len(f))
    for i in range(len(f)):
        w = full[bn, i] * area
        kk = 2 * math.pi * f[i] / 343.0
        W = np.exp(-1j * kk * (np.outer(ux, bx) + np.outer(uy, by_))) @ w
        omega = 2 * math.pi * f[i]
        etaRad[i] = omega * 1.2 / (2 * math.pi ** 2 * 343.0) * np.sum(np.abs(W) ** 2 * dOmega)
    print("radiation loss factor: median %.4f, 10-90%% %.4f-%.4f" % (np.median(etaRad), *np.percentile(etaRad, [10, 90])))

    with open(args.out, "w") as fo:
        fo.write("// Generated by tools/soundboard/build_board.py: do not edit by hand.\n")
        fo.write(f"// Soundboard modes below {F_TOP:.0f} Hz of a tapered, ribbed grand soundboard with bridges,\n// coupled to its rim, iron frame, braces, legs and the air in the case (Yamaha C5 body).\n")
        fo.write("// Per mode: frequency [Hz], shape at the 16 bridge zones (unit modal mass), radiation towards\n")
        fo.write("// the left and right listening points, net volume (both in arbitrary units), and the loss\n")
        fo.write("// factor from the sound it radiates (both faces).\n")
        fo.write(f"static constexpr int kBoardModeCount = {len(f)};\n")
        fo.write(f"static constexpr int kBoardZones = {N_BASS + N_TREBLE};\n")
        fo.write("struct BoardModeData { float f; float phi[%d]; float radLeft, radRight, volume, radiationLoss, rimShare, airShare, micLeft, micRight; };\n" % (N_BASS + N_TREBLE))
        fo.write("static const BoardModeData kBoardModes[kBoardModeCount] = {\n")
        for i in range(len(f)):
            fo.write("    { %.2ff, { %s }, %s, %s, %s, %s, %s, %s, %s, %s },\n" % (
                f[i], ", ".join(lit(v) for v in phi[:, i]), lit(radL[i]), lit(radR[i]), lit(vol[i]), lit(etaRad[i]),
                lit(rimShare[i] + braceShare[i]), lit(airE[i]), lit(micL[i]), lit(micR[i])))
        fo.write("};\n")
    print("wrote", args.out)

    if args.json:
        import json
        yy = np.linspace(0, B, 120)
        outline = [[A, 0.0], [0.0, 0.0]] + [[round(float(bentside(v)), 3), round(float(v), 3)] for v in yy] + [[A, B]]
        sel = [i for i in range(len(f)) if i < 40 or i % 6 == 0]
        shapes = []
        for i in sel:
            w = full[:, i].reshape(nx, ny)
            m = np.abs(w).max()
            shapes.append({"f": round(float(f[i]), 1), "eta": round(float(etaRad[i]), 4),
                           "z": [int(v) for v in np.round(w.T.ravel() / m * 99)]})
        json.dump({"A": A, "B": B, "nx": nx, "ny": ny, "dx": DX, "grain": math.degrees(GRAIN), "outline": outline,
                   "beams": beam_lines, "zones": [[round(float(a_), 3), round(float(b_), 3)] for a_, b_ in zip(zx, zy)],
                   "nBass": N_BASS, "modes": shapes, "allf": [round(float(v), 1) for v in f]},
                  open(args.json, "w"), separators=(",", ":"))
        print("wrote", args.json)

    ff, fphi, frad = frame_modes(np.array(zx), np.arange(N_BASS + N_TREBLE) < N_BASS, 3000.0)
    write_frame(args.frame_out, ff, fphi, frad)

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
