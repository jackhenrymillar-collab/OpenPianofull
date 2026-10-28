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
import argparse, math, os
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as sla
from scipy import ndimage

DX = 0.008                     # grid step [m]
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
# The bridges also twist (rocking across their width, under the strings' lengthwise pull) and
# tilt along their length (under the strings' sideways pull, applied at the bridge top, a lever
# of about its height above the board). Torsion of a square beech section: GJ = G 0.141 a^4.
BRIDGE_G = 1.0e9
TREBLE_W, BASS_W = 0.030, 0.035          # width = height [m]
# The bass bridge stands on an apron, a thin plank glued under its end that carries it out
# towards the rim: it adds mass and ties the bridge's tail end to the rim through the board edge
APRON_M, APRON_K = 0.35, 4.0e5           # [kg], [N/m]
F_TOP = 5000.0
# The body around the board (a Yamaha C5): the rim it is glued into, the braces under it,
# the iron frame bolted on the rim, the legs, and the air in the case above the board.
RIM_EI = 12e9 * 0.035 * 0.26 ** 3 / 12.0 * 3.0   # laminated maple, 26 cm deep; x3: the bolted iron frame ties it
RIM_M = 650.0 * 0.035 * 0.26                      # [kg/m]
FRAME_MASS, FRAME_BOLTS = 125.0, 10               # iron frame [kg], spread over its bolts on the rim
# The board's edge sits glued on the rim liner's ledge, a strip about 2.5 cm wide inside the
# rim. The joint holds the edge level with the rim (the ring nodes are the rim's); against
# rotation it is only as stiff as the ledge strip in compression across the grain, whose
# stress field spreads about one width deep: k_rot = E_perp w^2 / 12 per metre of edge. So the
# thinned edge rocks on its ledge rather than being clamped. Per ring node (spacing DX) the
# hinge is EDGE_K DX^2 = k_rot DX. LEDGE_SCALE (env) scales it for A/B tests.
LEDGE_W = float(os.environ.get("LEDGE_W", "0.025"))                      # ledge width [m]
LINER_E = 0.4e9                                                           # maple liner across the grain [Pa]
EDGE_K = LINER_E * LEDGE_W ** 2 / 12.0 / DX * float(os.environ.get("LEDGE_SCALE", "1"))   # [N/m per node]
TIE_K = 1.0e8                                     # brace ends let into the rim [N/m]
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
# Material losses (loss factor eta = energy lost per radian / energy stored), from measurements:
# spruce along the grain ~0.007, across it and in shear ~0.018 (Haines 1979; Ono & Norimoto 1983),
# beech ~0.012, maple ~0.010, cast iron ~0.002, glue lines ~0.05, castor cups on the floor ~0.25.
ETA = {"along": 0.007, "across": 0.018, "shear": 0.018, "rib": 0.007, "bridge": 0.012, "crown": 0.010,
       "rim": 0.010, "rimtwist": 0.010, "glue": 0.05, "leg": 0.25, "brace": 0.010, "tie": 0.05,
       "frame": 0.002, "bolt": 0.02, "nose": 0.02, "lid": 0.02, "pinblock": 0.05, "str_damped": float(os.environ.get("ETA_SD", "0.3")), "str_free": float(os.environ.get("ETA_SF", "0.02"))}
PLANK_W, PLANK_SPREAD = 0.12, 0.08   # the board is glued from ~12 cm planks; each plank's stiffness varies ~8%
RIM_GJ = 1.1e9 * 0.26 * 0.035 ** 3 / 3.0 * 0.92 * 3.0   # rim twisting stiffness (maple, x3 for the bolted frame) [N m^2]
RIM_IP = 650.0 * 0.035 * 0.26 * (0.035 ** 2 + 0.26 ** 2) / 12.0   # rim rotary inertia per length [kg m]
PIN_OFFSET = 0.008          # strings bear on the bridge's front pins, ~8 mm on the speaking-length side of its centre line
LISTEN = {"L": (1.40, -1.20, 1.00), "R": (0.00, -1.20, 1.00)}     # the player's seat, a stereo pair a little apart
MICS = {"L": (1.05, 0.55, 0.35), "R": (0.40, 0.25, 0.35)}         # close pair over the strings, near the hammers
LID_ANGLE = math.radians(42.0)   # full stick
RIM_TOP = 0.24                   # rim top above the board [m]; the lid rests there along the spine
LID_HINGE_X = A + 0.02           # the hinge line, just outside the spine
LID_GAP_MIN = 0.10               # the rim-to-lid gap along the bentside counts as open above this [m]
LID_FRONT = float(os.environ.get("LID_FRONT", "0.40"))   # the main lid's front edge: the front flap is folded back on full stick [m from the belly rail]
LID_TAIL = B + 0.10              # the lid's tail overhang
LID_REFL = 0.9                   # lacquered hardwood
LID_H, LID_E, LID_RHO, LID_NU = 0.019, 8e9, 550.0, 0.3   # the lid panel: veneered softwood core, 19 mm
LID_PROP = (0.55, 0.12)          # where the prop stick meets the lid: y [m], and this far in from the outer edge [m]
LID_PROP_K = 1.0e7               # the stick in compression [N/m]
FLOOR_REFL = 0.8
N_BASS, N_TREBLE = 4, 28       # bridge zones (bass bridge keys A0..E2, treble bridge F2..C8)


# Measured geometry (optional, --geometry FILE.json). Anything the file leaves out keeps the
# generic shape below. Coordinates are metres on the board: x across the strings (0 at the
# treble edge .. A at the spine), y along them (0 at the belly rail .. B at the tail).
#   "A", "B"               bounding box
#   "outline"              polygon [[x, y], ...] of the vibrating board (inside the rim, in front of the cut-off bar)
#   "ribs"                 [[[x0, y0], [x1, y1]], ...] rib centre lines, end to end
#   "treble_bridge"        polyline [[x, y], ...] from the C8 end to the F2 end
#   "bass_bridge"          polyline from the E2 end to the A0 end
#   "cutoff"               [[x0, y0], [x1, y1]] the cut-off bar across the bass corner
#   "grain_deg", "h_mid", "h_edge", "rib_w", "rib_t"
GEOM = {}


def load_geometry(path):
    import json
    g = json.load(open(path))
    GEOM.update(g)
    G = globals()
    for key, name in (("A", "A"), ("B", "B"), ("h_mid", "H_MID"), ("h_edge", "H_EDGE"), ("rib_w", "RIB_W"), ("rib_t", "RIB_T")):
        if key in g: G[name] = float(g[key])
    if "grain_deg" in g: G["GRAIN"] = math.radians(float(g["grain_deg"]))
    if "rib_w" in g or "rib_t" in g:
        G["RIB_EI"] = 2.0 * E_L * RIB_W * RIB_T ** 3 / 12.0
        G["RIB_M"] = RHO * RIB_W * RIB_T
    if "ribs" in g and len(g["ribs"]) > 1:
        # mean spacing between neighbouring ribs, for the membrane stiffness along them
        mids = np.array([[(r[0][0] + r[1][0]) / 2, (r[0][1] + r[1][1]) / 2] for r in g["ribs"]])
        d = np.sort(np.hypot(*(mids[1:] - mids[:-1]).T))
        G["RIB_SPACING"] = float(np.median(d))
    print(f"geometry from {path}: " + ", ".join(k for k in g))


def _polyline_fn(points):
    """A curve parametrised by arc length, t = 0 at the first point .. 1 at the last."""
    pts = np.asarray(points, float)
    seg = np.hypot(*np.diff(pts, axis=0).T)
    s = np.concatenate([[0.0], np.cumsum(seg)]) / seg.sum()
    return lambda t: (np.interp(t, s, pts[:, 0]), np.interp(t, s, pts[:, 1]))


def bentside(y):
    """Treble-side edge of the board: straight near the keyboard, then curving in to the tail."""
    if "outline" in GEOM:
        # the smallest x where the outline crosses this y
        pts = np.asarray(GEOM["outline"], float)
        y = np.asarray(y, float)
        out = np.full(y.shape, A, float)
        for (x0, y0), (x1, y1) in zip(pts, np.roll(pts, -1, axis=0)):
            if y0 == y1:
                continue
            u = (y - min(y0, y1)) / abs(y1 - y0)
            hit = (u >= 0) & (u <= 1)
            xc = x0 + (x1 - x0) * (y - y0) / (y1 - y0)
            out = np.where(hit, np.minimum(out, xc), out)
        return out
    t = np.clip((y - 0.25) / (B - 0.25), 0.0, 1.0)
    return 1.10 * t ** 1.8


def inside(x, y):
    if "outline" in GEOM:
        from matplotlib.path import Path
        x, y = np.asarray(x, float), np.asarray(y, float)
        return Path(np.asarray(GEOM["outline"], float)).contains_points(np.column_stack([x.ravel(), y.ravel()])).reshape(x.shape)
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
    if "treble_bridge" in GEOM:
        return _polyline_fn(GEOM["treble_bridge"])(t)
    k = 88 - t * (88 - 21)
    L = np.interp(k, np.arange(1, 89), LEN)
    x = 0.12 + 0.80 * t
    return x, CAPO_Y + FAN * L


def bass_bridge(t):
    """t = 0 at E2 .. 1 at A0: the cross-strung bass runs at about 25 degrees over the tenor,
    towards the spine and the tail."""
    if "bass_bridge" in GEOM:
        return _polyline_fn(GEOM["bass_bridge"])(t)
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
# move the strings' far ends; the shapes at the key zones along the front bar are what
# the engine needs.
IRON_E, IRON_RHO = 1.1e11, 7200.0
FRAME_A, FRAME_I = 0.0022, 2.6e-6          # bar section [m^2], bending [m^4] (8 cm deep ribbed bars)
FRONT_Y = -0.05                             # the front bar, just ahead of the belly rail
BOLT_K = 2.0e7                              # each bolt into the rim [N/m]
PINBLOCK_K = float(os.environ.get("PINBLOCK_K", "5e7"))   # the pin block under each front-bar node [N/m]


def frame_bars():
    """The frame's bar network: stiffness and mass matrices, node coordinates, bolt nodes and
    the front-bar nodes (where the strings end), ordered across."""
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
    front = sorted([k for k in range(n) if abs(C[k, 1] - FRONT_Y) < 0.01], key=lambda k: C[k, 0])
    return K, Mv, C, bolts, front


def lid_plate(f_top):
    """The lid: a 19 mm veneered panel (isotropic here), hinged along the spine (a pinned line),
    held up by the prop stick near its front bentside corner, free along its other edges. Its
    bending modes are solved as a plate on a grid in the lid's own plane (s up from the hinge,
    y along the piano). Returns f, shapes at the nodes (unit modal mass), the nodes' positions
    in space and over the board plane, the node area and each node's distance to the lid's edge."""
    d = 0.03
    ys = np.arange(LID_FRONT, LID_TAIL + 1e-9, d)
    smax = (LID_HINGE_X - (np.asarray(bentside(np.clip(ys, 0, B)), float) - 0.04)) / math.cos(LID_ANGLE)
    ss = np.arange(0.0, smax.max() + 1e-9, d)
    nsg, nyg = len(ss), len(ys)
    mask = ss[:, None] <= smax[None, :]                                     # (s, y)
    inner = ndimage.binary_erosion(mask, structure=np.ones((3, 3)), iterations=1)
    lin_ = lambda i, j: i * nyg + j
    def op(stencil):
        r, c, v = [], [], []
        pts = np.argwhere(inner)
        for k, (i, j) in enumerate(pts):
            for (di, dj, w) in stencil:
                ii, jj = i + di, j + dj
                if 0 <= ii < nsg and 0 <= jj < nyg:
                    r.append(k); c.append(lin_(ii, jj)); v.append(w / d ** 2)
        return sp.csr_matrix((v, (r, c)), shape=(len(pts), nsg * nyg))
    Css = op([(-1, 0, 1), (0, 0, -2), (1, 0, 1)])
    Cyy = op([(0, -1, 1), (0, 0, -2), (0, 1, 1)])
    Csy = op([(0, 0, 1), (1, 1, 1), (1, 0, -1), (0, 1, -1)])
    Dl = LID_E * LID_H ** 3 / (12.0 * (1.0 - LID_NU ** 2))
    K = (Css.T @ Css + Cyy.T @ Cyy + LID_NU * (Css.T @ Cyy + Cyy.T @ Css) + 2.0 * (1.0 - LID_NU) * Csy.T @ Csy) * (Dl * d * d)
    # The prop: a stiff point support
    yp, inset = LID_PROP
    jp = int(np.argmin(np.abs(ys - yp))); ip = int(np.argmin(np.abs(ss - (smax[jp] - inset))))
    K = K + sp.csr_matrix(([LID_PROP_K], ([lin_(ip, jp)], [lin_(ip, jp)])), shape=(nsg * nyg, nsg * nyg))
    # Free nodes: on the lid, not on the hinge line (pinned there)
    free = mask.copy(); free[0, :] = False
    fi = np.flatnonzero(free.ravel())
    m = LID_RHO * LID_H * d * d
    Kf = K.tocsr()[fi][:, fi].toarray() / m
    vals, vecs = np.linalg.eigh(Kf)
    keep = (vals > (2 * math.pi * 5.0) ** 2) & (vals < (2 * math.pi * f_top) ** 2)
    f = np.sqrt(vals[keep]) / (2 * math.pi)
    V = np.zeros((nsg * nyg, keep.sum())); V[fi] = vecs[:, keep] / math.sqrt(m)   # unit modal mass
    on = np.flatnonzero(mask.ravel())
    S_, Y_ = np.meshgrid(ss, ys, indexing="ij")
    sN, yN = S_.ravel()[on], Y_.ravel()[on]
    hinge = np.array([LID_HINGE_X, 0.0, RIM_TOP])
    tl = np.array([-math.cos(LID_ANGLE), 0.0, math.sin(LID_ANGLE)])
    P3 = hinge[None, :] + sN[:, None] * tl[None, :] + np.stack([np.zeros_like(yN), yN, np.zeros_like(yN)], 1)
    dE = (ndimage.distance_transform_edt(mask) * d).ravel()[on]
    print(f"lid: {len(on)} nodes, {LID_RHO * LID_H * len(on) * d * d:.1f} kg, {len(f)} modes below {f_top:.0f} Hz, lowest {f[:5].round(1)}; "
          f"prop at ({P3[np.argmin(np.hypot(sN - ss[ip], yN - ys[jp])), 0]:.2f}, {ys[jp]:.2f})")
    return f, V[on], P3, np.stack([sN, yN], 1), d * d, dE


def frame_modes(zone_x, zone_is_bass, f_top):
    """The frame alone on rigid bolts (kept for reference; the body model solves it with the rim)."""
    K, Mv, C, bolts, front = frame_bars()
    n = len(Mv); coords = C
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


def write_frame(path, f, phi, rad, eta=None):
    if eta is None:
        eta = np.full(len(f), 0.002)
    def lit(v):
        t = "%.5g" % v
        if "." not in t and "e" not in t:
            t += ".0"
        return t + "f"
    with open(path, "w") as fo:
        fo.write("// Generated by tools/soundboard/build_board.py: do not edit by hand.\n")
        fo.write("// Bending modes of the cast-iron frame (bars bolted on the rim). Per mode: frequency [Hz], shape at the\n")
        fo.write("// strings' far ends for each key zone (unit modal mass, 1/sqrt(kg)), its sound at the player's seat\n")
        fo.write("// (relative units, like the board's) and its structural loss factor. Solved with the rim and board.\n")
        fo.write(f"static constexpr int kFrameModeCount = {len(f)};\n")
        fo.write("struct FrameModeData { float f; float phi[%d]; float rad, loss; };\n" % phi.shape[1])
        fo.write("static const FrameModeData kFrameModes[kFrameModeCount] = {\n")
        for i in range(len(f)):
            fo.write("    { %.2ff, { %s }, %s, %s },\n" % (f[i], ", ".join(lit(v) for v in phi[i]), lit(rad[i]), lit(eta[i])))
        fo.write("};\n")
    print("wrote", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="engine/src/piano/BoardModes.inc")
    ap.add_argument("--frame-out", default="engine/src/piano/FrameModes.inc")
    ap.add_argument("--plot")
    ap.add_argument("--json", help="also write the board geometry and some mode shapes for the visualiser")
    ap.add_argument("--geometry", help="measured board geometry (JSON, see GEOM above); anything missing keeps the generic shape")
    args = ap.parse_args()
    if args.geometry:
        load_geometry(args.geometry)

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
    # Planks: strips along the grain, each a little stiffer or softer than its neighbours
    c0, s0 = math.cos(GRAIN), math.sin(GRAIN)
    plank = np.floor((-s0 * X + c0 * Y) / PLANK_W).astype(int)
    prng = np.random.default_rng(5)
    pf = 1.0 + PLANK_SPREAD * prng.standard_normal(plank.max() - plank.min() + 1).clip(-2, 2)
    pfac = pf[plank - plank.min()]
    D11 = E_L * pfac * h ** 3 / 12.0          # along the grain
    D22 = E_T * pfac * h ** 3 / 12.0          # across it
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
    # Strain energy is kept per part (wood along / across the grain, shear, ribs, bridges, rim,
    # glue, legs ...), so each mode's loss can be worked out from what it actually bends
    parts = {}
    def part(name, mat):
        parts[name] = parts[name] + mat if name in parts else mat
        return mat
    K = (part("along", C1.T @ g(D11) @ C1 + 0.5 * (C1.T @ g(D12) @ C2 + C2.T @ g(D12) @ C1))
         + part("across", C2.T @ g(D22) @ C2 + 0.5 * (C1.T @ g(D12) @ C2 + C2.T @ g(D12) @ C1))
         + part("shear", C3.T @ g(D66) @ C3))
    M = sp.diags(mass.ravel())

    # Bridges: beams along their curves, sampled every DX, bilinear onto the grid
    def interp_rows(xs, ys):
        r, c, v = [], [], []
        xs = np.clip(xs, 0.0, (nx - 2) * DX); ys = np.clip(ys, 0.0, (ny - 2) * DX)   # keep measured beams on the grid
        for k, (x, y) in enumerate(zip(xs, ys)):
            fx, fy = x / DX, y / DX
            i0, j0 = int(fx), int(fy)
            ax, ay = fx - i0, fy - j0
            for di, dj, w in ((0, 0, (1 - ax) * (1 - ay)), (1, 0, ax * (1 - ay)), (0, 1, (1 - ax) * ay), (1, 1, ax * ay)):
                r.append(k); c.append(lin(i0 + di, j0 + dj)); v.append(w)
        return sp.csr_matrix((v, (r, c)), shape=(len(xs), nx * ny))
    beams = [(treble_bridge, TREBLE_EI, TREBLE_M), (bass_bridge, BASS_EI, BASS_M)]
    # Ribs: straight beams across the grain, every RIB_SPACING along it, cut to the board
    # (or the measured rib lines, when a geometry file gives them)
    gdir = np.array([c, s_]); rdir = np.array([-s_, c])
    centre = np.array([0.85, 0.75])
    no_taper = set()
    for r_ in GEOM.get("ribs", []):
        beams.append((_polyline_fn(r_), RIB_EI, RIB_M))
    for off in (np.arange(-1.2, 1.2, RIB_SPACING) if "ribs" not in GEOM else []):
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
    if "cutoff" in GEOM:
        # The cut-off bar: a heavy spruce beam across the bass corner (about 25 x 100 mm), glued
        # under the board; the corner behind it is left out of the vibrating area by the outline
        no_taper.add(len(beams))
        beams.append((_polyline_fn(GEOM["cutoff"]), 2.0 * E_L * 0.025 * 0.10 ** 3 / 12.0, RHO * 0.025 * 0.10))
    beam_lines = []
    for fn, EI, mlin in beams:
        u = np.linspace(0, 1, 40)
        bxs, bys = fn(u)
        beam_lines.append([[round(float(a_), 3), round(float(b_), 3)] for a_, b_ in zip(bxs, bys)])
    for bi_, (fn, EI, mlin) in enumerate(beams):
        _, _, L = curve_points(fn, 2)
        n = max(3, int(L / DX) + 1)
        xs, ys, L = curve_points(fn, n)
        ds = L / (n - 1)
        Bm = interp_rows(xs, ys)
        D2 = sp.diags([np.ones(n - 2), -2 * np.ones(n - 2), np.ones(n - 2)], [0, 1, 2], shape=(n - 2, n)) / ds ** 2
        C = D2 @ Bm
        if bi_ >= 2:
            # Ribs are planed down towards their ends, where they run into the rim (about 25 mm
            # tall in the middle, ~18 mm at the ends: EI falls to ~40% over the last 15 cm)
            sc = np.linspace(0, L, n)[1:-1]
            dend = np.minimum(sc, L - sc)
            if not os.environ.get('NO_TAPER') and bi_ not in no_taper: C = sp.diags(np.sqrt(0.4 + 0.6 * np.clip(dend / 0.15, 0.0, 1.0))) @ C
        K = K + part("bridge" if bi_ < 2 else "rib", C.T @ C * (EI * ds))
        M = M + Bm.T @ sp.diags(np.full(n, mlin * ds)) @ Bm
        if bi_ < 2:
            # Twist: the slope across the bridge's width, its rate of change along it (stiffness)
            # and its rotary inertia
            w_ = TREBLE_W if bi_ == 0 else BASS_W
            Tw = (interp_rows(xs, ys + 0.5 * w_) - interp_rows(xs, ys - 0.5 * w_)) / w_
            D1 = sp.diags([-np.ones(n - 1), np.ones(n - 1)], [0, 1], shape=(n - 1, n)) / ds
            K = K + part("bridge", (D1 @ Tw).T @ (D1 @ Tw) * (BRIDGE_G * 0.141 * w_ ** 4 * ds))
            M = M + Tw.T @ Tw * (750.0 * w_ * w_ * (2 * w_ * w_) / 12.0 * ds)
            if bi_ == 1:
                # The apron at the bass bridge's tail end (A0, t = 1)
                tail = interp_rows(xs[-1:], ys[-1:])
                M = M + tail.T @ tail * APRON_M
                K = K + part("bridge", tail.T @ tail * APRON_K)

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
    K = K + part("crown", sp.csr_matrix(((Eh_r[ipts] * kappa ** 2) * area, (np.array([lin(i, j) for i, j in pts]),) * 2), shape=(nx * ny, nx * ny)))
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
    K = K + part("rim", Cr.T @ Cr * (RIM_EI * ds))
    rimMass = np.full(nr, RIM_M * ds)                      # (the iron frame is modelled on its own below)
    M = M + R.T @ sp.diags(rimMass) @ R
    # The glue joint: the board's edge is glued to the rim's inner face, so its slope there
    # follows the rim's twist (a stiff hinge spring between the edge slope and the twist)
    ring_index = {(int(i), int(j)): r_ for r_, (i, j) in enumerate(rp)}
    glue = []                                               # (ring node, inner neighbour, ring index)
    for (i, j) in np.argwhere(ring):
        nb = [(i + a, j + b) for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)) if 0 <= i + a < nx and 0 <= j + b < ny and interior[i + a, j + b]]
        for (ii, jj) in nb[:1]:
            glue.append((lin(i, j), lin(ii, jj), ring_index[(int(i), int(j))]))
    # Legs: two at the keyboard corners, one under the tail
    def nearest_ring(x, y):
        d = np.hypot(rp[:, 0] * DX - x, rp[:, 1] * DX - y)
        return int(np.argmin(d))
    legs = [nearest_ring(A, 0.0), nearest_ring(0.0, 0.0), nearest_ring(0.75, B)]
    for r_ in legs:
        K = K + part("leg", sp.csr_matrix(([LEG_K], ([rlin[r_]], [rlin[r_]])), shape=(nx * ny, nx * ny)))

    # --- Back braces: beams under the board (not touching it) from a hub towards the tail out
    # to the rim, each with its own nodes, the ends pinned onto the rim and the hub shared
    hub = np.array([0.70, 1.05])
    ends = [(A, 0.35), (A, 0.95), (0.95, B), (0.35, 1.05), (0.12, 0.45)]
    nb_nodes = 12
    nbr = len(ends) * nb_nodes + 1                          # + the shared hub node
    Kb = sp.lil_matrix((nbr, nbr)); Mb = np.zeros(nbr)
    tie = []                                                # (brace node, ring index)
    Cbr = np.zeros((nbr, 2)); Cbr[nbr - 1] = hub
    for e, (ex, ey) in enumerate(ends):
        L = float(np.hypot(ex - hub[0], ey - hub[1]))
        h_ = L / nb_nodes
        nodes = [nbr - 1] + [e * nb_nodes + k for k in range(nb_nodes)]   # hub .. rim end
        for k, n_ in enumerate(nodes):
            Cbr[n_] = hub + (np.array([ex, ey]) - hub) * k / nb_nodes
        for k in range(1, nb_nodes):
            a, b, c_ = nodes[k - 1], nodes[k], nodes[k + 1]
            v = np.zeros(nbr); v[a] += 1; v[b] -= 2; v[c_] += 1
            Kb += sp.lil_matrix(np.outer(v, v) * BRACE_EI / h_ ** 3)
        for k in nodes:
            Mb[k] += BRACE_M * h_
        tie.append((nodes[-1], nearest_ring(ex, ey)))
    Kb = Kb.tocsr()
    # --- The iron frame, bolted onto the rim (its bars and bolts, see frame_bars)
    Kfr, Mfr, Cfr, bolts, frontNodes = frame_bars()
    nfr = len(Mfr)
    # Global layout: [board nodes | brace nodes | rim twist (one per ring node) | frame nodes]
    oB, oT, oF = nf, nf + nbr, nf + nbr + nr
    NT = oF + nfr
    def emb(mat, at, size):
        mat = sp.coo_matrix(mat)
        return sp.csr_matrix((mat.data, (mat.row + at, mat.col + at)), shape=(NT, NT))
    def spring_rows(rows_):
        r, c, v = [], [], []
        for k_, terms in enumerate(rows_):
            for (col, val) in terms:
                r.append(k_); c.append(col); v.append(val)
        return sp.csr_matrix((v, (r, c)), shape=(len(rows_), NT))
    pos_of = -np.ones(nx * ny, int); pos_of[np.flatnonzero(free.ravel())] = np.arange(nf)
    gparts = {name: emb(P.T @ mat @ P, 0, nf) for name, mat in parts.items()}
    gparts["brace"] = emb(Kb, oB, nbr)
    gparts["frame"] = emb(Kfr, oF, nfr)
    # Rim twist: GJ between neighbouring ring positions (closed loop)
    tw = spring_rows([[(oT + r_, -1.0 / ds), (oT + (r_ + 1) % nr, 1.0 / ds)] for r_ in range(nr)])
    gparts["rimtwist"] = tw.T @ tw * (RIM_GJ * ds)
    # Glue: edge slope (inner neighbour minus ring node, over DX) follows the rim's twist
    gl = spring_rows([[(pos_of[a], -1.0 / DX), (pos_of[b], 1.0 / DX), (oT + r_, -1.0)] for a, b, r_ in glue])
    gparts["glue"] = gl.T @ gl * (EDGE_K * DX * DX)
    # Brace ends pinned onto the rim
    ti = spring_rows([[(pos_of[rlin[r_]], 1.0), (oB + bn_, -1.0)] for bn_, r_ in tie])
    gparts["tie"] = ti.T @ ti * TIE_K
    # Frame bolts onto the nearest point of the rim
    bo = spring_rows([[(pos_of[rlin[nearest_ring(*Cfr[b_])]], 1.0), (oF + int(b_), -1.0)] for b_ in bolts])
    gparts["bolt"] = bo.T @ bo * BOLT_K
    # Nose bolts: the plate is also held down through the board onto the back braces, by two
    # or three bolts through its web where the struts pass over a brace (a spacer sets the
    # plate's height; the board has a clearance hole and is not touched). They tie the plate's
    # middle, where the struts meet, to the stiff beams instead of leaving it hanging between
    # the rim bolts. NO_NOSE (env) leaves them out for A/B tests.
    nose = []
    if not os.environ.get("NO_NOSE"):
        inner = [k_ for k_ in range(nfr) if Cfr[k_, 1] > 0.0 and Cfr[k_, 1] < 1.35 and Cfr[k_, 0] < 1.30
                 and abs(Cfr[k_, 0] - (float(bentside(Cfr[k_, 1])) + 0.04)) > 0.05]
        cands = []
        for e in range(len(ends)):
            bn = [e * nb_nodes + k for k in range(2, nb_nodes)]
            d = np.hypot(Cfr[inner][:, None, 0] - Cbr[bn][None, :, 0], Cfr[inner][:, None, 1] - Cbr[bn][None, :, 1])
            i_, j_ = np.unravel_index(int(np.argmin(d)), d.shape)
            cands.append((float(d[i_, j_]), inner[i_], bn[j_]))
        for d_, fn_, bn_ in sorted(cands):
            if d_ < 0.06 and len(nose) < 3 and all(np.hypot(*(Cfr[fn_] - Cfr[f2])) > 0.25 for f2, _ in nose):
                nose.append((fn_, bn_))
        if nose:
            ns = spring_rows([[(oF + int(fn_), 1.0), (oB + int(bn_), -1.0)] for fn_, bn_ in nose])
            gparts["nose"] = ns.T @ ns * BOLT_K
        print("nose bolts at " + ", ".join(f"({Cfr[fn_, 0]:.2f}, {Cfr[fn_, 1]:.2f})" for fn_, _ in nose))
    # The pin block: a laminated maple block, glued to the stretcher across the whole front,
    # that the plate's front flange is screwed down onto every few cm. It holds the front bar
    # (and the strings' far ends on it) along its whole length, not just at the end bolts: a
    # stiff, lossy wooden support (screw joints and wood, eta ~0.05)
    if not os.environ.get('NO_PINBLOCK'):
        pb = spring_rows([[(pos_of[rlin[nearest_ring(*Cfr[fn])]], 1.0), (oF + int(fn), -1.0)] for fn in frontNodes])
        gparts["pinblock"] = pb.T @ pb * PINBLOCK_K
    # The strings themselves. Every string under tension T is a spring between the bridge and
    # the iron frame: T / L along its speaking length to the agraffe or capo bar on the front bar,
    # and T / Lb along its short back length to the hitch pin on the frame behind the bridge.
    # Together (~230 strings, ~3.6 MN/m) they stiffen the bridges about as much as the board
    # does at low frequencies, lifting the lowest modes and flattening their peaks (Conklin
    # 1996), and they tie the board to the frame. Dampers rest on the speaking lengths up to F6
    # and felt mutes the bass and tenor back lengths: those parts of the load are lossy.
    frontC = Cfr[frontNodes]
    def frame_node_near(x, y, cands=None):
        idx = np.arange(len(Cfr)) if cands is None else np.array(cands)
        d = np.hypot(Cfr[idx, 0] - x, Cfr[idx, 1] - y)
        return int(idx[np.argmin(d)]), float(d.min())
    rows_d, rows_f, kd, kf_ = [], [], [], []
    n_rim = 0
    for key in range(1, 89):
        ns = 1 if key <= 8 else 2 if key <= 20 else 3
        T = 1600.0 if key <= 20 else 1100.0 if key <= 30 else 750.0
        Lb = 0.14 if key <= 20 else 0.10 if key <= 40 else 0.06
        Ls = LEN[key - 1]
        if os.environ.get('REAL_T'):
            # tensions as in ScaleDesign (engine): bass 1800->1420 N, tenor 1150->1030, plain ~640
            T = 1800.0 - 380.0 * (key - 1) / 19 if key <= 20 else (1150.0 - 120.0 * (key - 21) / 7 if key <= 28 else 640.0)
            # back length bridge->hitch grows with the speaking length (rear duplex ~1/5 L)
            Lb = float(np.clip(0.2 * Ls, 0.03, 0.25))
        if key <= 20:
            bx_, by_ = bass_bridge((20 - key) / 19.0); ang = math.radians(25.0)
        else:
            bx_, by_ = treble_bridge((88 - key) / (88 - 21)); ang = 0.0
        bx_, by_ = float(bx_), float(by_)
        Br = interp_rows(np.array([bx_]), np.array([by_]))
        bterms = [(int(pos_of[c_]), float(v_)) for c_, v_ in zip(Br.indices, Br.data) if pos_of[c_] >= 0]
        # speaking length: to the front bar, along the string
        fn_, _ = frame_node_near(bx_ - Ls * math.sin(ang), FRONT_Y, frontNodes)
        # back length: to the frame behind the bridge, or the rim if no bar is near
        hx, hy = bx_ + Lb * math.sin(ang), by_ + Lb * math.cos(ang)
        hn_, hd = frame_node_near(hx, hy)
        hterm = [(oF + hn_, -1.0)] if hd < 0.15 else [(int(pos_of[rlin[nearest_ring(hx, hy)]]), -1.0)]
        if hd >= 0.15: n_rim += 1
        if hd >= 0.15 and os.environ.get('NO_RIMBACK'): hterm = [(oF + hn_, -1.0)]
        sp_row = bterms + [(oF + fn_, -1.0)]
        bk_row = bterms + hterm
        damped_speak = key <= 69
        felted_back = key <= 40
        (rows_d if damped_speak else rows_f).append(sp_row); (kd if damped_speak else kf_).append(ns * T / Ls)
        (rows_d if felted_back else rows_f).append(bk_row); (kd if felted_back else kf_).append(ns * T / Lb)
    Sd = spring_rows(rows_d); Sf = spring_rows(rows_f)
    gparts["str_damped"] = Sd.T @ sp.diags(np.array(kd)) @ Sd
    gparts["str_free"] = Sf.T @ sp.diags(np.array(kf_)) @ Sf
    print(f"back lengths to rim: {n_rim}; string load: {(sum(kd) + sum(kf_)) / 1e6:.2f} MN/m ({sum(kd) / 1e6:.2f} through dampers or felt)")
    Kf = sum(gparts.values())
    # The downbearing's compression along the ribs softens bending (stored energy, no loss of its own)
    Kf = Kf + emb(P.T @ (-Gr.T @ sp.diags(np.full(len(pts), Nc * area)) @ Gr) @ P, 0, nf)
    Mf = (emb(P.T @ M @ P, 0, nf) + emb(sp.diags(Mb), oB, nbr) + emb(sp.diags(np.full(nr, RIM_IP * ds)), oT, nr)
          + emb(sp.diags(Mfr), oF, nfr))
    Kf = Kf.tocsc(); Mf = Mf.tocsc()
    rimNodes = pos_of[rlin]
    lam_top = (2 * math.pi * F_TOP) ** 2
    nev = 2600
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
    # Energy shares of the structural modes (kinetic): rim, braces, frame
    rimE = np.array([np.sum(rimMass * vecs[rimNodes, i] ** 2) for i in range(vecs.shape[1])]) + RIM_IP * ds * (vecs[oT:oF] ** 2).sum(0)
    brE = (Mb[:, None] * vecs[oB:oT, :] ** 2).sum(0)
    frE = (Mfr[:, None] * vecs[oF:, :] ** 2).sum(0)
    # Structural loss of each mode: every part's loss factor, weighted by the strain energy it holds
    estr = {name: np.einsum("ij,ij->j", vecs, mat @ vecs) for name, mat in gparts.items()}
    etot = sum(np.maximum(e, 0) for e in estr.values())
    etaS = sum(ETA[name] * np.maximum(e, 0) for name, e in estr.items()) / etot
    print("structural loss factor: median %.4f, 10-90%% %.4f-%.4f; frame-heavy modes %d" %
          (np.median(etaS), *np.percentile(etaS, [10, 90]), (frE > 0.5).sum()))
    # The frame's shapes where the strings end (front bar), kept for the frame table
    frontX = Cfr[frontNodes, 0]
    vfront = vecs[oF + np.array(frontNodes), :]
    vtwist = vecs[oT:oF, :]
    vecs = vecs[:nf]

    # --- The air in the case above the board, coupled to the structure. With the lid on its
    # full stick the space over the board is a wedge: the board below, the rim walls round it,
    # the lid sloping up from its hinge on the spine (rim height there) to about 1.5 m over
    # the treble corner. It is open across the front (above the rim, the flap folded back and
    # the music desk apart) and through the gap between the rim and the lid all along the
    # bentside, closed at the spine and the tail where the lid sits on the rim. Pressure modes
    # of that wedge, with the layer's height in the equation (shallow-layer form of the wave
    # equation, d/dx (h dp/dx) + ... + k^2 h p = 0) and the pressure released on the open faces:
    # two families, one even through the height and one with a half wave between board and
    # lid (cos(pi z / h), whose local cut-on is c / 2h). The board pumps them and they push
    # back on it (vibro-acoustic coupling, solved together).
    bi = np.flatnonzero(board.ravel())
    nbd = len(bi)
    pos = -np.ones(nx * ny, int); pos[bi] = np.arange(nbd)
    xb_, yb_ = X.ravel()[bi], Y.ravel()[bi]
    hLid = RIM_TOP + np.maximum(0.0, LID_HINGE_X - xb_) * math.tan(LID_ANGLE)   # lid height over each node
    ringB = ring.ravel()[bi]
    gap = hLid - RIM_TOP                                                     # lid above the rim top
    openB = ringB & ((yb_ < 0.5 * DX) | ((xb_ < A - 0.06) & (gap > LID_GAP_MIN)))
    boxAir = bool(os.environ.get("BOX_AIR"))
    if boxAir:
        hLid = np.full(nbd, CASE_H); openB = np.zeros(nbd, bool)
    # The openings are not pressure-release planes at the rim: the air just outside them moves
    # with the mode (the end correction, which lowers the resonances) and only beyond that is
    # the pressure released. The air grid is the board grid padded outwards, and each open
    # edge node grows a tongue of outside air END_CORR times the opening's local height,
    # at that height, released at its far end.
    END_CORR = float(os.environ.get("END_CORR", "0.6"))
    pad = int(math.ceil(max(END_CORR, 0.0) * max(gap.max(), 0.01) / DX)) + 2 if not boxAir else 1
    nxa, nya = nx + 2 * pad, ny + 2 * pad
    maskA = np.zeros((nxa, nya), bool); maskA[pad:pad + nx, pad:pad + ny] = board
    hA = np.zeros((nxa, nya)); hA[pad:pad + nx, pad:pad + ny].flat[bi] = hLid
    ext = np.zeros((nxa, nya), bool)
    bij = np.argwhere(board)
    ri_, rj_ = np.meshgrid(np.arange(-pad, pad + 1), np.arange(-pad, pad + 1), indexing="ij")
    rr2 = ri_ ** 2 + rj_ ** 2
    for k_ in np.flatnonzero(openB):
        i, j = bij[k_] + pad
        rad = END_CORR * gap[k_] / DX
        disc = rr2 <= rad * rad
        sl = (slice(i - pad, i + pad + 1), slice(j - pad, j + pad + 1))
        grow = disc & ~maskA[sl]
        ext[sl] |= grow
        hA[sl] = np.where(grow, np.maximum(hA[sl], gap[k_]), hA[sl])
    maskA |= ext
    print(f"case air: {int(openB.sum())} open edge nodes, {int(ext.sum())} outside-air nodes for the end correction")
    aidx = np.flatnonzero(maskA.ravel()); na = len(aidx)
    posA = -np.ones(nxa * nya, int); posA[aidx] = np.arange(na)
    hAv = hA.ravel()[aidx]
    # Board nodes in air numbering
    bInA = posA[(bij[:, 0] + pad) * nya + (bij[:, 1] + pad)]
    # Pressure released at the outer edge of the outside air
    extv = ext.ravel()[aidx]
    dirich = np.zeros(na, bool)
    if END_CORR <= 0.0 and not boxAir:                                       # no end correction: released at the rim itself
        dirich[bInA[openB]] = True
    lr_, lc_, lv_ = [], [], []
    ur_, uc_, uv_ = [], [], []                                               # the plain Laplacian on the board, for the air below
    aij = np.argwhere(maskA)
    for k_, (i, j) in enumerate(aij):
        for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + a, j + b
            if 0 <= ii < nxa and 0 <= jj < nya and maskA[ii, jj]:
                q_ = posA[ii * nya + jj]
                he = 0.5 * (hAv[k_] + hAv[q_])
                lr_ += [k_, k_]; lc_ += [k_, q_]; lv_ += [he, -he]
            elif extv[k_] and END_CORR > 0.0:
                dirich[k_] = True
    LapA = sp.csr_matrix((lv_, (lr_, lc_)), shape=(na, na)) / DX ** 2       # -div(h grad), Neumann (positive)
    # The plain Laplacian on the board nodes (the layer under the board has one height)
    for k_, (i, j) in enumerate(bij):
        for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            ii, jj = i + a, j + b
            if 0 <= ii < nx and 0 <= jj < ny and board[ii, jj]:
                q_ = pos[lin(ii, jj)]
                ur_ += [k_, k_]; uc_ += [k_, q_]; uv_ += [1.0, -1.0]
    Lap = sp.csr_matrix((uv_, (ur_, uc_)), shape=(nbd, nbd)) / DX ** 2
    keepN = np.flatnonzero(~dirich)
    Abd = nbd * DX * DX
    aOpen = math.sqrt(max(1e-6, (gap[openB] * DX).sum()) / math.pi)        # the openings' equivalent radius (above the rim)
    def wedge_family(vertical):
        """Air modes of one vertical family: shapes at the board (mean square 1 there), frequencies,
        the volume integral of psi^2 (per unit pressure) and the radiation loss out of the openings."""
        Lh = LapA.tocsr()
        if vertical:
            Lh = Lh + sp.diags((math.pi / hAv) ** 2 * hAv)
        Lk = Lh[keepN][:, keepN].tocsc()
        Mk = sp.diags(hAv[keepN]).tocsc()
        nmode = 40 if not vertical else 30
        if boxAir:
            Lk = Lk + 1e-6 * sp.eye(len(keepN))
        ev, ew = sla.eigsh(Lk, k=nmode, M=Mk, sigma=0, which="LM")
        psiA = np.zeros((na, nmode)); psiA[keepN] = ew
        if boxAir and not vertical:
            ev = ev + (math.pi / (2 * CASE_H)) ** 2
        psi = psiA[bInA]
        scale = np.sqrt((psi ** 2).sum(0) / nbd)
        psi = psi / scale; psiA = psiA / scale
        fr = C_AIR * np.sqrt(np.maximum(ev, 0)) / (2 * math.pi)
        vfrac = 0.5 if vertical else 1.0
        lam_ = (psiA ** 2 * (hAv * vfrac)[:, None]).sum(0) * DX * DX       # integral of psi^2 over the volume, outside air included
        # Loss out of the released edge: the normal velocity there is grad psi / (i omega rho); power
        # rho c R |v|^2 / 2 over the opening's area, against the stored energy lam_ / (2 rho c^2)
        flux = np.zeros(nmode)
        for k_ in np.flatnonzero(dirich):
            i, j = aij[k_]
            for a, b in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                ii, jj = i + a, j + b
                if 0 <= ii < nxa and 0 <= jj < nya and maskA[ii, jj] and not dirich[posA[ii * nya + jj]]:
                    flux += (psiA[posA[ii * nya + jj]] / DX) ** 2 * (hAv[k_] * vfrac) * DX
                    break
        om = 2 * math.pi * np.maximum(fr, 1.0)
        R_ = np.minimum(1.0, (om / C_AIR * aOpen) ** 2 / 2.0)
        eta_ = C_AIR ** 3 * R_ * flux / (om ** 3 * lam_)
        if boxAir:                                                           # the old open-top box, for A/B
            eta_ = 2 * np.minimum(1.0, (om / C_AIR * math.sqrt(Abd / math.pi)) ** 2 / 2.0) * C_AIR / (om * CASE_H)
        return psi, fr, lam_, eta_
    psi0, fa0, lam0, eta0 = wedge_family(False)
    psi1, fa1, lam1, eta1 = wedge_family(True)
    aw = np.concatenate([psi0, psi1], 1); fa = np.concatenate([fa0, fa1]); LamV = np.concatenate([lam0, lam1])
    etaA_modes = np.concatenate([eta0, eta1])
    nair = aw.shape[1]
    print(f"case air (lid at {math.degrees(LID_ANGLE):.0f} deg, {openB.sum()} open edge nodes, opening radius {aOpen:.2f} m): "
          f"even family from {fa0[:4].round(0)} Hz (loss {eta0[:4].round(2)}), half-wave family from {fa1[:3].round(0)} Hz (loss {eta1[:3].round(2)}); "
          f"loss factors median {np.median(etaA_modes):.3f}")
    full_s = P @ vecs
    Cba = (full_s[bi, :].T @ aw) * DX * DX                                   # (structural modes, air modes)
    ns = len(f)
    Ob = np.diag((2 * math.pi * f) ** 2); Oa = np.diag((2 * math.pi * fa) ** 2)
    Dg = RHO_AIR * C_AIR ** 2 / LamV                                         # per air mode
    # --- The lid: a plate over the case air (see lid_plate). The air's pressure under it
    # pushes it out along its normal, and its motion pumps the air the opposite way to the
    # board's. Its sound reaches the listener from both its faces, with its own short circuit.
    if not os.environ.get("NO_LIDPLATE"):
        fl, Vl, P3l, syl, dAl, dEl = lid_plate(F_TOP)
    else:
        fl, Vl, P3l, syl, dAl, dEl = np.zeros(0), np.zeros((0, 0)), np.zeros((0, 3)), np.zeros((0, 2)), 1.0, np.zeros(0)
    nl = len(fl)
    psiFull = np.zeros((nx * ny, nair)); psiFull[bi] = aw
    onb = inside(P3l[:, 0], P3l[:, 1]) if nl else np.zeros(0, bool)
    Psi_l = interp_rows(np.clip(P3l[:, 0], 0.0, A), np.clip(P3l[:, 1], 0.0, B)) @ psiFull if nl else np.zeros((0, nair))
    if nl:
        Psi_l[~onb] = 0.0
        Psi_l[:, len(fa0):] *= -1.0                                          # the half-wave family: cos(pi) at the lid
    Cla = (Vl.T @ Psi_l) * dAl                                               # (lid modes, air modes)
    Ol = np.diag((2 * math.pi * fl) ** 2)
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
    nsl = ns + nl
    Kc = np.block([[Ob, Z_((ns, nl)), Cba, -Cbu],
                   [Z_((nl, ns)), Ol, -Cla, Z_((nl, nu_))],
                   [Z_((nair, nsl)), Oa, Z_((nair, nu_))],
                   [Z_((nu_, nsl)), Z_((nu_, nair)), Ou]])
    Mc = np.block([[np.eye(nsl), Z_((nsl, nair)), Z_((nsl, nu_))],
                   [-Dg[:, None] * Cba.T, Dg[:, None] * Cla.T, np.eye(nair), Z_((nair, nu_))],
                   [Du * Cbu.T, Z_((nu_, nl)), Z_((nu_, nair)), np.eye(nu_)]])
    import scipy.linalg as la
    lam, Vc = la.eig(Kc, Mc)
    lam, Vc = lam.real, Vc.real
    o = np.argsort(lam); lam, Vc = lam[o], Vc[:, o]
    keep = (lam > (2 * math.pi * 5.0) ** 2) & (lam < lam_top)
    lam, Vc = lam[keep], Vc[:, keep]
    Q = Vc[:ns]                                                              # structural part: the board and its body
    Ql = Vc[ns:nsl]                                                          # ... and the lid
    nq = np.sqrt((Q ** 2).sum(0) + (Ql ** 2).sum(0)); nq[nq == 0] = 1
    # Share of the potential energy held by the air (structure: w^2 q^2 / 2; air: Lam P^2 / (2 rho c^2))
    eS = ((2 * math.pi * f)[:, None] ** 2 * Q ** 2).sum(0) + ((2 * math.pi * fl)[:, None] ** 2 * Ql ** 2).sum(0)
    eA = ((LamV / (RHO_AIR * C_AIR ** 2))[:, None] * Vc[nsl:nsl + nair] ** 2).sum(0) + (LamU / (RHO_AIR * C_AIR ** 2)) * (Vc[nsl + nair:] ** 2).sum(0)
    airE = eA / (eA + eS)
    # Losses of the coupled modes. Structure: each structural mode's material loss, weighted by
    # its share of the stiffness energy. Air: an air mode loses what it radiates out of the case,
    # through the open top (above the board) or the open sides (below it): the opening's
    # radiation resistance, rho c min(1, (ka)^2 / 2) per area.
    w2s = (2 * math.pi * f) ** 2                                             # structural modes (before coupling)
    w2l = (2 * math.pi * fl) ** 2
    etaSc = ((Q ** 2 * (w2s * etaS)[:, None]).sum(0) + (Ql ** 2 * w2l[:, None]).sum(0) * ETA["lid"]) / \
        np.maximum((Q ** 2 * w2s[:, None]).sum(0) + (Ql ** 2 * w2l[:, None]).sum(0), 1e-30)
    frameShare = (Q ** 2 * frE[:, None]).sum(0) / np.maximum((Q ** 2).sum(0) + (Ql ** 2).sum(0), 1e-30)
    lamA = LamV / (RHO_AIR * C_AIR ** 2); lamU = LamU / (RHO_AIR * C_AIR ** 2)
    perim = ring.sum() * DX
    aTop, aSide = math.sqrt(Abd / math.pi), FLOOR_GAP
    def rnorm(fr, a_):
        ka = 2 * math.pi * fr / C_AIR * a_
        return np.minimum(1.0, ka ** 2 / 2.0)
    fc_ = np.sqrt(lam) / (2 * math.pi)
    etaU_modes = rnorm(fu, aSide) * (perim * FLOOR_GAP) * C_AIR / (2 * math.pi * np.maximum(fu, 1.0) * LamU)
    eAir = (lamA[:, None] * (Vc[nsl:nsl + nair] ** 2) * etaA_modes[:, None]).sum(0) + (lamU * (Vc[nsl + nair:] ** 2) * etaU_modes[:, None]).sum(0)
    etaAir = eAir / (eA + eS)                                               # already weighted by the air's share
    etaStruct = etaSc * (1.0 - airE)
    Q = Q / nq; Ql = Ql / nq                                                 # unit structural modal mass (board, body and lid)
    f = fc_
    rimShare = (Q ** 2 * rimE[:, None]).sum(0)
    braceShare = (Q ** 2 * brE[:, None]).sum(0)
    lidShare = (Ql ** 2).sum(0)
    print(f"coupled body: {len(f)} modes; lowest {f[:6].round(1)}; air modes from {fa[0]:.0f} Hz; "
          f"rim-heavy modes {(rimShare > 0.3).sum()}, air-heavy {(airE > 0.5).sum()}, frame-heavy {(frameShare > 0.5).sum()}, lid-heavy {(lidShare > 0.5).sum()}")
    full = full_s @ Q                                       # (nx*ny, n)
    lidV = Vl @ Ql if nl else np.zeros((0, len(f)))         # lid node motion along its normal, per coupled mode
    frontV = vfront @ Q                                     # frame front bar, per coupled mode
    twistQ = vtwist @ Q                                     # rim twist angle at each ring node, per coupled mode

    # Bridge zone points: the middle of each group of keys along each bridge, on the line of the
    # front bridge pins where the strings actually bear (a little in front of the bridge's
    # centre line), so the bridge's rocking is felt as well as its up-and-down motion
    zx, zy = [], []
    bx, by, _ = curve_points(bass_bridge, 2 * N_BASS + 1)
    for z in range(N_BASS):          # zone 0 = lowest keys (A0 end, t = 1)
        k = 2 * (N_BASS - 1 - z) + 1
        zx.append(bx[k]); zy.append(by[k] - PIN_OFFSET)
    tx, ty, _ = curve_points(treble_bridge, 2 * N_TREBLE + 1)
    for z in range(N_TREBLE):        # from F2 (t = 1) up to C8 (t = 0)
        k = 2 * (N_TREBLE - 1 - z) + 1
        zx.append(tx[k]); zy.append(ty[k] - PIN_OFFSET)
    Z = interp_rows(np.array(zx), np.array(zy))
    phi = (Z @ full)                                        # (zones, n)
    # Sideways: the strings' sideways pull acts at the bridge top, one bridge height above the
    # board, and tilts the bridge along its length: the board's slope along the bridge there
    hz = np.array([BASS_W] * N_BASS + [TREBLE_W] * N_TREBLE) + 0.5 * H_MID
    ex = 0.012
    phiH = hz[:, None] * ((interp_rows(np.array(zx) + ex, np.array(zy)) - interp_rows(np.array(zx) - ex, np.array(zy))) @ full) / (2 * ex)
    print(f"sideways bridge shapes: {np.abs(phiH).mean() / np.abs(phi).mean() * 100:.1f}% of the up-and-down ones")
    # The bridge's slope across its width at each zone (its rocking), for reference
    Zb = interp_rows(np.array(zx), np.array(zy) + 2 * PIN_OFFSET)
    rock = np.abs(Z @ full - Zb @ full).mean() / max(np.abs(phi).mean(), 1e-12)
    print(f"bridge rocking: the front and back pin lines differ by {rock * 100:.0f}% of the motion on average")

    # Radiation. The board is open on both faces: sound from its underside reaches a listener
    # above only by bending round the rim (a path longer by about the distance to the edge) or
    # off the floor, and it is in opposite phase. At low frequencies the two nearly cancel (the
    # acoustic short circuit), which the old baffled model missed. The open lid adds a reflection
    # of the top face (an image source behind the lid; it only reflects well above a few hundred Hz).
    xs, ys = X.ravel(), Y.ravel()
    bn = board.ravel()
    bx, by_ = xs[bn], ys[bn]
    dE = dist.ravel()[bn]                                   # distance to the rim [m]
    k = 2 * math.pi * f / C_AIR
    hinge = np.array([LID_HINGE_X, 0.0, RIM_TOP])
    nlid = np.array([math.sin(LID_ANGLE), 0.0, math.cos(LID_ANGLE)])
    tlid = np.array([-math.cos(LID_ANGLE), 0.0, math.sin(LID_ANGLE)])      # up the lid from the hinge
    def lid_image(P3):
        d_ = (P3 - hinge) @ nlid
        return P3 - 2 * d_[:, None] * nlid[None, :]
    nodes3 = np.stack([bx, by_, np.zeros_like(bx)], 1)
    img_lid = lid_image(nodes3)
    # The lid is finite: the reflection from a node reaches a listener only if the mirror
    # point on the lid's plane falls on the lid (hinge to its outer edge, which follows the
    # case's outline, and from the front flap's fold to the tail overhang). Near the lid's
    # edges the reflection fades over a Fresnel zone, sqrt(lambda r1 r2 / (r1 + r2)): long
    # waves spill round the lid, short ones reflect right up to its edge (Kirchhoff).
    def lid_window_setup(pt):
        pt = np.asarray(pt, float)
        dP = (img_lid - hinge) @ nlid; dQ = float((pt - hinge) @ nlid)
        t_ = dP / (dP - dQ + 1e-12)
        spec = img_lid + t_[:, None] * (pt[None, :] - img_lid)
        s_ = (spec - hinge) @ tlid                                           # up the lid [m]
        ys_ = spec[:, 1]
        s_max = (LID_HINGE_X - (np.asarray(bentside(np.clip(ys_, 0, B)), float) - 0.04)) / math.cos(LID_ANGLE)
        d_edge = np.minimum.reduce([s_, s_max - s_, ys_ - LID_FRONT, LID_TAIL - ys_])
        r1 = np.linalg.norm(spec - nodes3, axis=1); r2 = np.linalg.norm(pt[None, :] - spec, axis=1)
        return d_edge, r1 * r2 / np.maximum(r1 + r2, 1e-6)
    from scipy.special import erf
    def lid_window(setup, fi):
        d_edge, rr = setup
        rF = np.sqrt(C_AIR / max(fi, 1.0) * rr)
        return 0.5 * (1.0 + erf(d_edge * math.sqrt(2.0) / rF))
    img_floor = np.stack([bx, by_, np.full_like(bx, -2 * (FLOOR_GAP + 0.13))], 1)
    RIM_H = 0.26
    ds_ring = ds
    ringXY = np.stack([rp[:, 0] * DX, rp[:, 1] * DX, np.full(len(rp), 0.5 * 0.26 - 0.1)], 1)
    def pressure_at(pt, lid=True, floor=True):
        pt = np.array(pt)
        lidSetup = lid_window_setup(pt)
        # The lid panel's own sound: the listener hears the face turned towards them (its
        # motion is counted along the normal, away from the case), with the other face's
        # sound bending round the lid's edge in opposite phase
        sgnL = 1.0 if float((pt - hinge) @ nlid) > 0 else -1.0
        rl1 = np.linalg.norm(P3l - pt, axis=1) if nl else np.zeros(0); rl2 = rl1 + dEl
        rr1 = np.linalg.norm(ringXY - pt, axis=1)
        r1 = np.linalg.norm(nodes3 - pt, axis=1)
        r2 = r1 + dE
        r3 = np.linalg.norm(img_floor - pt, axis=1) + 0.5 * dE
        r4 = np.linalg.norm(img_lid - pt, axis=1)
        out = np.zeros(len(f))
        for i in range(len(f)):
            w = full[bn, i] * area
            kk = k[i]
            p_ = np.sum(w * (np.exp(-1j * kk * r1) / r1 - np.exp(-1j * kk * r2) / r2))
            # The rim's walls: twisting, a wall swings its top in and out (half its height times
            # the angle) and pushes air on both faces; outside and inside are in opposite phase,
            # the inside path longer by about the wall's height
            wr_ = twistQ[:, i] * (0.5 * RIM_H) * (RIM_H * ds_ring)
            if not os.environ.get('NO_RIMRAD'): p_ += np.sum(wr_ * (np.exp(-1j * kk * rr1) / rr1 - np.exp(-1j * kk * (rr1 + RIM_H)) / (rr1 + RIM_H)))
            if floor:
                p_ -= FLOOR_REFL * 0.5 * np.sum(w * np.exp(-1j * kk * r3) / r3)
            if lid:
                p_ += LID_REFL * np.sum(w * lid_window(lidSetup, f[i]) * np.exp(-1j * kk * r4) / r4)
            if nl:
                wl = lidV[:, i] * dAl
                p_ += sgnL * np.sum(wl * (np.exp(-1j * kk * rl1) / rl1 - np.exp(-1j * kk * rl2) / rl2))
            p_ *= kk
            out[i] = abs(p_) * (-1 if p_.real < 0 else 1)
        return out
    radL, radR = pressure_at(LISTEN["L"]), pressure_at(LISTEN["R"])
    micL, micR = pressure_at(MICS["L"], floor=False), pressure_at(MICS["R"], floor=False)
    vol = np.array([np.sum(full[bn, i]) for i in range(len(f))])

    def lit(v):
        t = "%.5g" % v
        if "." not in t and "e" not in t:
            t += ".0"
        return t + "f"

    # Radiation damping: power radiated over the whole sphere, both faces, with the short circuit
    th = (np.arange(8) + 0.5) / 8 * (math.pi / 2)
    ph = np.arange(24) / 24 * 2 * math.pi
    TH, PH = np.meshgrid(th, ph, indexing="ij")
    dOmega = (math.pi / 2 / 8) * (2 * math.pi / 24) * np.sin(TH).ravel()
    ux = (np.sin(TH) * np.cos(PH)).ravel(); uy = (np.sin(TH) * np.sin(PH)).ravel()
    etaRad = np.zeros(len(f))
    for i in range(len(f)):
        kk = k[i]
        w = full[bn, i] * area * (1.0 - np.exp(-1j * kk * dE))
        W = np.exp(-1j * kk * (np.outer(ux, bx) + np.outer(uy, by_))) @ w
        omega = 2 * math.pi * f[i]
        etaRad[i] = 2.0 * omega * RHO_AIR / (2 * math.pi ** 2 * C_AIR) * np.sum(np.abs(W) ** 2 * dOmega) / 2.0
        if nl:                                                               # the lid radiates too, in its own plane
            wl = lidV[:, i] * dAl * (1.0 - np.exp(-1j * kk * dEl))
            Wl = np.exp(-1j * kk * (np.outer(ux, syl[:, 0]) + np.outer(uy, syl[:, 1]))) @ wl
            etaRad[i] += 2.0 * omega * RHO_AIR / (2 * math.pi ** 2 * C_AIR) * np.sum(np.abs(Wl) ** 2 * dOmega) / 2.0
    print("radiation loss factor: median %.4f, 10-90%% %.4f-%.4f" % (np.median(etaRad), *np.percentile(etaRad, [10, 90])))

    # Modes that mostly move the frame go to the frame table; the rest to the board table
    isFrame = frameShare > 0.5
    # Modes that live almost entirely in the lid barely move the bridge, so the strings cannot
    # drive them: leave them out of the table (the lid's answer to the board is in the mixed modes)
    pnorm = np.sqrt((phi ** 2).sum(0) + (phiH ** 2).sum(0))
    inert = (lidShare > 0.5) & (pnorm < 0.1 * np.median(pnorm))
    keepB = ~isFrame & ~inert
    print(f"lid-only modes left out: {int(inert.sum())}; board table {int(keepB.sum())} modes")
    with open(args.out, "w") as fo:
        fo.write("// Generated by tools/soundboard/build_board.py: do not edit by hand.\n")
        fo.write(f"// Body modes below {F_TOP:.0f} Hz of a Yamaha C5: tapered, ribbed soundboard with bridges, crown and\n")
        fo.write("// downbearing, coupled to its rim (bending and twist), back braces, legs, iron frame and the air above\n")
        fo.write("// and below the board, all solved together. Per mode: frequency [Hz]; shape at each bridge zone\n")
        fo.write("// (front pin line; unit modal mass, 1/sqrt(kg)); sound pressure at the player's seat (left,\n")
        fo.write("// right) and at a close mic pair (left, right), with the open board's short circuit, the floor\n")
        fo.write("// and the lid (relative units); net volume; loss factors: radiation, structure (from the\n")
        fo.write("// materials it strains), air (leaving the open case); and the rim and air shares of the mode.\n")
        fo.write(f"static constexpr int kBoardModeCount = {int(keepB.sum())};\n")
        fo.write(f"static constexpr int kBoardZones = {N_BASS + N_TREBLE};\n")
        fo.write("struct BoardModeData { float f; float phi[%d]; float phiH[%d]; float radLeft, radRight, micLeft, micRight, volume, "
                 "radiationLoss, structLoss, airLoss, rimShare, airShare; };\n" % (N_BASS + N_TREBLE, N_BASS + N_TREBLE))
        fo.write("static const BoardModeData kBoardModes[kBoardModeCount] = {\n")
        for i in np.flatnonzero(keepB):
            fo.write("    { %.2ff, { %s }, { %s }, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s },\n" % (
                f[i], ", ".join(lit(v) for v in phi[:, i]), ", ".join(lit(v) for v in phiH[:, i]), lit(radL[i]), lit(radR[i]), lit(micL[i]), lit(micR[i]),
                lit(vol[i]), lit(etaRad[i]), lit(etaStruct[i]), lit(etaAir[i]), lit(rimShare[i] + braceShare[i]), lit(airE[i])))
        fo.write("};\n")
    print("wrote", args.out)
    # Frame table: the frame-heavy body modes, their shapes where each zone's strings end
    axz = np.where(np.arange(N_BASS + N_TREBLE) < N_BASS, 1.25 - 0.08 * np.arange(N_BASS + N_TREBLE), np.array(zx))
    fidx = np.flatnonzero(isFrame & (f < 3000.0))
    fphi = np.array([[np.interp(x_, frontX, frontV[:, i]) for x_ in axz] for i in fidx]).reshape(len(fidx), N_BASS + N_TREBLE)
    frad = np.array([radL[i] for i in fidx])
    write_frame(args.frame_out, f[fidx], fphi, frad, etaStruct[fidx])

    if args.json:
        import json
        yy = np.linspace(0, B, 120)
        outline = GEOM["outline"] if "outline" in GEOM else \
            [[A, 0.0], [0.0, 0.0]] + [[round(float(bentside(v)), 3), round(float(v), 3)] for v in yy] + [[A, B]]
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
