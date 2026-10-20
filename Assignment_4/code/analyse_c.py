"""Part C follow-up reductions of a silo trajectory (trajectories.xyz) to small files in data/.

    python code/analyse_c.py settling runs/D7/trajectories.xyz data/D7.dat data/C_settling_frames.csv
    python code/analyse_c.py crossing runs/D7/trajectories.xyz data/D7.dat data/C_crossing_edge.csv
    python code/analyse_c.py seeds data/C_edge_seeds.csv

The table (second argument) is only read for its first line, which records R_cyl, z_silo, t_open and the box. The format of the
trajectory is the one of extract_frame.py: per frame a line with N, a line with the box and `Time=...`, then `C x y z radius speed`.

settling  For every frame from --t0 (default 0.15 s) to t_open, the packing of the specified core of task C4 (centres closer than
          R_cyl - 2 d_p to the axis, between z_silo + 2 d_p and 3 d_p below the highest centre of that frame; d_p = 2.0 mm):
          phi = summed particle volume of the centres in the core / core volume, rho_b = phi * 2500 kg/m^3, the coordination
          number of the core particles (contacts |r_i - r_j| < R_i + R_j, all neighbours counted, no wall contacts) and of all
          particles, the highest centre z_top, the 99th percentile of z and the mean height z_cm of all centres (a robust measure
          of the bed height), and the number of centres in the core. Columns of the csv: frame time z_top z_p99 z_cm phi_core
          rho_b Z_core Z_all N_core.
crossing  For every grain that leaves the silo (first frame with z < z_silo): the time of that frame and the distance r_c of the grain
          from the axis at the moment its centre crossed the plane of the plate, found by linear interpolation of x, y and z between
          the last frame above and the first frame below the plate (frames are 0.022 s apart, a grain falls and drifts several mm in
          that time: r_c is an estimate good to about a radius at best). Columns: particle t_frame r_cross.
seeds     The silo runs with the edge of the hole on and off for five random seeds (12345: data/D7.dat and data/C_noedge_D7.dat;
          1 to 4: data/C_seed<s>_edge.dat and _noedge.dat, made by run_c_seeds.sh), reduced with the definition of Q of Part D2
          (the functions slope_rate, fit_window and rate are copied from cell 21 of A4.ipynb): straight-line fit of N_silo(t) from
          10 % to 60 % of the particles gone (the last row if the run ends first), error = standard error of the mean of the slopes
          of 4 equal pieces of the window. Columns: seed edge(1 = edge on) Q sQ t10 (time of 10 % gone) gone_end (fraction gone at the
          end of the window) N_end Z_open Ekin_open (the last table line before t_open).
"""
import argparse
import re
import sys

import numpy as np
from scipy.spatial import cKDTree

D_P, RHO_P = 2.0e-3, 2500.0      # mean particle diameter (task text) and material density


def frames(path):
    with open(path) as f:
        while (line := f.readline()):
            if not line.strip():
                continue
            n = int(line)
            hdr = f.readline()
            time = float(re.search(r"Time=(\S+)", hdr).group(1))
            yield time, np.loadtxt([f.readline() for _ in range(n)], usecols=(1, 2, 3, 4, 5), ndmin=2)


def geometry(table):
    hdr = open(table).readline()
    par = lambda key: float(re.search(key + r"\s+(\S+)", hdr).group(1))
    return dict(R_cyl=par("R_cyl"), z_silo=par("z_silo"), t_open=par("t_open"), axis=0.5 * par("box"))


def contacts_per_particle(xyz, R):
    pairs = cKDTree(xyz).query_pairs(2.0 * R.max(), output_type="ndarray")
    touch = np.linalg.norm(xyz[pairs[:, 0]] - xyz[pairs[:, 1]], axis=1) < R[pairs[:, 0]] + R[pairs[:, 1]]
    return np.bincount(pairs[touch].ravel(), minlength=len(R))


def settling(a, geo):
    rows = []
    for k, (t, d) in enumerate(frames(a.xyz)):
        if t < a.t0 or t > geo["t_open"]:
            continue
        x, y, z, R, _ = d.T
        rc = np.hypot(x - geo["axis"], y - geo["axis"])
        Vp = 4.0 / 3.0 * np.pi * R**3
        r_max, z_lo, z_hi = geo["R_cyl"] - 2 * D_P, geo["z_silo"] + 2 * D_P, z.max() - 3 * D_P
        sel = (rc < r_max) & (z > z_lo) & (z < z_hi)
        phi = Vp[sel].sum() / (np.pi * r_max**2 * (z_hi - z_lo))
        nc = contacts_per_particle(np.c_[x, y, z], R)
        rows.append([k, t, z.max(), np.percentile(z, 99), z.mean(), phi, phi * RHO_P, nc[sel].mean(), nc.mean(), sel.sum()])
    with open(a.out, "w", newline="") as fh:        # newline="": LF also on Windows
        np.savetxt(fh, rows, delimiter=",", fmt=["%d"] + ["%.8g"] * 9,
                   header="frame,time,z_top,z_p99,z_cm,phi_core,rho_b,Z_core,Z_all,N_core")
    print(f"{len(rows)} frames -> {a.out}", file=sys.stderr)


def crossing(a, geo):
    prev, crossed, out = None, None, []
    for k, (t, d) in enumerate(frames(a.xyz)):
        if prev is None:
            crossed = np.zeros(len(d), bool)
        new = (d[:, 2] < geo["z_silo"]) & ~crossed
        if prev is not None and new.any():
            z0, z1 = prev[new, 2], d[new, 2]
            f = (z0 - geo["z_silo"]) / (z0 - z1)             # fraction of the step at which the centre reaches the plate plane
            x = prev[new, 0] + f * (d[new, 0] - prev[new, 0])
            y = prev[new, 1] + f * (d[new, 1] - prev[new, 1])
            r = np.hypot(x - geo["axis"], y - geo["axis"])
            out += [[i, t, ri] for i, ri in zip(np.flatnonzero(new), r)]
        crossed |= new
        prev = d
    with open(a.out, "w", newline="") as fh:
        np.savetxt(fh, out, delimiter=",", fmt=["%d", "%.6g", "%.6g"], header="particle,t_frame,r_cross")
    print(f"{len(out)} grains crossed -> {a.out}", file=sys.stderr)


N0 = 7500


def slope_rate(t, N):
    """Particles per second leaving: minus the slope of a straight-line fit of N(t)."""
    return -np.polyfit(t, N, 1)[0]


def fit_window(tab, lo=0.10, hi=0.60):
    """First row with 10 % of the particles gone to the first row with 60 % gone (the last row if the run ends before)."""
    gone = N0 - tab[:, 6]
    i0 = int(np.argmax(gone >= lo * N0))
    i1 = int(np.argmax(gone >= hi * N0)) if np.any(gone >= hi * N0) else len(tab) - 1
    return i0, i1


def rate(tab, K=4):
    """Q from a fit of N(t) over the window; error = standard error of the mean of the slopes of K equal pieces."""
    i0, i1 = fit_window(tab)
    t, N = tab[:, 1], tab[:, 6]
    Q = slope_rate(t[i0:i1 + 1], N[i0:i1 + 1])
    edges = np.linspace(t[i0], t[i1], K + 1)
    parts = np.array([slope_rate(t[(t >= lo) & (t <= hi)], N[(t >= lo) & (t <= hi)]) for lo, hi in zip(edges[:-1], edges[1:])])
    return Q, parts.std(ddof=1) / np.sqrt(K), t[i0], (N0 - N[i1]) / N0


def seeds(a):
    files = {(12345, 1): "data/D7.dat", (12345, 0): "data/C_noedge_D7.dat"}
    for sd in (1, 2, 3, 4):
        files[(sd, 1)], files[(sd, 0)] = f"data/C_seed{sd}_edge.dat", f"data/C_seed{sd}_noedge.dat"
    rows = []
    for (sd, edge), path in files.items():
        tab = np.loadtxt(path)
        t_open = float(re.search(r"t_open\s+(\S+)", open(path).readline()).group(1))
        Q, sQ, t10, gone_end = rate(tab)
        k = np.flatnonzero(tab[:, 1] < t_open)[-1]
        rows.append([sd, edge, Q, sQ, tab[np.argmax(N0 - tab[:, 6] >= 0.10 * N0), 1], gone_end, tab[-1, 6], tab[k, 2], tab[k, 4]])
    with open(a.out, "w", newline="") as fh:
        np.savetxt(fh, rows, delimiter=",", fmt=["%d", "%d"] + ["%.8g"] * 7, header="seed,edge,Q,sQ,t10,gone_end,N_end,Z_open,Ekin_open")
    print(f"{len(rows)} runs -> {a.out}", file=sys.stderr)


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
sub = ap.add_subparsers(dest="what", required=True)
for name in ("settling", "crossing"):
    p = sub.add_parser(name)
    p.add_argument("xyz")
    p.add_argument("table")
    p.add_argument("out")
    p.add_argument("--t0", type=float, default=0.15, help="settling: first frame time (s)")
sub.add_parser("seeds").add_argument("out")
a = ap.parse_args()
if a.what == "seeds":
    seeds(a)
else:
    {"settling": settling, "crossing": crossing}[a.what](a, geometry(a.table))
