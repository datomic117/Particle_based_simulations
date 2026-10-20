"""Part D reductions of a silo trajectory (trajectories.xyz) to small files in data/.

    python code/analyse_d.py runs/D7/trajectories.xyz D7 --z-silo 0.0530997
    python code/analyse_d.py runs/D3/trajectories.xyz D3 --z-silo 0.0530997 --stats

Particle i is line i of every frame (the code never reorders the particles, and the radius column is checked to be identical
in all frames), so a particle can be followed through the file. Writes, with the tag given as second argument:
  data/D_exit_<tag>.csv      one row per particle: x0 y0 z0 (position in the frame nearest t_open = 0.3 s), radius,
                             t_exit = time of the first frame with z < z_silo (nan if it never leaves), t_exit is therefore
                             known to one frame spacing only
  data/D_speedmap_<tag>.csv  mean speed and mean downward velocity of the grains in bins of (r, z - z_silo), r = distance from
                             the silo axis, averaged over the frames of the time window --t0 to --t1 (default 0.58 to 1.92 s,
                             the fit window of D2); bins with fewer than 5 grain-frames are left out
  data/D_speedstats_<tag>.csv  per frame, only with --stats: time, number of grains above the plate, and the median, 99th percentile and
                             maximum of the speed of those grains (m/s) and the position of the fastest one; used for the blocked orifice of D4 (a grain below the plate
                             is falling and is not counted)
"""
import argparse
import re

import numpy as np


def frames(path):
    with open(path) as f:
        while (line := f.readline()):
            if not line.strip():
                continue
            n = int(line)
            hdr = f.readline()
            time = float(re.search(r"Time=(\S+)", hdr).group(1))
            lat = [float(x) for x in re.search(r'Lattice="([^"]+)"', hdr).group(1).split()]
            yield time, (lat[0], lat[4], lat[8]), np.loadtxt([f.readline() for _ in range(n)], usecols=(1, 2, 3, 4, 5), ndmin=2)


ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("xyz")
ap.add_argument("tag")
ap.add_argument("--z-silo", type=float, required=True)
ap.add_argument("--t-open", type=float, default=0.3)
ap.add_argument("--t0", type=float, default=0.58)
ap.add_argument("--t1", type=float, default=1.92)
ap.add_argument("--bin", type=float, default=1e-3, help="bin size in m (default 1 mm)")
ap.add_argument("--stats", action="store_true", help="write the per-frame speed statistics (D_speedstats_<tag>.csv) and nothing else")
a = ap.parse_args()

times, pos0, rad0, t_exit = [], None, None, None
stats = []
sum_v, sum_vz, cnt = {}, {}, {}
prev = None
best = 1e9
for k, (t, box, d) in enumerate(frames(a.xyz)):
    times.append(t)
    if rad0 is None:
        rad0 = d[:, 3].copy()
        t_exit = np.full(len(d), np.nan)
    assert np.array_equal(d[:, 3], rad0), "particle order changed"
    if abs(t - a.t_open) < best:
        best = abs(t - a.t_open)
        pos0 = d[:, :3].copy()
    inside = d[:, 2] > a.z_silo
    jf = np.flatnonzero(inside)[np.argmax(d[inside, 4])]          # the fastest grain above the plate and where it is
    stats.append((t, inside.sum(), np.median(d[inside, 4]), np.percentile(d[inside, 4], 99), d[jf, 4],
                  np.hypot(d[jf, 0] - 0.5 * box[0], d[jf, 1] - 0.5 * box[1]), d[jf, 2] - a.z_silo))
    gone = (d[:, 2] < a.z_silo) & np.isnan(t_exit)
    t_exit[gone] = t
    if a.t0 <= t <= a.t1 and prev is not None:
        # velocity from the displacement over one frame (centred on the mid-point), speed from the file
        vz = (d[:, 2] - prev[1][:, 2]) / (t - prev[0])
        mid = 0.5 * (d[:, :3] + prev[1][:, :3])
        r = np.hypot(mid[:, 0] - 0.5 * box[0], mid[:, 1] - 0.5 * box[1])
        ib = np.floor(r / a.bin).astype(int)
        jb = np.floor((mid[:, 2] - a.z_silo) / a.bin).astype(int)
        for i_, j_, s_, w_ in zip(ib, jb, d[:, 4], vz):
            key = (i_, j_)
            sum_v[key] = sum_v.get(key, 0.0) + s_
            sum_vz[key] = sum_vz.get(key, 0.0) + w_
            cnt[key] = cnt.get(key, 0) + 1
    prev = (t, d)
if a.stats:
    np.savetxt(f"data/D_speedstats_{a.tag}.csv", stats, delimiter=",", fmt=["%.6f", "%d", "%.6e", "%.6e", "%.6e", "%.6e", "%.6e"],
               header="t,n_silo,v_median,v_p99,v_max,r_fast,zrel_fast  (grains with z > z_silo; speed in m/s; position of the fastest grain: distance from the axis and height above the plate, m)", comments="# ")
    print(f"{len(stats)} frames written to data/D_speedstats_{a.tag}.csv")
    raise SystemExit
np.savetxt(f"data/D_exit_{a.tag}.csv", np.column_stack([pos0, rad0, t_exit]), delimiter=",", fmt="%.8e",
           header=f"x0,y0,z0,radius,t_exit  (x0,y0,z0 in the frame nearest t_open = {a.t_open} s; t_exit = first frame with z < z_silo = {a.z_silo} m, nan if none; SI units)", comments="# ")
keys = sorted(k for k in cnt if cnt[k] >= 5)
rows = [((i + 0.5) * a.bin, (j + 0.5) * a.bin, sum_v[(i, j)] / cnt[(i, j)], sum_vz[(i, j)] / cnt[(i, j)], cnt[(i, j)]) for i, j in keys]
np.savetxt(f"data/D_speedmap_{a.tag}.csv", rows, delimiter=",", fmt="%.6e",
           header=f"r,z_rel,speed,vz,count  (bin {a.bin} m in r and in z - z_silo; mean over the frames with {a.t0} <= t <= {a.t1} s; SI units; bins with < 5 grain-frames omitted)", comments="# ")
print(f"{len(times)} frames, {np.isfinite(t_exit).sum()} of {len(t_exit)} grains passed z_silo, {len(rows)} map bins")
