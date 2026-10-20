"""Cut one frame out of a trajectories.xyz of the DEM code and save it as a small text file.

    python code/extract_frame.py runs/B_rate20/trajectories.xyz data/B_frame_20deg.dat --time 1.3
    python code/extract_frame.py runs/B_rate20/trajectories.xyz out.dat --frame 59

The file format is the one of cell 5 of A4.ipynb: per frame a line with the number
of particles, a line with the box and `Time=...`, then one line per particle,
`C x y z radius speed` (SI units). The output has the columns x y z radius speed
(`np.savetxt`, header line starting with `#` gives the frame index, its time and
what was asked for). With --time the frame whose time is nearest is taken.
Frames are read one at a time, so a large trajectory is not held in memory.
"""
import argparse
import re
import sys

import numpy as np


def iter_frames(path):
    """Yield (time, array of shape (N, 5): x y z radius speed) for every frame in the file."""
    with open(path) as f:
        while (line := f.readline()):
            if not line.strip():
                continue
            n = int(line)
            time = float(re.search(r"Time=(\S+)", f.readline()).group(1))
            rows = [f.readline() for _ in range(n)]
            yield time, np.loadtxt(rows, usecols=(1, 2, 3, 4, 5), ndmin=2)


def read_frame(path, time=None, frame=None):
    """Return (index, time, array) of the frame with the given index, or the one nearest to `time`."""
    if (time is None) == (frame is None):
        raise ValueError("give exactly one of time and frame")
    best = None
    for i, (t, data) in enumerate(iter_frames(path)):
        if frame is not None:
            if i == frame:
                return i, t, data
        elif best is None or abs(t - time) < abs(best[1] - time):
            best = (i, t, data)
        elif t > time:      # times increase: the nearest frame has passed
            break
    if best is None:
        raise ValueError(f"frame {frame} not found in {path}")
    return best


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("xyz", help="trajectories.xyz")
    ap.add_argument("out", help="output text file")
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--time", type=float, help="take the frame nearest to this time (s)")
    g.add_argument("--frame", type=int, help="take the frame with this index (0-based)")
    a = ap.parse_args()
    i, t, data = read_frame(a.xyz, a.time, a.frame)
    ask = f"time {a.time} s" if a.time is not None else f"frame {a.frame}"
    np.savetxt(a.out, data, fmt="%.8e",
               header=f"frame {i} at t = {t:.6g} s of {a.xyz} (asked for {ask}); columns: x y z radius speed (SI)")
    print(f"frame {i}, t = {t:.6g} s, {len(data)} particles -> {a.out}", file=sys.stderr)


if __name__ == "__main__":
    main()
