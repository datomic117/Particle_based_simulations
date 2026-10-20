"""Render a sliced snapshot of a DEM trajectory (trajectories.xyz) with OVITO, coloured by speed.

    python code/render_ovito.py runs/D7/trajectories.xyz data/C_snapshot_settled.png --time 0.286 --vmax 0.002 --zrange 50 115
    python code/render_ovito.py runs/D7/trajectories.xyz out.png --frame 13

Needs the `ovito` Python module (pip install ovito). The file states the box, the radii and the speed itself, so no
column mapping is needed. What is drawn (the recipe of the assignment text):
  * Slice, normal (0, 1, 0), slab width `--slab` (default 4 mm), through the centre of the box: only the particles whose
    centre lies in that slab are drawn, so the inside of the silo is visible.
  * Color coding by `speed` (m/s) from 0 to `--vmax` (default: the largest speed in the slab of this frame), with a colour legend.
  * Orthographic view along y (x to the right, z up), Tachyon renderer, the box drawn as thin lines.
Frame choice: `--frame N` (0-based) or `--time T` (the frame nearest to T seconds); default is the last frame.
Crop: `--zrange ZMIN ZMAX` (mm) shows only that height range, `--xhalf W` (mm) only the band of half-width W around
the axis (default: the whole box).
Used by Part C (settled bed) and Part D (discharge snapshots).
"""
import argparse
import sys

import numpy as np
from ovito.io import import_file
from ovito.modifiers import ColorCodingModifier, SliceModifier
from ovito.qt_compat import QtCore
from ovito.vis import ColorLegendOverlay, TachyonRenderer, TextLabelOverlay, Viewport

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument("trajectory", help="trajectories.xyz of a run")
ap.add_argument("output", help="output png")
g = ap.add_mutually_exclusive_group()
g.add_argument("--frame", type=int, help="frame index (0-based)")
g.add_argument("--time", type=float, help="take the frame nearest to this time (s)")
ap.add_argument("--slab", type=float, default=4.0, help="slab width in mm (default 4)")
ap.add_argument("--vmax", type=float, help="upper end of the speed colour scale in m/s (default: largest speed in the slab)")
ap.add_argument("--zrange", type=float, nargs=2, metavar=("ZMIN", "ZMAX"), help="height range to show, mm")
ap.add_argument("--xhalf", type=float, help="half-width of the shown band around the axis, mm")
ap.add_argument("--size", type=int, nargs=2, default=(1000, 1000), metavar=("W", "H"), help="image size in pixels")
a = ap.parse_args()

pipeline = import_file(a.trajectory)
n_frames = pipeline.source.num_frames
times = [pipeline.source.compute(f).attributes["Time"] for f in range(n_frames)]
if a.frame is not None:
    frame = a.frame
elif a.time is not None:
    frame = int(np.argmin(np.abs(np.array(times) - a.time)))
else:
    frame = n_frames - 1
if not 0 <= frame < n_frames:
    sys.exit(f"frame {frame} not in 0..{n_frames - 1}")

box = pipeline.source.compute(frame).cell
Lx, Ly, Lz = box[0, 0], box[1, 1], box[2, 2]
pipeline.add_to_scene()
slab_m = a.slab * 1e-3
pipeline.modifiers.append(SliceModifier(normal=(0, 1, 0), distance=0.5 * Ly, slab_width=slab_m))

# speed range: given, or the largest speed of the particles that are drawn
if a.vmax is None:
    sliced = pipeline.compute(frame)
    vmax = float(np.max(sliced.particles["speed"]))
else:
    vmax = a.vmax
colour = ColorCodingModifier(property="speed", start_value=0.0, end_value=vmax, gradient=ColorCodingModifier.Viridis())
pipeline.modifiers.append(colour)

pipeline.source.data.cell_.vis.enabled = False  # the box lines add nothing to a cropped view

# view along +y: x to the right, z up; ortho `fov` is half the visible extent of the smaller image dimension
W, H = a.size
vp = Viewport(type=Viewport.Type.Ortho, camera_dir=(0, 1, 0))
z0, z1 = (0.0, Lz) if a.zrange is None else (a.zrange[0] * 1e-3, a.zrange[1] * 1e-3)
half_h = 0.5 * (z1 - z0)
half_w = 0.5 * Lx if a.xhalf is None else a.xhalf * 1e-3
# keep the whole width Lx (or the crop height) visible, whichever needs the larger fov
vp.camera_pos = (0.5 * Lx, 0.5 * Ly, 0.5 * (z0 + z1))
vp.fov = max(half_h, half_w * H / W) * 1.02
legend = ColorLegendOverlay(modifier=colour, title="speed (m/s)",
                            alignment=QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignTop,
                            orientation=QtCore.Qt.Orientation.Vertical, font_size=0.07, format_string="%.3g")
vp.overlays.append(legend)
vp.overlays.append(TextLabelOverlay(text="speed (m/s)", alignment=QtCore.Qt.AlignmentFlag.AlignRight | QtCore.Qt.AlignmentFlag.AlignTop,
                                    offset_y=-0.325, font_size=0.016, text_color=(0, 0, 0)))  # legend title (the legend's own title is not drawn in this layout)
vp.render_image(size=(W, H), filename=a.output, frame=frame, renderer=TachyonRenderer(), background=(1, 1, 1))
print(f"wrote {a.output} (frame {frame}, t = {times[frame]:.4f} s, slab {a.slab} mm at y = {0.5 * Ly * 1e3:.2f} mm, speed scale 0..{vmax:.4g} m/s)")
