"""Render a snapshot of a DPD trajectory (PDB) with OVITO, beads coloured by type.

    python code/render_ovito.py data/C2_a37_traj.pdb data/C2_ovito_snapshot.png

Needs the `ovito` Python module (pip install ovito). The last frame of the file is rendered, looking
along z at a slab of thickness 3 so that both interfaces (perpendicular to x) are visible.
Type A beads are red and type B beads blue (the PDB writer gives A the element H and B the element X).
"""
import sys

from ovito.io import import_file
from ovito.modifiers import SliceModifier
from ovito.vis import TachyonRenderer, Viewport

src, out = sys.argv[1], sys.argv[2]
slab_width = float(sys.argv[3]) if len(sys.argv) > 3 else 3.0

pipeline = import_file(src)
pipeline.add_to_scene()
pipeline.modifiers.append(SliceModifier(normal=(0, 0, 1), distance=slab_width / 2, slab_width=slab_width))

frame = pipeline.source.num_frames - 1
# colours and radii are set on the pipeline source, so every frame uses them
# (element H = type A, element X = type B, see record_trajectories_pdb)
src_types = pipeline.source.data.particles_.particle_types_
for ptype in src_types.types_:
    ptype.radius = 0.35
    ptype.color = (0.85, 0.15, 0.15) if ptype.name.upper().startswith("H") else (0.15, 0.35, 0.85)

pipeline.source.data.cell_.vis.line_width = 0.08
vp = Viewport(type=Viewport.Type.Ortho, camera_dir=(0, 0, -1))
vp.zoom_all((1600, 800))
vp.render_image(size=(1600, 800), filename=out, frame=frame, renderer=TachyonRenderer(), background=(1, 1, 1))
print(f"wrote {out} (frame {frame})")
