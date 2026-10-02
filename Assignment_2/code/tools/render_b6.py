from ovito.io import import_file
from ovito.vis import Viewport, TachyonRenderer

# Load the B6 starting configuration
pipeline = import_file("../b6_start.pdb")

# Add the pipeline to the OVITO scene
pipeline.add_to_scene()

# Create a perspective camera
vp = Viewport(type=Viewport.Type.Perspective)

# Choose a nice viewing direction
vp.camera_dir = (-1.0, -1.0, -0.8)

# Automatically frame the whole system
vp.zoom_all()

# Render to a PNG file
vp.render_image(
    filename="../b6_render.png",
    size=(1400, 1000),
    background=(1, 1, 1),
    renderer=TachyonRenderer()
)

print("Rendered image saved as code/b6_render.png")