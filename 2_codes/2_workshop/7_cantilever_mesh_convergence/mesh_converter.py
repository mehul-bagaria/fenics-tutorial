import sys
import json
import os
sys.path.append('scripts')

# Import the mesh converter function
from meshio_converter_med2xdmf import convert_med_to_xdmf

# Load simulation config
with open('simulation_config/study_config.json', 'r') as f:
    config = json.load(f)

# Get mesh file path from config
med_mesh_file = config['mesh']['med_mesh_file']
mesh_output_file = config['mesh']['mesh_file']

# Extract output directory from the mesh output file path
mesh_output_dir = os.path.dirname(mesh_output_file)

# Run conversion - only specify input and output directory
markers = convert_med_to_xdmf(
    mesh_input=med_mesh_file,
    mesh_output_dir=mesh_output_dir
)