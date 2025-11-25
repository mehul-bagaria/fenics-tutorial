# FEniCS Static Analysis

## How to Run

### Prerequisites
- Docker must be running in the background

### Workflow

1. **Create mesh using Salome:**
   - Open the `.hdf` file from the `mesh/` folder in Salome
   - Adjust the element size in meshing settings
   - Export the mesh as `.med` format
   - Save with a descriptive name in the `mesh/` folder (e.g., `mesh_coarse.med`, `mesh_fine.med`)

2. **Configure simulation parameters:**
   - Open `simulation_config/study_config.json`
   - Update the following paths:
     - `"med_mesh_file"`: Path to your mesh file (e.g., `"mesh/mesh_coarse.med"`)
     - `"mesh_file"`: Output folder for converted mesh (e.g., `"mesh/output-1/mesh.xdmf"`)
     - `"output_displacement_file"`: Results folder (e.g., `"post-processing/output-1/displacement.xdmf"`)
     - `"output_vonmises_file"`: Results folder (e.g., `"post-processing/output-1/vonmises.xdmf"`)

3. **Run the simulation:**
   ```
   python main.py
   ```

### Mesh Convergence Study

To compare different mesh densities:

1. Create multiple mesh files in Salome with different element sizes:
   - `mesh/mesh_coarse.med`
   - `mesh/mesh_medium.med`
   - `mesh/mesh_fine.med`

2. For each mesh, update the config file with unique output folders:
   ```json
   {
     "med_mesh_file": "mesh/mesh_coarse.med",
     "mesh_file": "mesh/output-1/mesh.xdmf",
     "output_displacement_file": "post-processing/output-1/displacement.xdmf",
     "output_vonmises_file": "post-processing/output-1/vonmises.xdmf"
   }
   ```

3. Run `python main.py` for each configuration

This approach keeps results from different mesh sizes organized in separate folders (`output-1`, `output-2`, etc.) for easy comparison.