# ESD Engine – Mini FEniCS Solver

## overview

This repository implements a lightweight, modular FEniCS-based solver framework designed to run structural or multiphysics studies inside a Docker environment. The solver reads study configurations from JSON files, loads domain meshes and marker data, applies boundary conditions, and exports results in XDMF format.



## folder structure

```
esd_engine/
├─ main.py                  # Main launcher and Docker entrypoint
├─ src/                     # Core solver modules
│  └─ solver.py             # FEniCS problem setup and solve
├─ study_config/
│  ├─ input.json            # Input configuration (physics, BCs, solver options)
│  └─ output.json           # Output configuration (fields, probes, integrals)
├─ mesh/
│  ├─ domain.xdmf           # Domain mesh
│  ├─ markers_facet.xdmf    # Facet markers
│  └─ markers.json          # Physical group mapping
└─ outputs/
   └─ run_001/              # Simulation results and summary
```



## workflow

1. **Prepare study configuration**

   * Define physics parameters, solver settings, and boundary conditions in `study_config/input.json`.
   * Define desired output fields and probes in `study_config/output.json`.

2. **Prepare mesh**

   * Place domain and marker files (`.xdmf`, `.h5`, and `markers.json`) inside the `mesh/` directory.
   * The `markers.json` file should map physical names (e.g., “fixed”, “loaded_edge”) to their numeric IDs.

3. **Run simulation**

   ```bash
   python3 main.py --input study_config/input.json --output study_config/output.json
   ```

   The launcher will automatically start a FEniCS Docker container and execute the solver.

4. **View results**

   * Displacement and stress fields are stored in `outputs/run_xxx/` as `.xdmf` files.
   * A `summary.json` file contains extracted probe values and integrated quantities (e.g., strain energy).



## docker requirements

* Docker must be installed and accessible from the command line.
* The project uses the **FEniCS 2019.1.0 stable** image:

  ```
  quay.io/fenicsproject/stable:2019.1.0
  ```


Command to run the docker image - `docker run -it -v D:\1_codes\0_github\fenics-tutorial\2_codes\2_workshop\5_esd_jason_scripts\:/root/ -w /root/ iitrabhi/fenics`