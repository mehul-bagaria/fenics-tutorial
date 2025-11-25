import meshio
import json
import os


def convert_med_to_xdmf(mesh_input, mesh_output_dir, markers_output_file=None, print_markers=True):
    """
    Convert MED mesh file to XDMF format and extract markers.
    
    Parameters:
    -----------
    mesh_input : str
        Path to the input MED mesh file
    mesh_output_dir : str
        Directory where output XDMF files will be saved
    markers_output_file : str, optional
        Path to save the markers JSON file. If None, will be saved as 'markers.json' in mesh_output_dir
    print_markers : bool, optional
        Whether to print the markers to console (default: True)
    
    Returns:
    --------
    dict
        Dictionary containing surface and solid markers
    """
    
    # Set default markers output file if not provided
    if markers_output_file is None:
        markers_output_file = os.path.join(mesh_output_dir, "markers.json")
    
    # Create output directory if it doesn't exist
    os.makedirs(mesh_output_dir, exist_ok=True)
    
    print(f"Reading mesh from: {mesh_input}")
    
    # Read the MED mesh file
    msh = meshio.read(mesh_input)
    
    # Extract mesh components
    points = msh.points
    cell = msh.cells
    cell_data = msh.cell_data
    
    print(f"Mesh loaded: {len(points)} points, cell types: {list(dict(cell).keys())}")
    
    # --- Write mesh files ---
    
    # Write plain volume mesh (no tags)
    if "tetra" in dict(cell):
        tetra_cells = dict(cell)["tetra"]
        mesh_file = os.path.join(mesh_output_dir, "mesh.xdmf")
        meshio.write(mesh_file, meshio.Mesh(
            points=points[:, :],
            cells={"tetra": tetra_cells},
        ))
        print(f"  ✓ Created: {mesh_file}")
    
    # Write surface mesh with cell data
    if "triangle" in dict(cell):
        triangle_cells = dict(cell)["triangle"]
        triangle_tags = cell_data["triangle"]["cell_tags"]
        surface_file = os.path.join(mesh_output_dir, "surface.xdmf")
        meshio.write(surface_file, meshio.Mesh(
            points=points[:, :],
            cells={"triangle": triangle_cells},
            cell_data={"triangle": {"marker": triangle_tags}}
        ))
        print(f"  ✓ Created: {surface_file}")
    
    # Write solid mesh with cell data
    if "tetra" in dict(cell):
        tetra_cells = dict(cell)["tetra"]
        tetra_tags = cell_data["tetra"]["cell_tags"]
        solid_file = os.path.join(mesh_output_dir, "solid.xdmf")
        meshio.write(solid_file, meshio.Mesh(
            points=points[:, :],
            cells={"tetra": tetra_cells},
            cell_data={"tetra": {"marker": tetra_tags}}
        ))
        print(f"  ✓ Created: {solid_file}")
    
    print(f"\nConversion complete. Mesh files saved in: {mesh_output_dir}")
    
    # --- Extract markers ---
    
    print("\nExtracting markers...")
    
    mesh = msh
    points, cell, cell_data, field_data = mesh.points, mesh.cells, mesh.cell_data, mesh.field_data
    cell_tags = mesh.cell_tags  # group → ids mapping
    
    # Mapping of cell types to geometric dimensions
    cell_type_dim = {
        "triangle": 2,
        "tetra": 3
    }
    
    # Build marker mapping (only for surface=2D and solid=3D)
    marker_map = {}
    for tag, names in cell_tags.items():
        tag_int = int(tag)
        found_dim = -1
        for cell_type, dim in cell_type_dim.items():
            tag_array = cell_data.get(cell_type, {}).get("cell_tags", None)
            if tag_array is not None and tag_int in tag_array:
                found_dim = dim
                break
        if found_dim in (2, 3):  # keep only surface and solid
            for name in names:
                marker_map[name] = {
                    "id": tag_int,
                    "dim": found_dim,
                    "name": name
                }
    
    # Separate surface vs solid
    surface_markers = {n: info for n, info in marker_map.items() if info["dim"] == 2}
    solid_markers   = {n: info for n, info in marker_map.items() if info["dim"] == 3}
    
    # Combine into one JSON
    all_markers = {"surface": surface_markers, "solid": solid_markers}
    
    # Save to JSON file
    with open(markers_output_file, "w") as f:
        json.dump(all_markers, f, indent=4)
    
    print(f"  ✓ Markers saved to: {markers_output_file}")
    
    # Print markers if requested
    if print_markers:
        print("\n" + "=" * 60)
        print("MARKERS:")
        print("=" * 60)
        print(json.dumps(all_markers, indent=4))
        print("=" * 60)
        print(f"\nSummary:")
        print(f"  - Surface markers: {len(surface_markers)}")
        for name, info in surface_markers.items():
            print(f"    • {name}: ID={info['id']}, Dim={info['dim']}")
        print(f"  - Solid markers: {len(solid_markers)}")
        for name, info in solid_markers.items():
            print(f"    • {name}: ID={info['id']}, Dim={info['dim']}")
    
    return all_markers