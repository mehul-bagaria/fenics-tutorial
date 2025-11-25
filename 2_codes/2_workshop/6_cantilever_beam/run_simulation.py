"""
Static Analysis of Cantilever Beam

This script performs static structural analysis using FEniCS on the imported mesh.
It reads configuration from a JSON file and executes the analysis.
"""

import sys
import json
import os

# Add scripts directory to path
sys.path.append('scripts')

# Import the static analysis function with material properties
from static_3d import static_analysis_with_materials


def main():
    """Main function to run the static analysis."""
    
    # Load analysis configuration from JSON file
    config_file = "simulation_config/study_config.json"
    
    print("=" * 70)
    print("CANTILEVER BEAM STATIC ANALYSIS")
    print("=" * 70)
    print(f"\nLoading configuration from: {config_file}")
    
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    # Get mesh file path and derive other file paths from its location
    mesh_file = config["mesh"]["mesh_file"]
    mesh_dir = os.path.dirname(mesh_file)
    
    # Parse boundaries - find the fixed boundary (Displacement Constraints)
    fixed_boundary = None
    for boundary in config["boundaries"]:
        if boundary["definition"]["type"] == "Displacement Constraints":
            fixed_boundary = boundary
            break
    
    if fixed_boundary is None:
        raise ValueError("No fixed boundary (Displacement Constraints) found in configuration")
    
    fixed_marker_name = fixed_boundary["assigned"][0]  # e.g., "Fixed end"
    
    # Parse loads - extract the force load definition
    force_load = None
    for load in config["loads"]:
        if "FX" in load["definition"]:  # Identify force load
            force_load = load
            break
    
    if force_load is None:
        raise ValueError("No force load found in configuration")
    
    forced_marker_name = force_load["assigned"][0]  # e.g., "Forced end"
    f_surface = (
        float(force_load["definition"]["FX"]),
        float(force_load["definition"]["FY"]),
        float(force_load["definition"]["FZ"])
    )
    
    # Parse materials - extract material properties
    material = None
    for mat in config["materials"]:
        if mat["properties"]["type"] == "Linear Elastic":
            material = mat
            break
    
    if material is None:
        raise ValueError("No Linear Elastic material found in configuration")
    
    material_marker_name = material["assigned"][0]  # e.g., "Domain"
    E = float(material["properties"]["E"])
    nu = float(material["properties"]["nu"])
    
    # Flatten into function parameters
    analysis_config = {
        # File paths
        "mesh_file": mesh_file,
        "surface_markers_file": os.path.join(mesh_dir, "surface.xdmf"),
        "solid_markers_file": os.path.join(mesh_dir, "solid.xdmf"),
        "markers_json_file": os.path.join(mesh_dir, "markers.json"),
        "output_file": "post-processing/output/cantilever_results.xdmf",
        
        # Boundary condition marker names
        "fixed_marker_name": fixed_marker_name,
        "forced_marker_name": forced_marker_name,
        "material_marker_name": material_marker_name,
        
        # Material properties
        "E": E,
        "nu": nu,
        
        # Loading conditions
        "f_body": (0.0, 0.0, 0.0),  # No body forces
        "f_surface": f_surface,
        
        # Output control
        "print_results": True
    }
    
    print(f"\nConfiguration loaded successfully:")
    print(f"  Fixed boundary: '{fixed_marker_name}'")
    print(f"  Forced boundary: '{forced_marker_name}'")
    print(f"  Material domain: '{material_marker_name}'")
    print(f"  Material: E={E}, nu={nu}")
    print(f"  Surface force: {f_surface}")
    
    print("\n" + "=" * 70)
    print("RUNNING STATIC ANALYSIS")
    print("=" * 70 + "\n")
    
    # Run static analysis with material properties using config dictionary
    results = static_analysis_with_materials(**analysis_config)
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nResults saved to: {analysis_config['output_file']}")
    
    if results:
        print("\nSummary:")
        for key, value in results.items():
            print(f"  {key}: {value}")
    
    return results


if __name__ == "__main__":
    try:
        results = main()
        sys.exit(0)
    except Exception as e:
        print(f"\n{'='*70}")
        print("ERROR OCCURRED")
        print(f"{'='*70}")
        print(f"Error: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)
