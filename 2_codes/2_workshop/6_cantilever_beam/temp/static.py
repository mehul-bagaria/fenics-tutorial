"""
FEniCS Static Analysis Module

This module performs static structural analysis on imported meshes using FEniCS.
It supports boundary conditions via marker names and exports results in XDMF format.
"""

from dolfin import *
import numpy as np
import json
import os


def static_analysis(
    mesh_file,
    surface_markers_file,
    markers_json_file,
    output_file,
    fixed_marker_name,
    forced_marker_name,
    E=2.1e11,
    nu=0.3,
    f_body=(0.0, 0.0, 0.0),
    f_surface=(0.0, 0.0, -1e6),
    print_results=True
):
    """
    Perform static structural analysis using FEniCS on an imported mesh.
    
    Parameters:
    -----------
    mesh_file : str
        Path to the mesh XDMF file
    surface_markers_file : str
        Path to the surface markers XDMF file
    markers_json_file : str
        Path to the markers JSON file
    output_file : str
        Path to save the results XDMF file
    fixed_marker_name : str
        Name of the fixed boundary marker (e.g., "Fixed")
    forced_marker_name : str
        Name of the forced/loaded boundary marker (e.g., "Forced")
    E : float, optional
        Young's modulus (Pa), default: 2.1e11 (steel)
    nu : float, optional
        Poisson's ratio, default: 0.3
    f_body : tuple of float, optional
        Body force (fx, fy, fz) in N/m³, default: (0, 0, 0)
    f_surface : tuple of float, optional
        Surface traction (fx, fy, fz) in Pa, default: (0, 0, -1e6)
    print_results : bool, optional
        Whether to print analysis results, default: True
    
    Returns:
    --------
    dict
        Dictionary containing:
        - 'displacement': displacement Function
        - 'von_mises_stress': von Mises stress Function
        - 'max_displacement': maximum displacement value
        - 'max_stress': maximum von Mises stress value
        - 'mesh': the mesh object
    """
    
    # Create output directory if it doesn't exist
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    
    if print_results:
        print("=" * 70)
        print("STARTING STATIC ANALYSIS")
        print("=" * 70)
        import sys
        sys.stdout.flush()
    
    # ===== LOAD MARKERS INFORMATION =====
    with open(markers_json_file, 'r') as f:
        markers_info = json.load(f)
    
    # Get marker IDs from names
    if fixed_marker_name not in markers_info['surface']:
        raise ValueError(f"Fixed marker '{fixed_marker_name}' not found in surface markers")
    if forced_marker_name not in markers_info['surface']:
        raise ValueError(f"Forced marker '{forced_marker_name}' not found in surface markers")
    
    FIXED_MARKER = markers_info['surface'][fixed_marker_name]['id']
    FORCED_MARKER = markers_info['surface'][forced_marker_name]['id']
    
    if print_results:
        print(f"\n[1] Markers loaded:")
        print(f"    Fixed marker: '{fixed_marker_name}' = {FIXED_MARKER}")
        print(f"    Forced marker: '{forced_marker_name}' = {FORCED_MARKER}")
        import sys
        sys.stdout.flush()
    
    # ===== LOAD MESH =====
    mesh = Mesh()
    with XDMFFile(mesh_file) as infile:
        infile.read(mesh)
    
    if print_results:
        print(f"\n[2] Mesh loaded:")
        print(f"    Vertices: {mesh.num_vertices()}")
        print(f"    Cells: {mesh.num_cells()}")
        print(f"    Facets: {mesh.num_facets()}")
        import sys
        sys.stdout.flush()
    
    # ===== LOAD SURFACE MARKERS =====
    mvc_surface = MeshValueCollection("int", mesh, 2)
    with XDMFFile(surface_markers_file) as infile:
        infile.read(mvc_surface, "marker")
    
    surface_markers = cpp.mesh.MeshFunctionInt(mesh, mvc_surface)
    
    if print_results:
        print(f"\n[3] Surface markers loaded:")
        print(f"    Unique marker values: {set(surface_markers.array())}")
        import sys
        sys.stdout.flush()
    
    # ===== SETUP FUNCTION SPACE =====
    V = VectorFunctionSpace(mesh, 'P', 1)
    
    if print_results:
        print(f"\n[4] Function space created:")
        print(f"    DOFs: {V.dim()}")
        import sys
        sys.stdout.flush()
    
    # ===== MATERIAL PROPERTIES =====
    # Calculate Lamé parameters
    mu = E / (2.0 * (1.0 + nu))
    lambda_ = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    
    if print_results==True:
        print(f"\n[5] Material properties:")
        print(f"    E = {E:.2e} Pa, nu = {nu}")
        print(f"    μ = {mu:.2e} Pa, λ = {lambda_:.2e} Pa")
        import sys
        sys.stdout.flush()
    
    # ===== BOUNDARY CONDITIONS SETUP =====
    boundaries = MeshFunction("size_t", mesh, mesh.topology().dim()-1, 0)
    
    # Debug: print marker values to verify
    if print_results:
        print(f"\n[6] Setting up boundary conditions:")
        print(f"    Looking for: Fixed={FIXED_MARKER}, Forced={FORCED_MARKER}")
        import sys
        sys.stdout.flush()
    
    fixed_count = 0
    forced_count = 0
    
    # Use the surface_markers directly to set boundaries
    for facet in facets(mesh):
        marker_val = surface_markers[facet]
        if marker_val == FIXED_MARKER:
            boundaries[facet] = 1  # Mark fixed boundary
            fixed_count += 1
        elif marker_val == FORCED_MARKER:
            boundaries[facet] = 2  # Mark forced boundary
            forced_count += 1
    
    if print_results:
        print(f"    Fixed facets found: {fixed_count}")
        print(f"    Forced facets found: {forced_count}")
        import sys
        sys.stdout.flush()
        
    if fixed_count == 0:
        print("WARNING: No fixed boundary facets found! Problem will be singular.")
    if forced_count == 0:
        print("WARNING: No forced boundary facets found! No load will be applied.")
    
    # Apply zero displacement on fixed boundary
    bc = DirichletBC(V, Constant((0, 0, 0)), boundaries, 1)
    
    # ===== LOADING CONDITIONS =====
    # ===== VARIATIONAL FORMULATION =====
    if print_results:
        print(f"\n[7] Solving elasticity problem...")
        print(f"    Body force: {f_body}")
        print(f"    Surface traction: {f_surface}")
        import sys
        sys.stdout.flush()
    
    def epsilon(u):
        """Strain tensor"""
        return 0.5 * (nabla_grad(u) + nabla_grad(u).T)
    
    def sigma(u):
        """Stress tensor"""
        d = u.geometric_dimension()
        return lambda_ * div(u) * Identity(d) + 2 * mu * epsilon(u)
    
    # Define variational problem
    u = TrialFunction(V)
    v = TestFunction(V)
    
    # Bilinear form
    a = inner(sigma(u), epsilon(v)) * dx
    
    # Linear form (loads)
    f_body_const = Constant(f_body)
    f_surface_const = Constant(f_surface)
    
    ds = Measure('ds', domain=mesh, subdomain_data=boundaries)
    L = dot(f_body_const, v) * dx + dot(f_surface_const, v) * ds(2)
    
    # ===== SOLVE THE SYSTEM =====
    u_sol = Function(V, name="Displacement")
    solve(a == L, u_sol, bc)
    
    if print_results:
        print(f"    ✓ Solution computed!")
        import sys
        sys.stdout.flush()
    
    # ===== POST-PROCESSING: STRESS CALCULATIONS =====
    if print_results:
        print(f"\n[8] Post-processing...")
        import sys
        sys.stdout.flush()
    
    def von_mises_stress(u):
        """Calculate von Mises stress from displacement"""
        s = sigma(u) - (1.0/3.0) * tr(sigma(u)) * Identity(3)  # Deviatoric stress
        von_mises = sqrt(3.0/2.0 * inner(s, s))
        return von_mises
    
    V_scalar = FunctionSpace(mesh, 'P', 1)
    von_mises_expr = von_mises_stress(u_sol)
    von_mises_sol = project(von_mises_expr, V_scalar)
    von_mises_sol.rename("VonMisesStress", "von_mises_stress")
    
    # Also compute stress components
    sigma_expr = sigma(u_sol)
    
    # ===== CALCULATE RESULT STATISTICS =====
    u_magnitude = sqrt(dot(u_sol, u_sol))
    u_magnitude_proj = project(u_magnitude, V_scalar)
    
    max_displacement = u_sol.vector().max()
    max_stress = von_mises_sol.vector().max()
    
    # ===== EXPORT RESULTS =====
    if print_results:
        print(f"\n[9] Saving results to: {output_file}")
        import sys
        sys.stdout.flush()
    
    with XDMFFile(output_file) as xdmf:
        xdmf.parameters["flush_output"] = True
        xdmf.parameters["functions_share_mesh"] = True
        xdmf.write(u_sol, 0.0)
        xdmf.write(von_mises_sol, 0.0)
        xdmf.write(u_magnitude_proj, 0.0)
    
    # ===== RESULTS SUMMARY =====
    if print_results:
        print(f"\n[10] RESULTS:")
        print(f"     Max displacement: {max_displacement:.6e} m")
        print(f"     Max stress: {max_stress:.6e} Pa ({max_stress/1e6:.2f} MPa)")
        print("=" * 70)
        print("ANALYSIS COMPLETE!")
        print("=" * 70)
        import sys
        sys.stdout.flush()
    
    # Return results
    return {
        'displacement': u_sol,
        'von_mises_stress': von_mises_sol,
        'max_displacement': max_displacement,
        'max_stress': max_stress,
        'mesh': mesh,
        'markers_info': markers_info
    }
