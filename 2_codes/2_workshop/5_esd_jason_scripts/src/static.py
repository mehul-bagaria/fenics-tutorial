# %%
from dolfin import *
import json, os

def run_plane_stress_simulation(config_path):
    """
    Runs a 2D plane stress problem in FEniCS using parameters from a JSON config.
    Supports configurable element type and polynomial degree.
    Writes displacement & stress fields to XDMF and a summary JSON.
    """

    # --- Load configuration ---
    with open(config_path) as f:
        cfg = json.load(f)

    mesh_file = cfg["mesh_files"]["mesh"]
    facet_file = cfg["mesh_files"]["facet_markers"]
    point_file = cfg["mesh_files"]["point_markers"]
    out_dir = cfg["output"]["results_dir"]
    out_json = cfg["output"]["output_json"]

    os.makedirs(out_dir, exist_ok=True)

    # --- Load mesh ---
    mesh = Mesh()
    with XDMFFile(mesh_file) as f:
        f.read(mesh)

    # --- Load markers ---
    mvc_facet = MeshValueCollection("size_t", mesh, 1)
    with XDMFFile(facet_file) as f:
        f.read(mvc_facet, "name_to_read")
    facet_markers = cpp.mesh.MeshFunctionSizet(mesh, mvc_facet)

    mvc_point = MeshValueCollection("size_t", mesh, 0)
    with XDMFFile(point_file) as f:
        f.read(mvc_point, "name_to_read")
    point_markers = cpp.mesh.MeshFunctionSizet(mesh, mvc_point)

    # --- Material parameters ---
    E = cfg["material"]["E"]
    nu = cfg["material"]["nu"]
    t = cfg["material"].get("thickness", 1.0)

    mu = E / (2.0 * (1.0 + nu))
    lmbda = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu))
    lmbda_eff = (2 * mu * lmbda) / (lmbda + 2 * mu)  # plane stress correction

    # --- Element family and degree ---
    elem_family = cfg["elements"].get("family", "CG")
    elem_degree = cfg["elements"].get("degree", 1)

    # --- Function space ---
    V = VectorFunctionSpace(mesh, elem_family, elem_degree)
    u = Function(V)
    du = TrialFunction(V)
    v = TestFunction(V)

    # --- Constitutive relations ---
    def eps(u):
        return sym(grad(u))

    def sigma(u):
        return lmbda_eff * tr(eps(u)) * Identity(2) + 2 * mu * eps(u)

    # --- Measures ---
    dx_ = Measure("dx", domain=mesh)
    ds_ = Measure("ds", domain=mesh, subdomain_data=facet_markers)
    dP_ = Measure("ds", domain=mesh, subdomain_data=point_markers)

    # --- Variational form ---
    a = inner(sigma(du), eps(v)) * dx_
    L = Constant(0.0) * v[0] * dx_

    # --- Boundary conditions and loads ---
    bcs = []
    for entity in cfg["entities"]:
        dim = entity["dimension"]
        marker = entity["marker"]

        if entity["type"] == "Dirichlet":
            val = Constant(tuple(entity["value"]))
            if dim == 1:
                bc = DirichletBC(V, val, facet_markers, marker)
            elif dim == 0:
                bc = DirichletBC(V, val, point_markers, marker)
            bcs.append(bc)

        elif entity["type"] == "Neumann":
            if dim == 1:
                t_vec = Constant(tuple(entity["traction"]))
                L += dot(t_vec, v) * ds_(marker)
            elif dim == 0 and "force" in entity:
                f_vec = Constant(tuple(entity["force"]))
                L += dot(f_vec, v) * dP_(marker)

    # --- Solve system ---
    solve(a == L, u, bcs)

    # --- Compute stress field ---
    Vsig = TensorFunctionSpace(mesh, "DG", 0)
    sig = project(sigma(u), Vsig)

    # --- Save results ---
    with XDMFFile(os.path.join(out_dir, "displacement.xdmf")) as f:
        f.write(u)
    with XDMFFile(os.path.join(out_dir, "stress.xdmf")) as f:
        f.write(sig)

    # --- Post-processing: magnitudes and extrema ---
    u_mag = sqrt(dot(u, u))
    u_mag_func = project(u_mag, FunctionSpace(mesh, "DG", 0))
    max_disp = float(u_mag_func.vector().get_local().max())

    sxx, syy, sxy = sig[0,0], sig[1,1], sig[0,1]
    von_mises = sqrt(sxx**2 - sxx*syy + syy**2 + 3*sxy**2)
    vm_func = project(von_mises, FunctionSpace(mesh, "DG", 0))
    max_vm = float(vm_func.vector().get_local().max())

    # --- Write summary JSON ---
    summary = {
        "max_displacement": max_disp,
        "max_von_mises_stress": max_vm,
        "mesh_info": {
            "num_cells": mesh.num_cells(),
            "num_vertices": mesh.num_vertices(),
            "element_family": elem_family,
            "element_degree": elem_degree
        },
        "output_files": {
            "displacement": os.path.join(out_dir, "displacement.xdmf"),
            "stress": os.path.join(out_dir, "stress.xdmf")
        }
    }

    with open(out_json, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"✅ Simulation complete. Results written to {out_dir}")
    return summary