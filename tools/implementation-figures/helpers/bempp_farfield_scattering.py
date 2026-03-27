import argparse
import numpy as np
import bempp_cl.api as bempp
from bempp_cl.api.operators.boundary import helmholtz as boundary
from bempp_cl.api.operators.far_field import helmholtz as far_field
from bempp_cl.api.linalg import gmres as bempp_gmres
from scipy.sparse.linalg import gmres as scipy_gmres


def parse_args():
    parser = argparse.ArgumentParser(
        description = "BEMPP far-field scattering for single-target pressure-release or fluid cases."
    )
    parser.add_argument("--mesh", required = True, help = "Input Gmsh 2.2 surface mesh.")
    parser.add_argument(
        "--boundary",
        required = True,
        choices = ["pressure_release", "liquid_filled"],
        help = "Boundary type handled by this BEMPP runner."
    )
    parser.add_argument("--incidence-deg", type = float, default = 90.0)
    parser.add_argument("--frequency", type = float, required = True)
    parser.add_argument("--sound-speed-sw", type = float, required = True)
    parser.add_argument("--density-sw", type = float, default = 1026.8)
    parser.add_argument("--sound-speed-body", type = float)
    parser.add_argument("--density-body", type = float)
    parser.add_argument("--n-angles", type = int, default = 181)
    parser.add_argument(
        "--plane",
        choices = ["xy"],
        default = "xy",
        help = "Receive-angle plane. Only xy is currently implemented."
    )
    return parser.parse_args()


def incident_direction(theta_deg):
    theta = np.pi / 180.0 * theta_deg
    return np.array([np.cos(theta), np.sin(theta), 0.0])


def xy_plane_directions(n_angles):
    phi = np.linspace(0.0, 2.0 * np.pi, num = n_angles, endpoint = True)
    directions = np.array([[np.cos(v), np.sin(v), 0.0] for v in phi]).T
    return phi, directions


def solve_pressure_release(space, k0, d):
    @bempp.complex_callable
    def pinc(x, n, idx, res):
        res[0] = np.exp(1j * k0 * d.dot(x))

    f = bempp.GridFunction(space, fun = pinc)
    S = boundary.single_layer(space, space, space, k0)
    u, info = bempp_gmres(S, -f)
    if info != 0:
        raise RuntimeError(f"Pressure-release GMRES did not converge cleanly (info={info}).")

    return {"surface_density": u}


def solve_liquid_filled(space, k0, k1, rho0, rho1, d):
    @bempp.complex_callable
    def pinc(x, n, idx, res):
        res[0] = np.exp(1j * k0 * d.dot(x))

    @bempp.complex_callable
    def dpinc(x, n, idx, res):
        res[0] = 1j * k0 * np.exp(1j * k0 * d.dot(x)) * d.dot(n)

    f = bempp.GridFunction(space, fun = pinc)
    g = bempp.GridFunction(space, fun = dpinc)

    A0 = boundary.multitrace_operator(space.grid, k0)
    A1 = boundary.multitrace_operator(space.grid, k1)
    A1[0, 1] *= rho1 / rho0
    A1[1, 0] *= rho0 / rho1

    A = (A0 + A1).strong_form()
    rhs = np.concatenate([f.coefficients, g.coefficients])
    u, info = scipy_gmres(A.dot(A), A.dot(rhs))
    if info != 0:
        raise RuntimeError(f"Fluid GMRES did not converge cleanly (info={info}).")

    return {
        "psi": bempp.GridFunction(space, coefficients = u[: space.global_dof_count]),
        "phi": bempp.GridFunction(space, coefficients = u[space.global_dof_count :])
    }


def main():
    args = parse_args()
    d = incident_direction(args.incidence_deg)
    k0 = 2.0 * np.pi * args.frequency / args.sound_speed_sw

    grid = bempp.import_grid(args.mesh)
    space = bempp.function_space(grid, "P", 1)

    phi, recv_dirs = xy_plane_directions(args.n_angles)

    if args.boundary == "pressure_release":
        solution = solve_pressure_release(space, k0, d)
        FF = far_field.single_layer(space, recv_dirs, k0)
        scattered = -FF.evaluate(solution["surface_density"])
    else:
        if args.sound_speed_body is None or args.density_body is None:
            raise ValueError("Fluid-filled runs require --sound-speed-body and --density-body.")

        k1 = 2.0 * np.pi * args.frequency / args.sound_speed_body
        solution = solve_liquid_filled(
            space = space,
            k0 = k0,
            k1 = k1,
            rho0 = args.density_sw,
            rho1 = args.density_body,
            d = d
        )

        SL = far_field.single_layer(space, recv_dirs, k0)
        DL = far_field.double_layer(space, recv_dirs, k0)
        scattered = DL.evaluate(solution["psi"]) - SL.evaluate(solution["phi"])

    print("angle_deg\tphi_scatter_rad\tmod_f")
    for i in range(scattered.shape[1]):
        print(f"{phi[i] * 180.0 / np.pi}\t{phi[i]}\t{abs(scattered[0][i])}")


if __name__ == "__main__":
    main()
