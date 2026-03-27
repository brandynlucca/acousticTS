from pathlib import Path
import math


def write_cylinder_msh(path,
                       length_body,
                       radius_body,
                       n_circ = 24,
                       n_ax = 12,
                       n_rad = 1):
    nodes = []

    for i in range(n_ax + 1):
        x = -length_body / 2 + length_body * i / n_ax
        for j in range(n_circ):
            ang = 2 * math.pi * j / n_circ
            y = radius_body * math.cos(ang)
            z = radius_body * math.sin(ang)
            nodes.append((x, y, z))

    def add_cap_rings(x_cap):
        cap_rings = []
        # Keep the outer boundary on the sidewall ring, then add interior
        # concentric rings so the flat endcaps are not represented by a single
        # fan of triangles.
        for ir in range(n_rad - 1, 0, -1):
            radius_i = radius_body * ir / n_rad
            start_idx = len(nodes) + 1
            for j in range(n_circ):
                ang = 2 * math.pi * j / n_circ
                y = radius_i * math.cos(ang)
                z = radius_i * math.sin(ang)
                nodes.append((x_cap, y, z))
            cap_rings.append(start_idx)
        center_idx = len(nodes) + 1
        nodes.append((x_cap, 0.0, 0.0))
        return cap_rings, center_idx

    left_cap_rings, left_center = add_cap_rings(-length_body / 2)
    right_cap_rings, right_center = add_cap_rings(length_body / 2)

    def ring_idx(i, j):
        return i * n_circ + (j % n_circ) + 1

    def cap_ring_idx(cap_rings, ring_number, j):
        return cap_rings[ring_number] + (j % n_circ)

    triangles = []

    for i in range(n_ax):
        for j in range(n_circ):
            a = ring_idx(i, j)
            b = ring_idx(i + 1, j)
            c = ring_idx(i + 1, j + 1)
            d = ring_idx(i, j + 1)
            triangles.append((a, b, c))
            triangles.append((a, c, d))

    def add_cap_triangles(outer_start, cap_rings, center_idx, reverse = False):
        rings = [outer_start] + cap_rings
        for ring_number in range(len(rings) - 1):
            outer_ring_start = rings[ring_number]
            inner_ring_start = rings[ring_number + 1]
            for j in range(n_circ):
                if ring_number == 0:
                    a = outer_ring_start(j)
                    d = outer_ring_start(j + 1)
                else:
                    a = cap_ring_idx(cap_rings, ring_number - 1, j)
                    d = cap_ring_idx(cap_rings, ring_number - 1, j + 1)
                b = cap_ring_idx(cap_rings, ring_number, j)
                c = cap_ring_idx(cap_rings, ring_number, j + 1)

                if reverse:
                    triangles.append((a, c, b))
                    triangles.append((a, d, c))
                else:
                    triangles.append((a, b, c))
                    triangles.append((a, c, d))

        if cap_rings:
            inner_most_ring = len(cap_rings) - 1
            for j in range(n_circ):
                a = cap_ring_idx(cap_rings, inner_most_ring, j)
                b = cap_ring_idx(cap_rings, inner_most_ring, j + 1)
                if reverse:
                    triangles.append((center_idx, a, b))
                else:
                    triangles.append((center_idx, b, a))
        else:
            for j in range(n_circ):
                a = outer_start(j)
                b = outer_start(j + 1)
                if reverse:
                    triangles.append((center_idx, a, b))
                else:
                    triangles.append((center_idx, b, a))

    add_cap_triangles(
        outer_start = lambda j: ring_idx(0, j),
        cap_rings = left_cap_rings,
        center_idx = left_center,
        reverse = False
    )
    add_cap_triangles(
        outer_start = lambda j: ring_idx(n_ax, j),
        cap_rings = right_cap_rings,
        center_idx = right_center,
        reverse = True
    )

    path = Path(path)
    with path.open("w", encoding = "ascii") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write("$Nodes\n")
        f.write(f"{len(nodes)}\n")
        for idx, (x, y, z) in enumerate(nodes, start = 1):
            f.write(f"{idx} {x:.12g} {y:.12g} {z:.12g}\n")
        f.write("$EndNodes\n")
        f.write("$Elements\n")
        f.write(f"{len(triangles)}\n")
        for idx, (a, b, c) in enumerate(triangles, start = 1):
            # Gmsh 2.2 triangle record with two dummy tags so mm-bem parsing matches.
            f.write(f"{idx} 2 2 0 0 {a} {b} {c}\n")
        f.write("$EndElements\n")

    return len(nodes), len(triangles)


if __name__ == "__main__":
    out = Path(r"c:\Users\Brandyn\Desktop\acousticTS_upd\scratch\cylinder_axis_x_pr_70x10_24x12.msh")
    n_nodes, n_triangles = write_cylinder_msh(
        out,
        length_body = 0.07,
        radius_body = 0.01,
        n_circ = 24,
        n_ax = 12,
        n_rad = 6
    )
    print(out)
    print(n_nodes, n_triangles)
