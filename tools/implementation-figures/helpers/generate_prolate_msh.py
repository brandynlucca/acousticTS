from pathlib import Path
import math


def write_prolate_msh(path, semimajor_length, semiminor_length, n_u = 40, n_v = 48):
    nodes = []

    # North pole on +x semimajor tip
    nodes.append((semimajor_length, 0.0, 0.0))

    # Interior latitude rings
    for i in range(1, n_u):
        u = math.pi * i / n_u
        x = semimajor_length * math.cos(u)
        r = semiminor_length * math.sin(u)
        for j in range(n_v):
            v = 2 * math.pi * j / n_v
            y = r * math.cos(v)
            z = r * math.sin(v)
            nodes.append((x, y, z))

    # South pole on -x semimajor tip
    south_pole = len(nodes) + 1
    nodes.append((-semimajor_length, 0.0, 0.0))

    def ring_start(i):
        # i = 1..n_u-1
        return 2 + (i - 1) * n_v

    def ring_idx(i, j):
        return ring_start(i) + (j % n_v)

    triangles = []

    # Top cap
    north_pole = 1
    for j in range(n_v):
        triangles.append((north_pole, ring_idx(1, j + 1), ring_idx(1, j)))

    # Middle bands
    for i in range(1, n_u - 1):
        for j in range(n_v):
            a = ring_idx(i, j)
            b = ring_idx(i + 1, j)
            c = ring_idx(i + 1, j + 1)
            d = ring_idx(i, j + 1)
            triangles.append((a, b, c))
            triangles.append((a, c, d))

    # Bottom cap
    last_ring = n_u - 1
    for j in range(n_v):
        triangles.append((south_pole, ring_idx(last_ring, j), ring_idx(last_ring, j + 1)))

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
            f.write(f"{idx} 2 2 0 0 {a} {b} {c}\n")
        f.write("$EndElements\n")

    return len(nodes), len(triangles)


if __name__ == "__main__":
    out = Path(r"c:\Users\Brandyn\Desktop\acousticTS_upd\scratch\prolate_axis_x_pr_a35_b10_40x48.msh")
    n_nodes, n_triangles = write_prolate_msh(
        out,
        semimajor_length = 0.035,
        semiminor_length = 0.01,
        n_u = 40,
        n_v = 48
    )
    print(out)
    print(n_nodes, n_triangles)
