import numpy as np


def read_vtk_unstructured(filename):
    """
    Read an ASCII legacy VTK file with DATASET UNSTRUCTURED_GRID,
    POINTS, CELLS, CELL_TYPES and CELL_DATA (SCALARS p, rho, Mach, VECTORS velocity).

    Returns
    -------
    points : (N, 3) float array
        Node coordinates.
    cells : (Nc, 4) int array
        Quad cell connectivity (node indices).
    scalars : dict(str -> (Nc,) float array)
        Cell-centered scalar fields (e.g. 'p', 'rho', 'Mach').
    vectors : dict(str -> (Nc, 3) float array)
        Cell-centered vector fields (e.g. 'velocity').
    """
    with open(filename, "r") as f:
        # Header (4 lines)
        _ = f.readline()  # "# vtk DataFile Version 3.0"
        _ = f.readline()  # title
        _ = f.readline()  # "ASCII"
        line = f.readline()  # "DATASET UNSTRUCTURED_GRID"
        if not line.strip().startswith("DATASET UNSTRUCTURED_GRID"):
            raise ValueError("VTK file is not an UNSTRUCTURED_GRID dataset.")

        # Skip blank lines until POINTS
        line = f.readline()
        while line and not line.strip():
            line = f.readline()

        parts = line.split()
        if parts[0] != "POINTS":
            raise ValueError("Expected 'POINTS' section in VTK file.")

        n_points = int(parts[1])
        # parts[2] is data type (e.g. "double"), ignore

        # Read all point coordinates (3 * n_points numbers)
        coords = []
        while len(coords) < 3 * n_points:
            line = f.readline()
            if not line:
                raise EOFError("Unexpected end of file while reading POINTS.")
            tokens = line.split()
            if not tokens:
                continue
            coords.extend(map(float, tokens))

        points = np.array(coords, dtype=float).reshape((n_points, 3))

        # Next: CELLS
        line = f.readline()
        while line and not line.strip():
            line = f.readline()
        parts = line.split()
        if parts[0] != "CELLS":
            raise ValueError("Expected 'CELLS' section in VTK file.")

        n_cells = int(parts[1])
        n_entries = int(parts[2])  # = n_cells * (1 + nverts) for quads

        cell_list = []
        read_cells = 0
        entries_read = 0
        while read_cells < n_cells:
            line = f.readline()
            if not line:
                raise EOFError("Unexpected end of file while reading CELLS.")
            tokens = line.split()
            if not tokens:
                continue
            idx = 0
            while idx < len(tokens):
                nverts = int(tokens[idx])
                idx += 1
                if nverts != 4:
                    raise ValueError("This post-process assumes VTK_QUAD cells (4 vertices).")
                node_ids = list(map(int, tokens[idx:idx + nverts]))
                idx += nverts
                cell_list.append(node_ids)
                read_cells += 1
                entries_read += 1 + nverts
        cells = np.array(cell_list, dtype=int)

        # CELL_TYPES
        line = f.readline()
        while line and not line.strip():
            line = f.readline()
        parts = line.split()
        if parts[0] != "CELL_TYPES":
            raise ValueError("Expected 'CELL_TYPES' section.")
        n_types = int(parts[1])

        read_types = 0
        while read_types < n_types:
            line = f.readline()
            if not line:
                raise EOFError("Unexpected end of file while reading CELL_TYPES.")
            tokens = line.split()
            if not tokens:
                continue
            read_types += len(tokens)
        # We assume all types are 9 (VTK_QUAD) and ignore them.

        # CELL_DATA
        line = f.readline()
        while line and not line.strip():
            line = f.readline()
        parts = line.split()
        if parts[0] != "CELL_DATA":
            raise ValueError("Expected 'CELL_DATA' section.")
        n_cell_data = int(parts[1])

        scalars = {}
        vectors = {}

        # Now parse all remaining fields
        for line in f:
            if not line.strip():
                continue
            parts = line.split()
            key = parts[0]
            if key == "SCALARS":
                # Format: SCALARS name type num_components
                name = parts[1]
                # type_name = parts[2]
                # ncomp = int(parts[3])
                # Next line: LOOKUP_TABLE default
                lut_line = next(f)

                values = []
                while len(values) < n_cell_data:
                    l = next(f)
                    toks = l.split()
                    if not toks:
                        continue
                    values.extend(map(float, toks))
                scalars[name] = np.array(values, dtype=float)

            elif key == "VECTORS":
                # Format: VECTORS name type
                name = parts[1]
                # type_name = parts[2]
                comps = []
                while len(comps) < 3 * n_cell_data:
                    l = next(f)
                    toks = l.split()
                    if not toks:
                        continue
                    comps.extend(map(float, toks))
                vectors[name] = np.array(comps, dtype=float).reshape((n_cell_data, 3))
            else:
                # Ignore unknown sections
                continue

    return points, cells, scalars, vectors


def _build_boundary_loops(cells):
    """
    From quad connectivity, build boundary edges and closed loops.

    Parameters
    ----------
    cells : (Nc, 4) int array

    Returns
    -------
    segments : list of dict
        Each dict has keys: 'start', 'end', 'cell', 'edge'.
    loops_nodes : list of list[int]
        Node indices (with last possibly equal to first) for each loop.
    loops_segments : list of list[int]
        Indices into `segments` for each loop, in the same order as loops_nodes.
    """
    edges_map = {}  # (min(i,j), max(i,j)) -> list of (cell_idx, local_edge, start, end)
    for c_idx, nodes in enumerate(cells):
        for e in range(4):
            ni = nodes[e]
            nj = nodes[(e + 1) % 4]
            key = (min(ni, nj), max(ni, nj))
            edges_map.setdefault(key, []).append((c_idx, e, ni, nj))

    segments = []
    for key, usage in edges_map.items():
        # Boundary edges appear only once
        if len(usage) == 1:
            c_idx, e, ni, nj = usage[0]
            segments.append(
                {
                    "start": ni,
                    "end": nj,
                    "cell": c_idx,
                    "edge": e,
                }
            )

    # Map from start node to all segment indices starting there
    start_map = {}
    for k, seg in enumerate(segments):
        start_map.setdefault(seg["start"], []).append(k)

    visited = [False] * len(segments)
    loops_nodes = []
    loops_segments = []

    for k in range(len(segments)):
        if visited[k]:
            continue

        seg = segments[k]
        visited[k] = True

        loop_nodes = [seg["start"], seg["end"]]
        loop_segs = [k]
        current_end = seg["end"]

        while True:
            candidates = [idx for idx in start_map.get(current_end, []) if not visited[idx]]
            if not candidates:
                break
            next_idx = candidates[0]
            seg_next = segments[next_idx]
            visited[next_idx] = True
            loop_segs.append(next_idx)
            loop_nodes.append(seg_next["end"])
            current_end = seg_next["end"]
            # Close loop
            if current_end == loop_nodes[0]:
                break

        if len(loop_segs) > 0:
            loops_nodes.append(loop_nodes)
            loops_segments.append(loop_segs)

    return segments, loops_nodes, loops_segments


def _pick_inner_airfoil_loop(points, loops_nodes):
    """
    Pick the smallest-area closed loop -> interpreted as the airfoil boundary.

    Parameters
    ----------
    points : (N, 3)
    loops_nodes : list of list[int]

    Returns
    -------
    loop_index : int
        Index of chosen loop in loops_nodes.
    """
    if not loops_nodes:
        raise RuntimeError("No boundary loops found in mesh.")

    best_idx = None
    best_area = None

    for i, node_list in enumerate(loops_nodes):
        if len(node_list) < 4:
            continue
        idxs = np.array(node_list, dtype=int)
        # Remove closing duplicate, if present
        if idxs[0] == idxs[-1]:
            idxs = idxs[:-1]
        coords = points[idxs, :2]
        x = coords[:, 0]
        y = coords[:, 1]
        area = 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)
        if best_idx is None or abs(area) < best_area:
            best_idx = i
            best_area = abs(area)

    if best_idx is None:
        raise RuntimeError("Could not determine inner airfoil boundary loop.")

    return best_idx


def compute_coeff_from_vtk(
    vtk_filename,
    Mach_inf,
    alpha,
    T_inf,
    p_inf,
    chord=1.0,
    gamma=1.4,
    R=287.0,
):
    """
    Compute Cp distribution and C_L, C_D, C_M from an unstructured VTK file
    produced by writeSolutionVTKPhysical (quad cells, cell-centered p, rho, Mach, velocity).

    Parameters
    ----------
    vtk_filename : str
        Path to the VTK solution file (ASCII UNSTRUCTURED_GRID).
    Mach_inf : float
        Freestream Mach number.
    alpha : float
        Angle of attack in radians.
    T_inf : float
        Freestream static temperature [K].
    p_inf : float
        Freestream static pressure [Pa].
    chord : float, optional
        Reference chord length [m], default 1.0.
    gamma : float, optional
        Ratio of specific heats, default 1.4.
    R : float, optional
        Gas constant [J/(kg·K)], default 287.

    Returns
    -------
    x_surf : (Ns,) float array
        x-coordinates of airfoil surface nodes (loop order).
    y_surf : (Ns,) float array
        y-coordinates of airfoil surface nodes (loop order).
    Cp_surf : (Ns,) float array
        Cp distribution at surface nodes.
    C_L : float
        Lift coefficient.
    C_D : float
        Drag coefficient.
    C_M : float
        Pitching moment coefficient about x = chord/4, y = 0 (z out of plane).
    """
    points, cells, scalars, vectors = read_vtk_unstructured(vtk_filename)

    if "p" not in scalars:
        raise KeyError("VTK file has no scalar field named 'p'.")
    if "rho" not in scalars:
        # not strictly needed for Cp (we use freestream rho), so we don't enforce it
        rho_cells = None
    else:
        rho_cells = scalars["rho"]

    p_cells = scalars["p"]
    n_nodes = points.shape[0]
    n_cells = cells.shape[0]

    # Build boundary loops
    segments, loops_nodes, loops_segments = _build_boundary_loops(cells)
    loop_idx = _pick_inner_airfoil_loop(points, loops_nodes)

    node_loop = np.array(loops_nodes[loop_idx], dtype=int)
    seg_loop = loops_segments[loop_idx]

    # Remove closing duplicate if present
    if node_loop[0] == node_loop[-1]:
        node_loop = node_loop[:-1]

    # Map each node to the list of cells that touch it
    node_to_cells = [[] for _ in range(n_nodes)]
    for c_idx, conn in enumerate(cells):
        for n in conn:
            node_to_cells[n].append(c_idx)

    # Average cell-centered pressure to nodal pressure (only need boundary nodes)
    p_nodes = {}
    for n in node_loop:
        attached = node_to_cells[n]
        if not attached:
            raise RuntimeError(f"No cells attached to boundary node {n}.")
        p_nodes[n] = float(np.mean(p_cells[attached]))

    # Freestream quantities
    a_inf = np.sqrt(gamma * R * T_inf)
    U_inf = Mach_inf * a_inf
    rho_inf = p_inf / (R * T_inf)
    q_inf = 0.5 * rho_inf * U_inf ** 2

    # Coordinates and Cp on surface nodes (loop order)
    x_surf = points[node_loop, 0]
    y_surf = points[node_loop, 1]
    p_surf = np.array([p_nodes[n] for n in node_loop])
    Cp_surf = (p_surf - p_inf) / q_inf

    # Force and moment integration over the boundary loop
    # Reference point for moment: quarter-chord
    x_ref = chord / 4.0
    y_ref = 0.0

    # Approximate body centroid from boundary nodes (for outward normals)
    centroid = np.array([np.mean(x_surf), np.mean(y_surf)])

    Fx = 0.0
    Fy = 0.0
    Mz = 0.0

    n_edges = len(node_loop)
    for k in range(n_edges):
        i0 = node_loop[k]
        i1 = node_loop[(k + 1) % n_edges]

        x0, y0 = points[i0, 0], points[i0, 1]
        x1, y1 = points[i1, 0], points[i1, 1]

        # Local edge geometry
        dx = x1 - x0
        dy = y1 - y0
        ds = np.hypot(dx, dy)
        if ds == 0.0:
            continue

        # Node pressures and mid-point pressure
        p0 = p_nodes[i0]
        p1 = p_nodes[i1]
        p_mid = 0.5 * (p0 + p1)

        # Outward normal from body into fluid:
        # Start with n_raw = (dy, -dx); flip so that it points away from centroid.
        n_raw = np.array([dy, -dx])
        n_norm = n_raw / np.linalg.norm(n_raw)

        mid = np.array([(x0 + x1) * 0.5, (y0 + y1) * 0.5])
        r_centroid = mid - centroid
        if np.dot(n_norm, r_centroid) < 0.0:
            n_norm = -n_norm

        # Pressure force on body: F = -p * n_body * ds
        dFx = -p_mid * n_norm[0] * ds
        dFy = -p_mid * n_norm[1] * ds

        Fx += dFx
        Fy += dFy

        # Moment about (x_ref, y_ref): Mz = (r x F)_z
        rx = mid[0] - x_ref
        ry = mid[1] - y_ref
        Mz += rx * dFy - ry * dFx

    # Rotate to lift/drag relative to freestream
    # Same transformation as your original compute_coeff:
    #   L = Fy * cos(alpha) - Fx * sin(alpha)
    #   D = Fy * sin(alpha) + Fx * cos(alpha)
    L = Fy * np.cos(alpha) - Fx * np.sin(alpha)
    D = Fy * np.sin(alpha) + Fx * np.cos(alpha)

    C_L = L / (q_inf * chord)
    C_D = D / (q_inf * chord)
    C_M = Mz / (q_inf * chord ** 2)

    return x_surf, y_surf, Cp_surf, C_L, C_D, C_M


# Optional: quick example usage
if __name__ == "__main__":
    # Example (replace with your actual values / file):
    # alpha in radians, Mach, T_inf, p_inf must match your run.
    vtk_file = "solution_iter973.vtk"
    Mach = 0.7
    alpha = np.deg2rad(2.0)
    T_inf = 288.0
    p_inf = 1.0e5
    chord = 1.0

    x, y, Cp, CL, CD, CM = compute_coeff_from_vtk(
        vtk_file,
        Mach_inf=Mach,
        alpha=alpha,
        T_inf=T_inf,
        p_inf=p_inf,
        chord=chord,
    )

    # Minimal check: print integrated coefficients
    print(f"C_L = {CL:.6f}, C_D = {CD:.6f}, C_M = {CM:.6f}")