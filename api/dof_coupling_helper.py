from collections import defaultdict, deque

class DOFCouplingHelper:
    """
    Helper class for building coupling constraints from a VTK unstructured grid.
    Interfaces with the main APEC_M class to find connected components and generate constraints.
    """

    @staticmethod
    def build_dof_coupling_from_vtk(grid, array_name="CouplingID", master_node=None):
        """
        Builds coupling constraints from a VTK unstructured grid based on a point data array.
        Allows imposing a master node for each connected component.

        Parameters:
            grid (vtk.vtkUnstructuredGrid): The mesh to process.
            array_name (str): The name of the point data array containing coupling IDs (default "CouplingID").
            master_node (int, optional): Node ID to use as master for its connected component. If not in component, defaults to smallest ID.

        Returns:
            constraints (list of tuple): List of (slave, master) pairs representing constraints u[slave] - u[master] = 0.
        """
        # --- get coupling ID array ---
        cid_array = grid.GetPointData().GetArray(array_name)
        if cid_array is None:
            raise ValueError(f"PointData array '{array_name}' not found.")

        npoints = grid.GetNumberOfPoints()
        cids = [int(cid_array.GetValue(i)) for i in range(npoints)]

        # --- build adjacency graph for nodes with cid=1 ---
        adj = defaultdict(list)

        for i in range(grid.GetNumberOfCells()):
            cell = grid.GetCell(i)
            ids = [cell.GetPointId(j) for j in range(cell.GetNumberOfPoints())]
            # only keep nodes with cid=1
            active = [nid for nid in ids if cids[nid] == 1]
            for a in range(len(active)):
                for b in range(a+1, len(active)):
                    u, v = active[a], active[b]
                    adj[u].append(v)
                    adj[v].append(u)

        # --- find connected components ---
        visited = set()
        components = []

        for nid in range(npoints):
            if cids[nid] == 1 and nid not in visited:
                comp = []
                queue = deque([nid])
                visited.add(nid)
                while queue:
                    u = queue.popleft()
                    comp.append(u)
                    for v in adj[u]:
                        if v not in visited:
                            visited.add(v)
                            queue.append(v)
                components.append(comp)

        # --- build constraints for each component ---
        constraints = []
        for comp in components:
            # Use master_node if provided and in component, else use smallest ID
            if master_node is not None and master_node in comp:
                master = master_node
            else:
                master = min(comp)
            comp_set = set(comp)
            seen = {master}
            queue = deque([master])
            while queue:
                u = queue.popleft()
                for v in adj[u]:
                    if v in comp_set and v not in seen:
                        constraints.append((v, u))  # u[v] - u[u] = 0
                        seen.add(v)
                        queue.append(v)

        return constraints

    @staticmethod
    def print_constraints(constraints, master_node=None):
        """
        Prints the master node and all master-slave pairs in the constraints list in a readable format.

        Parameters:
            constraints (list of tuple): List of (slave, master) pairs.
            master_node (int, optional): The master node ID, if known.
        """
        if not constraints:
            print("No constraints found.")
            return

        masters = set(master for _, master in constraints)
        slaves = set(slave for slave, _ in constraints)
        all_nodes = masters | slaves

        if master_node is not None:
            print(f"Master node: {master_node}")
        else:
            # Guess master node(s) as those that appear only as master, never as slave
            possible_masters = masters - slaves
            if possible_masters:
                print(f"Master node(s): {', '.join(str(m) for m in sorted(possible_masters))}")
            else:
                print("Master node(s): Could not be determined.")

        print("Constraints (slave, master):")
        for slave, master in constraints:
            print(f"u[{slave}] - u[{master}] = 0")