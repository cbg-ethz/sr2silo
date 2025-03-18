"""Floyd-Warshall algorithm implementation for tree traversal."""

from __future__ import annotations

from typing import List, Union

import numpy as np

from sr2silo.tree.matrix_tree import MatrixTree


class FloydWarshall:
    """
    Implementation of the Floyd-Warshall algorithm to find shortest paths in a graph.

    While this is usually used for general graphs, we apply it here to traverse trees
    represented as adjacency matrices.
    """

    def __init__(self, tree: MatrixTree):
        """
        Initialize with a MatrixTree.

        Args:
            tree: MatrixTree object on which to perform the algorithm
        """
        self.tree = tree
        self.distances = np.copy(tree.get_matrix())
        self.predecessors: np.ndarray = np.full(
            (tree.num_nodes, tree.num_nodes), -1, dtype=int
        )

        # Initialize predecessor matrix
        for i in range(tree.num_nodes):
            for j in range(tree.num_nodes):
                if i != j and self.distances[i, j] < float("inf"):
                    self.predecessors[i, j] = i

    def compute_shortest_paths(self) -> None:
        """
        Execute the Floyd-Warshall algorithm to find shortest paths between all node pairs.

        This updates the internal distance and predecessor matrices.
        """
        n = self.tree.num_nodes

        # Floyd-Warshall algorithm
        for k in range(n):
            for i in range(n):
                for j in range(n):
                    if self.distances[i, k] < float("inf") and self.distances[
                        k, j
                    ] < float("inf"):
                        new_distance = self.distances[i, k] + self.distances[k, j]

                        if new_distance < self.distances[i, j]:
                            self.distances[i, j] = new_distance
                            self.predecessors[i, j] = self.predecessors[k, j]

    def get_shortest_distance(
        self, source: Union[int, str], dest: Union[int, str]
    ) -> float:
        """
        Get the shortest distance between two nodes.

        Args:
            source: Source node index or label
            dest: Destination node index or label

        Returns:
            Shortest distance between source and destination
        """
        # Convert labels to indices if needed
        src_idx = source if isinstance(source, int) else self.tree.labels.index(source)
        dst_idx = dest if isinstance(dest, int) else self.tree.labels.index(dest)

        return self.distances[src_idx, dst_idx]

    def get_shortest_path(
        self, source: Union[int, str], dest: Union[int, str]
    ) -> List[int]:
        """
        Get the shortest path between two nodes.

        Args:
            source: Source node index or label
            dest: Destination node index or label

        Returns:
            List of node indices representing the shortest path
        """
        # Convert labels to indices if needed
        src_idx = source if isinstance(source, int) else self.tree.labels.index(source)
        dst_idx = dest if isinstance(dest, int) else self.tree.labels.index(dest)

        if self.distances[src_idx, dst_idx] == float("inf"):
            return []  # No path exists

        path = [dst_idx]
        while path[0] != src_idx:
            pred = self.predecessors[src_idx, path[0]]
            if pred == -1:
                return []  # No path exists
            path.insert(0, pred)

        return path

    def print_shortest_path(
        self, source: Union[int, str], dest: Union[int, str]
    ) -> None:
        """
        Print the shortest path between two nodes with distances.

        Args:
            source: Source node index or label
            dest: Destination node index or label
        """
        # Convert labels to indices if needed
        src_idx = source if isinstance(source, int) else self.tree.labels.index(source)
        dst_idx = dest if isinstance(dest, int) else self.tree.labels.index(dest)

        path = self.get_shortest_path(src_idx, dst_idx)

        if not path:
            print(
                f"No path exists between {self.tree.labels[src_idx]} and {self.tree.labels[dst_idx]}"
            )
            return

        # Convert indices to labels for printing
        labeled_path = [self.tree.labels[idx] for idx in path]

        # Calculate segment distances
        total_distance = 0
        path_str = labeled_path[0]

        for i in range(1, len(path)):
            edge_distance = self.tree.adjacency_matrix[path[i - 1], path[i]]
            total_distance += edge_distance
            path_str += f" --({edge_distance:.1f})--> {labeled_path[i]}"

        print(f"Shortest path from {labeled_path[0]} to {labeled_path[-1]}:")
        print(f"  {path_str}")
        print(f"  Total distance: {total_distance:.1f}")

    def print_distance_matrix(self) -> None:
        """Print the distance matrix after Floyd-Warshall algorithm."""
        print("Distance matrix after Floyd-Warshall algorithm:")

        # Format header row with labels
        header = "    " + " ".join(f"{label:^5}" for label in self.tree.labels)
        print(header)

        # Format each row
        for i in range(self.tree.num_nodes):
            row = f"{self.tree.labels[i]:^4}"
            for j in range(self.tree.num_nodes):
                val = self.distances[i, j]
                if val == float("inf"):
                    row += " ∞   "
                else:
                    row += f"{val:^5.1f}"
            print(row)

    def print_all_paths(self) -> None:
        """Print all shortest paths between all pairs of nodes."""
        print("All shortest paths:")

        for i in range(self.tree.num_nodes):
            for j in range(self.tree.num_nodes):
                if i != j:
                    self.print_shortest_path(i, j)

        print()
