"""Tree representation using an adjacency matrix."""

from __future__ import annotations

import numpy as np
from typing import Dict, List, Optional, Set, Tuple, Union


class MatrixTree:
    """
    A class to represent a tree structure using an adjacency matrix.

    This implementation uses a 2D numpy array to represent connections between nodes.
    Each entry [i,j] represents the weight/distance from node i to node j.
    """

    def __init__(self, num_nodes: int, labels: Optional[List[str]] = None):
        """
        Initialize an empty tree with specified number of nodes.

        Args:
            num_nodes: Number of nodes in the tree
            labels: Optional labels for the nodes. If not provided, nodes will be labeled 0...n-1
        """
        # Initialize adjacency matrix with infinity (representing no direct connection)
        # numpy.inf represents infinity in floating point
        self.adjacency_matrix = np.full((num_nodes, num_nodes), float('inf'))

        # Set diagonal to 0 (distance from node to itself is 0)
        np.fill_diagonal(self.adjacency_matrix, 0)

        self.num_nodes = num_nodes

        # Use provided labels or create default ones (0...n-1)
        if labels and len(labels) == num_nodes:
            self.labels = labels
        else:
            self.labels = [str(i) for i in range(num_nodes)]

        # Keep track of parent-child relationships to maintain tree structure
        self.parent: Dict[int, Optional[int]] = {i: None for i in range(num_nodes)}
        self.children: Dict[int, Set[int]] = {i: set() for i in range(num_nodes)}

    def add_edge(self, parent: Union[int, str], child: Union[int, str], weight: float = 1.0) -> bool:
        """
        Add an edge between parent and child nodes with the given weight.

        Args:
            parent: Parent node index or label
            child: Child node index or label
            weight: Edge weight (default: 1.0)

        Returns:
            True if edge added successfully, False otherwise
        """
        # Convert labels to indices if needed
        parent_idx = parent if isinstance(parent, int) else self.labels.index(parent)
        child_idx = child if isinstance(child, int) else self.labels.index(child)

        # Check if edge would create a cycle
        if self._would_create_cycle(parent_idx, child_idx):
            return False

        # Add edge in both directions (undirected graph)
        self.adjacency_matrix[parent_idx, child_idx] = weight
        self.adjacency_matrix[child_idx, parent_idx] = weight

        # Update parent-child relationships
        if self.parent[child_idx] is not None:
            # Remove from previous parent's children
            self.children[self.parent[child_idx]].remove(child_idx)

        self.parent[child_idx] = parent_idx
        self.children[parent_idx].add(child_idx)

        return True

    def _would_create_cycle(self, parent_idx: int, child_idx: int) -> bool:
        """
        Check if adding an edge would create a cycle.

        Args:
            parent_idx: Parent node index
            child_idx: Child node index

        Returns:
            True if adding edge would create a cycle, False otherwise
        """
        # If child is already an ancestor of parent, adding this edge would create a cycle
        current = parent_idx
        while current is not None:
            if current == child_idx:
                return True
            current = self.parent[current]
        return False

    def get_matrix(self) -> np.ndarray:
        """
        Get the adjacency matrix representation of the tree.

        Returns:
            Adjacency matrix as numpy array
        """
        return self.adjacency_matrix

    def get_neighbors(self, node: Union[int, str]) -> List[Tuple[int, float]]:
        """
        Get all neighbors of a node with their edge weights.

        Args:
            node: Node index or label

        Returns:
            List of (neighbor_idx, weight) tuples
        """
        # Convert label to index if needed
        node_idx = node if isinstance(node, int) else self.labels.index(node)

        neighbors = []
        for j in range(self.num_nodes):
            if self.adjacency_matrix[node_idx, j] not in [0, float('inf')]:
                neighbors.append((j, self.adjacency_matrix[node_idx, j]))

        return neighbors

    def print_tree(self, root: Optional[Union[int, str]] = None, indent: str = "") -> None:
        """
        Print the tree structure starting from the specified root.

        Args:
            root: Root node to start printing from (default: find a root)
            indent: String used for indentation (used for recursion)
        """
        # Find a root if not specified
        if root is None:
            for i in range(self.num_nodes):
                if self.parent[i] is None:
                    root = i
                    break
            if root is None and self.num_nodes > 0:
                # If no root found, use node 0
                root = 0

        # Convert label to index if needed
        root_idx = root if isinstance(root, int) else self.labels.index(root)

        # Print current node
        print(f"{indent}{self.labels[root_idx]}")

        # Print children
        for child in self.children[root_idx]:
            weight = self.adjacency_matrix[root_idx, child]
            print(f"{indent}├─({weight})─", end="")
            self.print_tree(child, indent + "│  ")

    def __str__(self) -> str:
        """String representation showing the adjacency matrix."""
        result = "MatrixTree with adjacency matrix:\n"

        # Format header row with labels
        header = "    " + " ".join(f"{label:^5}" for label in self.labels)
        result += header + "\n"

        # Format each row
        for i in range(self.num_nodes):
            row = f"{self.labels[i]:^4}"
            for j in range(self.num_nodes):
                val = self.adjacency_matrix[i, j]
                if val == float('inf'):
                    row += " ∞   "
                else:
                    row += f"{val:^5.1f}"
            result += row + "\n"

        return result

