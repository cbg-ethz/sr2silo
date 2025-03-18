"""Examples demonstrating the tree representation and traversal."""

from __future__ import annotations

from sr2silo.tree import FloydWarshall, MatrixTree


def example_family_tree():
    """Create and analyze a family tree example."""
    # Create a family tree
    print("Creating a family tree...")
    tree = MatrixTree(
        7, ["Grandparent", "Parent1", "Parent2", "Child1", "Child2", "Child3", "Child4"]
    )

    # Add edges (relationships)
    tree.add_edge("Grandparent", "Parent1", 1.0)
    tree.add_edge("Grandparent", "Parent2", 1.0)
    tree.add_edge("Parent1", "Child1", 1.0)
    tree.add_edge("Parent1", "Child2", 1.0)
    tree.add_edge("Parent2", "Child3", 1.0)
    tree.add_edge("Parent2", "Child4", 1.0)

    # Print the tree structure
    print("\nTree structure:")
    tree.print_tree()

    # Print the adjacency matrix
    print("\n" + str(tree))

    # Calculate shortest paths
    fw = FloydWarshall(tree)
    fw.compute_shortest_paths()

    # Print the distance matrix
    print()
    fw.print_distance_matrix()

    # Print some specific paths
    print()
    fw.print_shortest_path("Child1", "Child4")
    fw.print_shortest_path("Child2", "Child3")

    return tree, fw


def example_weighted_tree():
    """Create and analyze a weighted tree example."""
    print("Creating a weighted tree...")
    tree = MatrixTree(5, ["A", "B", "C", "D", "E"])

    # Add weighted edges
    tree.add_edge("A", "B", 2.5)
    tree.add_edge("A", "C", 1.5)
    tree.add_edge("B", "D", 3.0)
    tree.add_edge("C", "E", 2.0)

    # Print the tree structure
    print("\nTree structure:")
    tree.print_tree()

    # Calculate shortest paths
    fw = FloydWarshall(tree)
    fw.compute_shortest_paths()

    # Print all shortest paths
    print()
    fw.print_all_paths()

    return tree, fw


if __name__ == "__main__":
    example_family_tree()
    print("\n" + "=" * 50 + "\n")
    example_weighted_tree()
