"""Tests for the MatrixTree class."""

import numpy as np
import pytest

from sr2silo.tree import MatrixTree


def test_matrix_tree_initialization():
    """Test that MatrixTree is initialized correctly."""
    # Test with default labels
    tree = MatrixTree(5)
    assert tree.num_nodes == 5
    assert tree.labels == ['0', '1', '2', '3', '4']
    assert np.array_equal(np.diag(tree.adjacency_matrix), np.zeros(5))

    # Test with custom labels
    custom_labels = ['A', 'B', 'C', 'D', 'E']
    tree = MatrixTree(5, labels=custom_labels)
    assert tree.labels == custom_labels


def test_add_edge():
    """Test adding edges to the tree."""
    tree = MatrixTree(4, ['A', 'B', 'C', 'D'])

    # Add edges and verify they're added correctly
    assert tree.add_edge('A', 'B', 2.0) == True
    assert tree.adjacency_matrix[0, 1] == 2.0
    assert tree.adjacency_matrix[1, 0] == 2.0  # Undirected graph

    assert tree.add_edge(1, 2) == True  # Using indices
    assert tree.adjacency_matrix[1, 2] == 1.0  # Default weight is 1.0

    # Verify parent-child relationships
    assert tree.parent[1] == 0  # B's parent is A
    assert tree.parent[2] == 1  # C's parent is B
    assert 1 in tree.children[0]  # B is child of A


def test_prevent_cycles():
    """Test that adding edges that create cycles is prevented."""
    tree = MatrixTree(3, ['A', 'B', 'C'])

    # Create a path A -> B -> C
    tree.add_edge('A', 'B')
    tree.add_edge('B', 'C')

    # Attempting to create C -> A should fail (would create cycle)
    assert tree.add_edge('C', 'A') == False


def test_get_neighbors():
    """Test retrieving neighbors of a node."""
    tree = MatrixTree(4, ['A', 'B', 'C', 'D'])
    tree.add_edge('A', 'B', 1.5)
    tree.add_edge('A', 'C', 2.0)
    tree.add_edge('B', 'D', 0.5)

    # Get neighbors by label
    neighbors = tree.get_neighbors('A')
    assert len(neighbors) == 2
    assert (1, 1.5) in neighbors  # B with weight 1.5
    assert (2, 2.0) in neighbors  # C with weight 2.0

    # Get neighbors by index
    neighbors = tree.get_neighbors(1)  # B
    assert len(neighbors) == 2
    assert (0, 1.5) in neighbors  # A with weight 1.5
    assert (3, 0.5) in neighbors  # D with weight 0.5


def test_print_tree(capsys):
    """Test tree printing functionality."""
    tree = MatrixTree(3, ['A', 'B', 'C'])
    tree.add_edge('A', 'B', 2.0)
    tree.add_edge('A', 'C', 3.0)

    tree.print_tree()
    captured = capsys.readouterr()

    # Check if output contains node labels
    assert 'A' in captured.out
    assert 'B' in captured.out
    assert 'C' in captured.out

    # Check if weights are displayed
    assert '(2.0)' in captured.out
    assert '(3.0)' in captured.out


def test_string_representation():
    """Test string representation of the tree."""
    tree = MatrixTree(3, ['A', 'B', 'C'])
    tree.add_edge('A', 'B', 1.0)

    str_rep = str(tree)

    # Check if the representation contains the expected elements
    assert 'MatrixTree' in str_rep
    assert 'A' in str_rep
    assert 'B' in str_rep
    assert 'C' in str_rep
    assert '1.0' in str_rep
    assert '∞' in str_rep  # Should contain infinity symbol for disconnected nodes

