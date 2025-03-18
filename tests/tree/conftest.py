"""Common fixtures for tree tests."""

import pytest

from sr2silo.tree import MatrixTree


@pytest.fixture
def basic_tree():
    """Create a basic tree with three nodes for testing."""
    tree = MatrixTree(3, ['A', 'B', 'C'])
    tree.add_edge('A', 'B', 1.0)
    tree.add_edge('A', 'C', 2.0)
    return tree


@pytest.fixture
def complex_tree():
    """Create a more complex tree for testing."""
    tree = MatrixTree(6, ['Root', 'Child1', 'Child2', 'Leaf1', 'Leaf2', 'Leaf3'])
    tree.add_edge('Root', 'Child1', 1.0)
    tree.add_edge('Root', 'Child2', 2.0)
    tree.add_edge('Child1', 'Leaf1', 1.5)
    tree.add_edge('Child1', 'Leaf2', 2.5)
    tree.add_edge('Child2', 'Leaf3', 3.0)
    return tree

