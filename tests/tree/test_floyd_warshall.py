"""Tests for the FloydWarshall class."""

import numpy as np
import pytest

from sr2silo.tree import MatrixTree, FloydWarshall


@pytest.fixture
def sample_tree():
    """Create a sample tree for testing."""
    tree = MatrixTree(5, ['A', 'B', 'C', 'D', 'E'])
    tree.add_edge('A', 'B', 2.0)
    tree.add_edge('A', 'C', 3.0)
    tree.add_edge('B', 'D', 1.0)
    tree.add_edge('C', 'E', 2.0)
    return tree


@pytest.fixture
def computed_floyd_warshall(sample_tree):
    """Create and compute FloydWarshall for the sample tree."""
    fw = FloydWarshall(sample_tree)
    fw.compute_shortest_paths()
    return fw


def test_floyd_warshall_initialization(sample_tree):
    """Test initialization of FloydWarshall."""
    fw = FloydWarshall(sample_tree)

    # Check that distances are initially copied from the tree
    assert np.array_equal(fw.distances, sample_tree.adjacency_matrix)

    # Check predecessors initialization
    assert fw.predecessors[0, 1] == 0  # A->B has A as predecessor
    assert fw.predecessors[0, 2] == 0  # A->C has A as predecessor
    assert fw.predecessors[1, 0] == 1  # B->A has B as predecessor


def test_compute_shortest_paths(computed_floyd_warshall, sample_tree):
    """Test computation of shortest paths."""
    fw = computed_floyd_warshall

    # Direct connections should remain the same
    assert fw.distances[0, 1] == 2.0  # A->B
    assert fw.distances[0, 2] == 3.0  # A->C

    # Transitive connections should be computed
    assert fw.distances[0, 3] == 3.0  # A->B->D
    assert fw.distances[0, 4] == 5.0  # A->C->E
    assert fw.distances[1, 4] == 7.0  # B->A->C->E


def test_get_shortest_distance(computed_floyd_warshall):
    """Test retrieving shortest distances."""
    fw = computed_floyd_warshall

    # Test using indices
    assert fw.get_shortest_distance(0, 3) == 3.0  # A->D

    # Test using labels
    assert fw.get_shortest_distance('A', 'D') == 3.0
    assert fw.get_shortest_distance('B', 'E') == 7.0


def test_get_shortest_path(computed_floyd_warshall):
    """Test retrieving shortest paths."""
    fw = computed_floyd_warshall

    # Path from A to D should be [A, B, D]
    path = fw.get_shortest_path('A', 'D')
    assert [fw.tree.labels[idx] for idx in path] == ['A', 'B', 'D']

    # Path from B to E should be [B, A, C, E]
    path = fw.get_shortest_path('B', 'E')
    assert [fw.tree.labels[idx] for idx in path] == ['B', 'A', 'C', 'E']


def test_print_shortest_path(computed_floyd_warshall, capsys):
    """Test printing shortest path."""
    fw = computed_floyd_warshall

    fw.print_shortest_path('A', 'E')
    captured = capsys.readouterr()

    # Check if output contains the expected elements
    assert 'A' in captured.out
    assert 'C' in captured.out
    assert 'E' in captured.out
    assert '5.0' in captured.out  # Total distance


def test_print_distance_matrix(computed_floyd_warshall, capsys):
    """Test printing distance matrix."""
    fw = computed_floyd_warshall

    fw.print_distance_matrix()
    captured = capsys.readouterr()

    # Check if output contains the expected elements
    for label in ['A', 'B', 'C', 'D', 'E']:
        assert label in captured.out

    # Check if some distances are shown
    assert '2.0' in captured.out
    assert '3.0' in captured.out
    assert '5.0' in captured.out


def test_print_all_paths(computed_floyd_warshall, capsys):
    """Test printing all paths."""
    fw = computed_floyd_warshall

    fw.print_all_paths()
    captured = capsys.readouterr()

    # Check that outputs for all node pairs are included
    for source in ['A', 'B', 'C', 'D', 'E']:
        for dest in ['A', 'B', 'C', 'D', 'E']:
            if source != dest:
                assert f"from {source} to {dest}" in captured.out or \
                       f"No path exists between {source} and {dest}" in captured.out

