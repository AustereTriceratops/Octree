from octree import *

# size as exponent of 2, deg as number of degenerate points
def generate_test_points(size, deg):
    num = 2 ** size
    points = []

    for i in range(num):
        for j in range(num):
            for k in range(num):
                for _ in range(deg):
                    points.append([i, j, k])

    return points


def test_node_counting():
    points_array = [generate_test_points(i, 1) for i in range(3)]
    values = [1, 9, 73]

    for i in range(3):
        tree = Octree(points_array[i], 1).tree

        assert tree.count_subnodes() == values[i]

    print('counting test passed')


def test_classification(deg):
    points = generate_test_points(2, deg)

    tree = Octree(points, deg, scale=2.01).tree

    assert len(tree.leaves) == 120
    assert len(tree.exterior_cells) == 56
    assert len(tree.boundary_cells) == 56

    print("classification test passed")


if __name__ == '__main__':
    print("testing")
    test_classification(1)
    test_node_counting()
