from copy import deepcopy


def encode_octant(octant):  # int from 0-7
    encoding = list(format(octant, "03b"))
    result = [0, 0, 0]

    for i in range(3):
        result[i] = int(encoding[i])

    return result


# (0,x,x) -> (1,x,x):  dim = 0, end_state = 1
# generates all possible (1,x,x) encodings

def generate_closest_octants(dim, end_state):  # dim: int in [0, 2], end_state: int in [0, 1]
    template = [0, 0, 0]
    template[dim] = (end_state + 1) % 2
    closest_octants = []

    for i in range(2):
        for j in range(2):
            octant = template.copy()

            octant[(dim + 1) % 3] = i
            octant[(dim + 2) % 3] = j

            closest_octants.append(octant)

    return closest_octants


def invert_path_at_index(path, index, dim=0):
    inverted_path = deepcopy(path)

    if index is None:
        return inverted_path

    n = len(path)
    for i in range(index, n):
        inverted_path[i][dim] = (inverted_path[i][dim] + 1) % 2

    return inverted_path


# each node is the tree has a unique path
# e.g. [[0,0,1],[1,0,1],[1,1,1]]
# as well as up to 6* neighboring nodes with similar paths
# *(up to 6 only if considering cells at the same depth or below)
#
# Our example path has neighboring paths:
# [[0,0,1],[1,0,1],[0,1,1]]  (the last octant with any of the 3 entries flipped)
# [[0,0,1],[1,1,1],[1,0,1]]  (entry [1] from the last two octants inverted)
# [[1,0,1],[0,0,1],[0,1,1]]  (sliceing the [0] indices gives a sequence
# {0, 1, 1} which is inverted to {1, 0, 0} keeping everything else fixed)
#
# So our path has 5 neighobring cells in the octree

def find_neighbors_at_same_depth(path):
    n = len(path)
    neighbors = []

    head_prev = [None, None, None]  # can just be None but I want to be consistent
    head = [None, None, None]

    for i in range(n):
        index = n - i - 1
        head_prev = head.copy()
        head = path[index].copy()

        for j in range(3):
            if head[j] != head_prev[j]:
                neighbor_path = invert_path_at_index(path, index, dim=j)
                neighbors.append(neighbor_path)

    return neighbors