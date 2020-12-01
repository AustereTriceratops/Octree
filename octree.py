from path import *

from copy import deepcopy


class Octree:
    def __init__(self, atom_positions, atoms_per_cell, scale=1.0):
        self.tree = Node(atom_positions)
        self.atoms_per_cell = atoms_per_cell
        self.scale = scale

        self.tree.scale_bounding_box(scale)

        self.build_octree(atoms_per_cell)

    def build_octree(self, atoms_per_cell):
        tree = self.tree
        tree.span(atoms_per_cell)
        tree.collect_leaf_nodes()
        tree.index_leaf_nodes()
        tree.find_connectivity()

        tree.classify_exterior()
        tree.classify_boundary()
        tree.classify_interior()


class Node:
    def __init__(self, points, parent=None, octant=None):
        self.points = points
        self.parent = parent
        self.octant = octant
        self.children = []
        self.path = []

        self.dimensions = None
        self.center = None
        self.depth = 0

        # unsure if this should be exclusive to root node or include all internal nodes
        self.leaves = []
        self.exterior_cells = []
        self.boundary_cells = []
        self.interior_cells = []

        # leaf node data
        self.index = None
        self.category = None
        self.neighbor_indices = []

        # this initializes the root node
        if parent is None:
            self.initialize_root()

        # initialize child nodes
        else:
            self.initialize_child()

    # ========= INITIALIZATION SUBROUTINES =========

    def initialize_root(self):
        minima_position, maxima_position = None, None

        # should always be true for root node, but just in case
        if self.points:
            minima_position = find_minima(self.points)
            maxima_position = find_maxima(self.points)

            self.dimensions = [maxima_position[i] - minima_position[i] for i in range(3)]
            self.center = [minima_position[i] + self.dimensions[i] / 2.0 for i in range(3)]

    def initialize_child(self):
        self.depth = self.parent.depth + 1
        self.dimensions = [self.parent.dimensions[i] / 2.0 for i in range(3)]

        octant_encoding = encode_octant(self.octant)  # convert octant to binary
        self.path = deepcopy(self.parent.path)
        self.path.append(octant_encoding)

        offset = [self.dimensions[i] / 2.0 for i in range(3)]

        self.center = [self.parent.center[i] + (1 - 2 * octant_encoding[i]) * offset[i] for i in range(3)]

    # ========= RECURSIVE METHODS =========

    # build octree with recursion
    def span(self, atoms_per_cell):
        self.create_children(atoms_per_cell)

        for child in self.children:
            child.span(atoms_per_cell)

    # returns total number of nodes
    def count_subnodes(self):
        count = 1

        for child in self.children:
            count += child.count_subnodes()

        return count

    def count_leaf_nodes(self):
        count = 0

        if not self.children:
            return 1

        for child in self.children:
            count += child.count_leaf_nodes()

        return count

    def collect_leaf_nodes(self):
        leaves = []

        if not self.children:
            leaves.append(self)
            return leaves

        for child in self.children:
            new_leaves = child.collect_leaf_nodes()
            leaves.extend(new_leaves)

        if self.depth == 0:
            self.leaves = leaves

        return leaves

    # ========= PATH METHODS =========

    # returns a descendant node with a unique path relative to this node
    # If there is no node with that path in the tree, it terminates
    # at the latest leaf node
    def follow_encoded_path(self, path, start=0):
        if self.children:
            n = len(path)

            if start < n:
                enc = path[start]
                ind = 4 * enc[0] + 2 * enc[1] + enc[2]

                return self.children[ind].follow_encoded_path(path, start=start + 1)
            else:
                return self
        else:
            return self

    def follow_path(self, path, start=0):
        if self.children:
            ind = path[start]
            return self.children[ind].follow_path(path, start=start + 1)
        else:
            return self

    # ========= OTHER METHODS =========

    def scale_bounding_box(self, fac):
        for i in range(3):
            self.dimensions[i] *= fac

    # generates 8 child nodes according to the octree algorithm.
    # this will only execute if the number of atoms in this cell
    # is greater than the minimum number of atoms allowed in a cell

    def create_children(self, atoms_per_cell):

        # exit function early if there are too few atoms in the cell
        if len(self.points) <= atoms_per_cell:
            self.children = []
            return

        # divide into octants
        divided_points = [[] for _ in range(8)]

        for i in range(len(self.points)):
            p = self.points[i]
            x_pol = (p[0] <= self.center[0])
            y_pol = (p[1] <= self.center[1])
            z_pol = (p[2] <= self.center[2])

            # abuse of Python's type coercion
            # basically take the sequence of x,y,z_pol as
            # a binary number and convert to base 10
            ind = x_pol * 4 + y_pol * 2 + z_pol

            divided_points[ind].append(p)

        for i in range(8):
            self.children.append(Node(divided_points[i], parent=self, octant=i))

    def index_leaf_nodes(self):
        n = len(self.leaves)

        for i in range(n):
            self.leaves[i].index = i

    def leaf_indices_along(self, dim, side):
        octants = generate_closest_octants(dim, side)
        leaf_indices = []

        if not self.children:
            return [self.index]

        for i in range(4):  # 4 octants to check
            ind = 4 * octants[i][0] + 2 * octants[i][1] + octants[i][2]
            inds = self.children[ind].leaf_indices_along(dim, side)
            leaf_indices.extend(inds)

        return leaf_indices

    def find_connectivity(self):
        for i in range(len(self.leaves)):
            leaf = self.leaves[i]
            leaf.neighbor_indices = []

            path = leaf.path
            neighbor_paths = find_neighbors_at_same_depth(path)

            for j in range(len(neighbor_paths)):
                neighbor_path = neighbor_paths[j]
                neighbor = self.follow_encoded_path(neighbor_path)

                if not neighbor.children:  # if a leaf node
                    leaf.neighbor_indices.append(neighbor.index)
                else:
                    main_octant = path[-1]
                    other_octant = neighbor_path[-1]

                    dim = None
                    end_state = None

                    for k in range(3):
                        if main_octant[k] != other_octant[k]:
                            dim = k
                            end_state = other_octant[k]
                            break

                    neighbor_indices = neighbor.leaf_indices_along(dim, end_state)
                    leaf.neighbor_indices.extend(neighbor_indices)

    # ========= CLASSIFICATION =========

    def classify_exterior(self):
        self.exterior_cells = []

        for leaf in self.leaves:
            if not leaf.points:
                leaf.category = 'exterior'
                self.exterior_cells.append(leaf)

    # loop through the empty cells to find all neighboring cells
    # classify neighboring cells appropriately
    def classify_boundary(self):  # this is a mess
        self.boundary_cells = []

        for i in range(len(self.exterior_cells)):
            leaf = self.exterior_cells[i]

            for index in leaf.neighbor_indices:
                neighbor = self.leaves[index]
                if neighbor.points and neighbor.category != 'boundary':
                    neighbor.category = 'boundary'
                    self.boundary_cells.append(neighbor)

    def classify_interior(self):
        self.interior_cells = []

        for leaf in self.leaves:
            if leaf.category != 'exterior' and leaf.category != 'boundary':
                leaf.category = 'interior'
                self.interior_cells.append(leaf)


# ============ FUNCTIONS ============
# not done with np because I plan to move this to C++

def find_minima(atom_centers):
    minima_position = [atom_centers[0][i] for i in range(3)]

    for pos in atom_centers:
        for i in range(3):
            if pos[i] < minima_position[i]:
                minima_position[i] = pos[i]
    return minima_position


def find_maxima(atom_centers):
    maxima_position = [atom_centers[0][i] for i in range(3)]

    for pos in atom_centers:
        for i in range(3):
            if pos[i] > maxima_position[i]:
                maxima_position[i] = pos[i]
    return maxima_position

