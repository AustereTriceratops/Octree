def atom_positions_from_pdb(filepath):
    atom_positions = []

    with open(filepath, "r") as file:
        line = file.readline()

        # loops over all lines in file until the end
        while line != '':
            if line[:4] == "ATOM":  # only read ATOM data
                split_line = line.split(' ')  # split by space
                # split will have to be self-implemented in C++
                entries = []

                for string in split_line:  # remove empty strings
                    if string != '':
                        entries.append(string)
                        # in C++ this is vector::push_back()

                position_vec = []
                for i in entries[6:9]:  # columns 6-8 contain XYZ coords
                    position_vec.append(float(i))

                atom_positions.append(position_vec)

            line = file.readline()

    return atom_positions