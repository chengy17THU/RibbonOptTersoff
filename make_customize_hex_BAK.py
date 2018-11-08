# -*- coding: utf-8 -*-


def make_hex(filename1, bond_len, filename2, r, atom_type=0, base_num=0, astyle='atomic'):
    """
    Ver 20180624: add positional astyle='atomic'
    """
    """
    Input:
        "radius", i.e. half of edge-to-edge distance of hexagon, Angstrom
        "coordinate origin", center of hexagon, Angstrom
        "base number", number of atoms in the substrate
        and,
        load original graphene data files("crude.data", in "atomic" style)
        x-axis // zigzag edge, y-axis // armchair edge
        i.e., Zigzag Graphene from VMD
    Function:
        print coordinates within defined hex into a "atomic" style data file
        re-arrange atom IDs
    """

    import numpy as np

    # load original data file
    if type(filename1) == str:
        crude = np.loadtxt(filename1, dtype=float)
        print("filename1 is str")
    elif type(filename1) == np.ndarray:
        print("filename1 is np.ndarray")
        crude = filename1
    # Finding the CO, bond length as input

    '''
    Find the "center of hex-void" nearest to (crude_x_mean, crude_y_mean),
    that's the CO to use
    '''

    # inter-facial layer
    lowest_z = min(crude[:, 4])
    facial = [crude[i, :] for i in range(crude.shape[0]) if crude[i, 4] == lowest_z]
    facial = np.array(facial)
    facial_x_mean = (facial.max(axis=0)[2] + facial.min(axis=0)[2])/2
    facial_y_mean = (facial.max(axis=0)[3] + facial.min(axis=0)[3])/2

    # Type 1
    init_x1 = facial.min(axis=0)[2] + bond_len
    init_y1 = facial.min(axis=0)[3] + bond_len*(3**0.5/2)

    co_x1 = init_x1 + bond_len*3*(round((facial_x_mean-init_x1)/(bond_len*3)))
    co_y1 = init_y1 + bond_len*(3**0.5)*(round((facial_y_mean-init_y1)/(bond_len*3**0.5)))
    # Type 2
    init_x2 = facial.min(axis=0)[2] + bond_len*2.5
    init_y2 = facial.min(axis=0)[3] + bond_len*(3**0.5)

    co_x2 = init_x2 + bond_len*3*(round((facial_x_mean-init_x2)/(bond_len*3)))
    co_y2 = init_y2 + bond_len*(3**0.5)*(round((facial_y_mean-init_y2)/(bond_len*3**0.5)))

    if abs(co_x1 - facial_x_mean) < abs(co_x2 - facial_x_mean):
        co_x = co_x1
        co_y = co_y1
    else:
        co_x = co_x2
        co_y = co_y2

    print("CO of hex slider is: ", co_x, co_y)
    CO = [co_x, co_y]

    # define Inputs:
    # r = 19
    # Choose VERY CAREFULLY and REPEATEDLY to avoid extra atoms on EDGES
    if not base_num:
        base_num = int(crude.min(axis=0)[0])

    # boundary: 6 lines
    k1 = 3**0.5/3
    k2 = -3**0.5/3

    i = 0
    j = 0
    data = np.zeros(crude.shape, dtype=float)
    while i < crude.shape[0]:
        if abs(co_x - crude[i, 2]) < r:
            r_x = crude[i, 2] - co_x
            r_y = crude[i, 3] - co_y
            dist1 = abs(r_y-r_x*k1)/(1+k1**2)**0.5  # distance to "y = k1 * x"
            dist2 = abs(r_y-r_x*k2)/(1+k2**2)**0.5  # distance to "y = k2 * x"
            if (dist1 < r) & (dist2 < r):
                data[j, :] = crude[i, :]
                data[j, 0] = base_num + j
                j = j + 1
        i = i + 1

    f = open(filename2, 'w')

    j = 0
    if astyle == 'atomic':
        # atom ID, atom type, x, y, z
        while int(data[j, 0]) != 0:
            if atom_type:
                data[j, 1] = atom_type
            line = tuple(data[j, :])
            print("%d   %d   %.6f   %.6f   %.6f" % line, file=f)
            j = j + 1
    elif astyle == 'full':
        print('Writing in style \"full\"')
        # atom ID, molecule ID, atom type, charge, x, y, z
        while int(data[j, 0]) != 0:
            if atom_type:
                data[j, 1] = atom_type
            line = tuple([data[j, 0], data[j, 1], data[j, 1], 0, data[j, 2], data[j, 3], data[j, 4]])
            print("%d   %d   %d   %.6f   %.6f   %.6f   %.6f" % line, file=f)
            j = j + 1
    f.close()

    tail_ID = base_num + j - 1
    print("Tail ID is %d, height is %f" % (tail_ID, data[0, 4]))
    return tail_ID, data, CO


# Test this function
# make_hex('crude.data', 1.418, 'make_hex.dat', 19, 2)

def make_O(spacing, xbounds, ybounds, z, init_ID, atom_type, filename, astyle='atomic'):
    """
    Ver 20180625: add positional astyle='atomic'
    """
    init_x = xbounds[0] + spacing/2
    init_y = ybounds[0] + spacing/2

    mut_x = init_x
    mut_y = init_y
    mut_ID = init_ID

    f = open(filename, 'w')
    if astyle == 'atomic':
        while mut_x < xbounds[1]:
            while mut_y < ybounds[1]:
                line = (mut_ID, atom_type, mut_x, mut_y, z)
                print("%d   %d   %.6f   %.6f   %.6f" % line, file=f)
                mut_y = mut_y + spacing
                mut_ID = mut_ID + 1
            mut_y = init_y
            mut_x = mut_x + spacing
    elif astyle == 'full':
        m_ID = atom_type
        while mut_x < xbounds[1]:
            while mut_y < ybounds[1]:
                m_ID += 1
                line = (mut_ID, m_ID, atom_type, 0.0, mut_x, mut_y, z)
                print("%d   %d   %d   %.6f   %.6f   %.6f   %.6f" % line, file=f)
                mut_y = mut_y + spacing
                mut_ID = mut_ID + 1
            mut_y = init_y
            mut_x = mut_x + spacing

    print("Tail ID is %d" % (mut_ID - 1))
    return 0
