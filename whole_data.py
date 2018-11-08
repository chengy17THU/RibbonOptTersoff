# -*- coding: utf-8 -*-


import numpy as np
from matplotlib import pyplot as plt

from make_customize_hex import make_hex, make_O

'''
C-C = 1.438 for Tersoff, June 29, 2018
in 'FULL' style
NO FUCKING O INCLUSION THIS TIME.
'''
'''
A 6-layer A-B stacking graphene sheet is needed
    x-axis // armchair edge, y-axis // zigzag edge
    i.e., Zigzag Graphene from VMD
'''

lo = np.loadtxt('./lo.dat', dtype=float)
hi = np.loadtxt('./hi.dat', dtype=float)

bond_len = 1.43876

f = open('./lo_cut.dat', 'w')
i = 0
while i < lo.shape[0]:
    if lo[i, 0] <= 3936:
        lo[i, 1] = 1
    elif lo[i, 0] <= 7872:
        lo[i, 1] = 2
        lo[i, 4] = 3.34
    elif lo[i, 0] <= 11808:
        lo[i, 1] = 3
        lo[i, 4] = 3.34*2
    line = tuple([lo[i, 0], lo[i, 1], lo[i, 1], 0, lo[i, 2], lo[i, 3], lo[i, 4]])
    print("%d   %d   %d   %.6f   %.6f   %.6f   %.6f" % line, file=f)
    i = i + 1

f.close()

i = 0
while i < hi.shape[0]:
    if hi[i, 0] <= 15744:
        hi[i, 1] = 4
        hi[i, 4] = 3.34*3
    elif hi[i, 0] <= 19680:
        hi[i, 1] = 5
        hi[i, 4] = 3.34*4
    elif hi[i, 0] <= 23616:
        hi[i, 1] = 6
        hi[i, 4] = 3.34*5
    i += 1

l4 = [hi[i, :] for i in range(hi.shape[0]) if hi[i, 1] == 4]
l4 = np.array(l4)
l5 = [hi[i, :] for i in range(hi.shape[0]) if hi[i, 1] == 5]
l5 = np.array(l5)
l6 = [hi[i, :] for i in range(hi.shape[0]) if hi[i, 1] == 6]
l6 = np.array(l6)

print(hi.shape, l4.shape, l5.shape, l6.shape)
r = int((hi.max(axis=0)[2] - hi.min(axis=0)[2])/4)
tail_l4, l4, CO4 = make_hex(l4, bond_len, 'l4_cut.dat', r+7, base_num=0, astyle='full')
tail_l5, l5, CO5 = make_hex(l5, bond_len, 'l5_cut.dat', r+7, base_num=tail_l4+1, astyle='full')
tail_l6, l6, CO6 = make_hex(l6, bond_len, 'l6_cut.dat', r+7, base_num=tail_l5+1, astyle='full')

# plot l4/type4
ax1 = plt.subplot(1, 3, 1)
ax1.set_aspect('equal')
plt.scatter(l4[:, 2], l4[:, 3], s=3)
plt.scatter(CO4[0], CO4[1], s=5)

ax2 = plt.subplot(1, 3, 2)
ax2.set_aspect('equal')
plt.scatter(l5[:, 2], l5[:, 3], s=3)
plt.scatter(CO5[0], CO5[1], s=5)

ax3 = plt.subplot(1, 3, 3)
ax3.set_aspect('equal')
plt.scatter(l6[:, 2], l6[:, 3], s=3)
plt.scatter(CO6[0], CO6[1], s=5)

plt.savefig('viewCO.png', dpi=600)
plt.close(0)

# Calculate box bounds
lay1 = [lo[i, :] for i in range(lo.shape[0]) if lo[i, 4] == 0.0]
lay1 = np.array(lay1)
x_lo, y_lo = lay1.min(axis=0)[2: 4]
x_hi = lay1.max(axis=0)[2] + bond_len
y_hi = lay1.max(axis=0)[3] + bond_len*(3**0.5)/2

xbounds = [x_lo, x_hi]
ybounds = [y_lo, y_hi]
print(xbounds, ybounds)

# Form the O array
# spacing = 5  # x and y spacing should be the same
# make_O(spacing, xbounds, ybounds, 3.34*3, ID_type2 + 1, 7, 'o.dat', astyle='full')
