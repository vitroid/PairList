import pairlist as pl
import numpy as np


def SimpleCubic(N):
    cell = np.eye(3) * N
    lattice = np.array([(x / N, y / N, z / N)
                        for x in range(N)
                        for y in range(N)
                        for z in range(N)])
    return lattice, cell


N = 2
lattice, cell = SimpleCubic(N)
# bipartile BCC lattice
lattice2 = lattice + 1 / (N * 2)

for i, j in pl.pairs_iter(lattice, rc=1.1 * 3**0.5 / 2,
                          cell=cell, pos2=lattice2, distance=False):
    print(i, j)
