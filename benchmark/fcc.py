import pairlist as pl
import numpy as np


def FaceCenteredCubic(N):
    N2 = N * 2
    ea = np.array([2**0.5, 0, 0])
    eb = np.array([0, 2**0.5, 0])
    ec = np.array([0, 0, 2**0.5])
    cell = np.array([ea, eb, ec])
    lattice = np.array([(x / N2, y / N2, z / N2)
                        for x in range(N2)
                        for y in range(N2)
                        for z in range(N2)
                        if (x + y + z) % 2 == 0])
    return lattice, cell * N
