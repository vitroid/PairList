import pairlist as pl
from fcc import FaceCenteredCubic
import time
import matplotlib.pyplot as plt
import numpy as np


def pairlist_c(lattice, cell, rc=1.1):
    "Neighboring pair list by pairlist in c."
    count = 0
    for i, j, d in pl.pairs_iter(lattice, maxdist=rc, cell=cell):
        count += 1
    return count


nvertex = []
durations = []
for i in (2, 3, 4, 6, 8, 12, 16, 24, 32, 48, 64):
    lattice, cell = FaceCenteredCubic(i)
    nvertex.append(len(lattice))
    start = time.time()
    pairlist_c(lattice, cell)
    duration = time.time() - start
    durations.append(duration)

nvertex = np.array(nvertex, dtype=float)
plt.loglog(nvertex, durations, "o-")
plt.loglog(nvertex, nvertex / 100000, "-")
plt.xlabel("# of vertices")
plt.ylabel("Time / s")

plt.show()
