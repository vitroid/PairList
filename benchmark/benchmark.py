import pairlist as pl
from fcc import FaceCenteredCubic
from logging import getLogger, basicConfig, INFO, DEBUG
from decorator import timeit, banner
import numpy as np
from pairlist import pairs_py, pairs2_py


basicConfig(level=INFO, format="%(levelname)s %(message)s")
logger = getLogger()
logger.debug("Debug mode.")


@banner
@timeit
def crude(lattice, cell, rc=1.1):
    "Neighboring pair list by a crude double loop."
    rc2 = rc**2
    count = 0
    for i in range(len(lattice)):
        for j in range(i):
            d = lattice[i] - lattice[j]
            d -= np.floor(d + 0.5)
            d = d @ cell
            if d @ d < rc2:
                count += 1
    return count


@banner
@timeit
def numpyish(lattice, cell, rc=1.1):
    "Neighboring pair list by numpy fancy array."
    # cross-differences
    M = lattice[:, None, :] - lattice[None, :, :]
    # wrap
    M -= np.floor(M + 0.5)
    # in absolute coordinate
    M = M @ cell
    d = (M * M).sum(2)
    return d[(d < rc**2) & (0 < d)].shape[0] / 2


@banner
@timeit
def pairlist_py(lattice, cell, rc=1.1):
    "Neighboring pair list by pairlist in pure python."
    count = 0
    for i, j, d in pl.pairs_iter(
        lattice, maxdist=rc, cell=cell, _engine=(pairs_py, pairs2_py)
    ):
        count += 1
    return count


@timeit
@banner
def pairlist_c(lattice, cell, rc=1.1):
    "Neighboring pair list by pairlist in c."
    count = 0
    for i, j, d in pl.pairs_iter(lattice, maxdist=rc, cell=cell):
        count += 1
    return count


lattice, cell = FaceCenteredCubic(10)

print(crude(lattice, cell), "pairs")
print(numpyish(lattice, cell), "pairs")
print(pairlist_py(lattice, cell), "pairs")
print(pairlist_c(lattice, cell), "pairs")
