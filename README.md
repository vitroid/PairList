# PairList
Generates the pair list of atoms that are closer to each other than the
given threshold under the periodic boundary conditions.

version 0.2.12.4

## Usage

See `pairlist.h` for the function definition and `pairlist-test.c` for usage.

Python API is served in pairlist.py. Here is a sample code to use it.

<!-- ```python

```

```python

``` -->

<!-- ## Benchmark tests -->

To find the neighbors in a face-centered cubic lattice of size 10x10x10 on a MacBook Air 2021 (Apple Silicon),

```shell
$ python benchmark.py
INFO crude: Neighboring pair list by a crude double loop.
INFO crude: 18024 ms
INFO crude: end.
24000 pairs
INFO numpyish: Neighboring pair list by numpy fancy array.
INFO numpyish: 741 ms
INFO numpyish: end.
24000.0 pairs
INFO pairlist_py: Neighboring pair list by pairlist in pure python.
INFO pairlist_py: 125 ms
INFO pairlist_py: end.
24000 pairs
INFO pairlist_c: Neighboring pair list by pairlist in c.
INFO pairlist_c: end.
INFO pairlist_c: 16 ms
24000 pairs
```

```python
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
            lattice, maxdist=rc, cell=cell, engine=(
            pairs_py, pairs2_py)):
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

```

![benchmark](https://github.com/vitroid/PairList/raw/master/benchmark/benchmark.png)

## Algorithm

A simple cell division algorithm is implemented.

## Demo

It requires [GenIce](https://github.com/vitroid/GenIce) to make the test data.

```shell
% make test
% make test2
% make test3
```

## Requirements

* numpy


## Bugs

