# PairList
Generates the pair list of atoms that are closer to each other than the
given threshold under the periodic boundary conditions.

version 0.2.10

## Usage

See `pairlist.h` for the function definition and `pairlist-test.c` for usage.

Python API is served in pairlist.py. Here is a sample code to use it.

```python
import pairlist as pl
import numpy as np

def FaceCenteredCubic(N):
    N2 = N*2
    ea = np.array([2**0.5,0,0])
    eb = np.array([0,2**0.5,0])
    ec = np.array([0,0,2**0.5])
    cell = np.array([ea,eb,ec])
    lattice = np.array([(x/N2, y/N2, z/N2)
                        for x in range(N2)
                        for y in range(N2)
                        for z in range(N2)
                        if (x+y+z)%2==0 ])
    return lattice, cell*N

lattice, cell = FaceCenteredCubic(2)
for i,j,d in pl.pairs_iter(lattice, rc=1.1, cell=cell):
    print(i,j,d)

```

```python
import pairlist as pl
import numpy as np

def SimpleCubic(N):
    cell = np.eye(3) * N
    lattice = np.array([(x/N, y/N, z/N)
                        for x in range(N)
                        for y in range(N)
                        for z in range(N)])
    return lattice, cell

N=2
lattice, cell = SimpleCubic(N)
# bipartile BCC lattice
lattice2 = lattice + 1/(N*2)

for i,j in pl.pairs_iter(lattice, rc=1.1*3**0.5/2, cell=cell, pos2=lattice2, distance=False):
    print(i,j)
```

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

