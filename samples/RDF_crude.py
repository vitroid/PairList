#!/usr/bin/env python

# crude RDF, without Numpy

import sys
from math import floor, pi
import gromacs


# round int


def rint(x):
    return floor(x + 0.5)


def main():
    atoms, cell = gromacs.load(sys.stdin)
    os = atoms['O']
    assert len(cell) == 3  # assume rect cell.
    density = len(os) / (cell[0] * cell[1] * cell[2])

    # histogram
    intv = 0.003
    maxbin = 300
    histo = [0.0 for i in range(maxbin)]

    for i in range(len(os)):
        xi = os[i]
        for j in range(i + 1, len(os)):
            xj = os[j]
            # relative position
            dx = xi[0] - xj[0]
            dy = xi[1] - xj[1]
            dz = xi[2] - xj[2]
            # treatment for the periodic boundary
            dx -= rint(dx / cell[0]) * cell[0]
            dy -= rint(dy / cell[1]) * cell[1]
            dz -= rint(dz / cell[2]) * cell[2]
            # distance
            r = (dx**2 + dy**2 + dz**2)**0.5
            # accumulate
            ir = int(r / intv)
            if ir < maxbin:
                histo[ir] += 1

    for ir, h in enumerate(histo):
        r = ir * intv
        # volume of the onion shell
        dv = ((r + intv)**3 - r**3) * 4. * pi / 3.
        # average number of particles in the shell
        dn = dv * density
        print(r, h / dn / len(os) * 2)


if __name__ == "__main__":
    main()
