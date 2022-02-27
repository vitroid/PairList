#!/usr/bin/env python

# Find the hydrogen bonds.

import sys
import numpy as np
import pairlist as pl
import gromacs
import networkx as nx

def main():
    atoms, cell = gromacs.load(sys.stdin)
    Os = np.array(atoms['O'])
    Hs = np.array(atoms['H'])
    assert len(cell) == 3  # assume rect cell.
    cell = np.array(cell)
    # cell shape matrix
    cellmat = np.diag(cell)

    g = nx.DiGraph()

    # hydrogen-oxygen pairs.
    for i, j, r in pl.pairs_iter(Hs, 0.245, cellmat, pos2=Os, distance=True, fractional=False):
        # if the distance is shorter than 0.11 nm, it is a covalent bond.
        # if the distance is longer than 0.245 nm, they are not bonded.
        if 0.12 < r: 
            # HB points from a hydrogen atom to an oxygen atom
            g.add_edge(i,j)
    print(f"{len(Os)} oxygen atoms")
    print(f"{len(Hs)} hyrogen atoms")
    print(f"{len(g.edges())} hydrogen bonds")


if __name__ == "__main__":
    main()
