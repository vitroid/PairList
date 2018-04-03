#!/usr/bin/env python

# crude RDF, without Numpy

import sys
from math import floor, pi
import numpy as np
import pairlist as pl
from test1_rdf import load_gro_o

def main():
    os, cell = load_gro_o(sys.stdin)
    assert len(cell) == 3 #assume rect cell.
    os = np.array(os)
    cell = np.array(cell)
    cellmat = np.diag(cell)
    # relative positions
    rpos = os / cell 
    density = os.shape[0] / np.product(cell)

    # histogram
    intv = 0.003
    maxbin = 300

    grid = pl.determine_grid(cellmat, intv*maxbin)
    i,j,r = pl.pairs_fine(rpos, intv*maxbin, cellmat, grid, distance=True, raw=True)
    # 対距離を整数化する
    delta = np.floor(r/intv)
    # 同じ値が何度出現したかをカウントする=ヒストグラムを作る。
    # この作業もnumpyにまかせてしまう。
    histo = dict(zip(*np.unique(delta, return_counts=True)))
    # 結果は距離をキーとし回数を値とする辞書。

    for ir in range(maxbin):
        h = 0.0
        r = float(ir)
        if r in histo:
            h = histo[r]
        # volume of the onion shell
        dv = ((r+intv)**3 - r**3)*4.*pi/3.
        # average number of particles in the shell
        dn = dv*density
        print(r,h/dn/len(os)*2)

if __name__ == "__main__":
    main() 
