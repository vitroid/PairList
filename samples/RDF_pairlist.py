#!/usr/bin/env python

import sys
from math import floor, pi
import numpy as np
import pairlist as pl
import gromacs


def main():
    atoms, cell = gromacs.load(sys.stdin)
    os = atoms["O"]
    assert len(cell) == 3  # assume rect cell.
    # pairlistはnumpy arrayのみ受け入れる。
    os = np.array(os)
    cell = np.array(cell)
    # セルの形状を表現する対角行列
    cellmat = np.diag(cell)
    # 相対位置に換算する(直方体なので簡単)
    rpos = os / cell
    density = os.shape[0] / np.prod(cell)

    # histogram
    intv = 0.003
    maxbin = 300  # int(min(min(cell[0],cell[1]),cell[2])/2 / intv)

    histo = [0.0 for i in range(maxbin)]

    # 対とその距離を計算し、ヒストグラムにする。
    # (Pure pythonの場合)
    for i, j, r in pl.pairs_iter(
        rpos, intv * maxbin, cellmat, distance=True, _engine=(pl.pairs_py, None)
    ):
        # accumulate
        ir = int(r / intv)
        if ir < maxbin:
            histo[ir] += 1

    for ir, h in enumerate(histo):
        r = ir * intv
        # volume of the onion shell
        dv = ((r + intv) ** 3 - r**3) * 4.0 * pi / 3.0
        # average number of particles in the shell
        dn = dv * density
        print(r, h / dn / len(os) * 2)


if __name__ == "__main__":
    main()
