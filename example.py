#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pairlist as pl

def test():
    file = open("pairlist-test3.gro")
    line = file.readline()
    natom = int(file.readline())
    Os = []
    Hs = []
    for i in range(natom):
        line = file.readline()
        cols = line.split()
        xyz = np.array([float(x) for x in cols[3:6]])
        if cols[1][0:2] == "OW":
            Os.append(xyz)
        else:
            Hs.append(xyz)
    Hs = np.array(Hs)
    Os = np.array(Os)
    #the last line is cell shape
    line = file.readline()
    cols = line.split()
    dimen = np.array([float(x) for x in cols])
    celli = 1.0 / dimen
    cell = np.diag(dimen)
    #Atoms must be in relative coordinate
    Hs[:] *= celli
    Os[:] *= celli
    assert len(Hs) == 864
    assert len(Os) == 432
    #find O-O pairs in 0.3 nm
    pairs = pl.pairlist(Os, 0.3, cell)
    assert len(pairs) == 864
    #find OH covalent bonds
    pairs = pl.pairlist(Os, 0.12, cell, Hs)
    assert len(pairs) == 864
    #find hydrogen bonds
    pairs = pl.pairlist(Os, 0.245, cell, Hs)
    pairs = pairs[pairs[:,2]> 0.12] #not covalent bond
    assert len(pairs) == 864
    for i,j,d in pairs:
        print("{0}--{1}\t{2}".format(int(i),int(j),d))
    print("{0} Number of HBs".format(len(pairs)))

test()
