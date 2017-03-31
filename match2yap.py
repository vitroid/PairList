#!/usr/bin/env python3

import numpy as np
import sys
from genice import yaplotlib as yp

def LoadGRO(file):
    line = file.readline()
    natom = int(file.readline())
    Os = []
    for i in range(natom):
        line = file.readline()
        cols = line[20:].split()[:3]
        xyz = [float(x) for x in cols]
        if " O" in line:
            Os.append(xyz)
    Os = np.array(Os)*10 #in angstrom
    #the last line is cell shape
    line  = file.readline()
    cols  = line.split()
    dimen = np.array([float(x) for x in cols])
    cell  = np.diag(dimen)*10 #in angstrom
    return cell, Os


def LoadAR3R(file):
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        if "@BOX3" == line[:5]:
            line  = file.readline()
            cols  = line.split()
            dimen = np.array([float(x) for x in cols])
            cell  = np.diag(dimen)
        elif "@AR3R" == line[:5]:
            natom = int(file.readline())
            Os = []
            for i in range(natom):
                line = file.readline()
                cols = line.split()
                xyz = np.array([float(x) for x in cols])
                xyz -= np.floor( xyz + 0.5 )
                Os.append(xyz)
            Os = np.array(Os)
    #the last line is cell shape
    return cell, Os


cell, atoms = LoadGRO(open(sys.argv[1]))
unitcell, unitatoms = LoadAR3R(open(sys.argv[2]))

dmin = 1e99
orig = None
for a in unitatoms:
    L = np.dot(a,a)
    if L < dmin:
        dmin = L
        orig = a
    
    
s = ""
for line in sys.stdin:
    cols = line.split()
    N = int(cols[0])
    members = [int(x) for x in cols[1:N+1]]
    msd     = float(cols[N+1])
    origin  = atoms[int(cols[N+2])]  #atom at the matching center
    rotmat  = np.array([float(x) for x in cols[N+3:N+12]]).reshape((3,3))
    #draw matched box
    boxorigin = np.dot(orig, rotmat)
    origin   -= boxorigin #corner of the cell
    box       = np.dot(unitcell, rotmat)
    print(box)
    for v in (np.zeros(3), box[1], box[2], box[1]+box[2]):
        s += yp.Line(v+origin, v+origin+box[0])
    for v in (np.zeros(3), box[0], box[2], box[0]+box[2]):
        s += yp.Line(v+origin, v+origin+box[1])
    for v in (np.zeros(3), box[0], box[1], box[0]+box[1]):
        s += yp.Line(v+origin, v+origin+box[2])
print(s)
