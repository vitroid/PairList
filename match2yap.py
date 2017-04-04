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


def drawbox(origin, box):
    #print(box)
    center = (box[0] + box[1] + box[2])/2
    ori = origin - center
    for v in (np.zeros(3), box[1], box[2], box[1]+box[2]):
        print(yp.Line(v+ori, v+ori+box[0]), end="")
        #s += yp.Line(v+ori, v+ori+box[0])
    for v in (np.zeros(3), box[0], box[2], box[0]+box[2]):
        print(yp.Line(v+ori, v+ori+box[1]), end="")
        #s += yp.Line(v+ori, v+ori+box[1])
    for v in (np.zeros(3), box[0], box[1], box[0]+box[1]):
        print(yp.Line(v+ori, v+ori+box[2]), end="")
        #s += yp.Line(v+ori, v+origin+box[2])

def drawbox2(origin, box):
    center = (box[0] + box[1] + box[2])/2
    ori = origin - center
    print(yp.Color(2), end="")
    print(yp.Polygon([ori,ori+box[0],ori+box[1]]), end="")
    print(yp.Color(3), end="")
    print(yp.Polygon([ori,ori+box[1],ori+box[2]]), end="")
    print(yp.Color(4), end="")
    print(yp.Polygon([ori,ori+box[2],ori+box[0]]), end="")

def drawatoms(atoms, members=None):
    if members is None:
        for a in atoms:
            print(yp.Circle(a),end="")
    else:
        for i in members:
            print(yp.Circle(atoms[i]),end="")


cell, atoms = LoadGRO(open(sys.argv[1]))
unitcell, unitatoms = LoadAR3R(open(sys.argv[2]))
#in absolute coord
unitatoms = np.dot(unitatoms, unitcell)
#print(unitatoms)
dmin = 1e99
orig = None
for a in unitatoms:
    L = np.dot(a,a)
    if L < dmin:
        dmin = L
        orig = a
    
#s = ""
matched = set()
for line in sys.stdin:
    #parse the line
    cols = line.split()
    N = int(cols[0])
    members = [int(x) for x in cols[1:N+1]]
    matched |= set(members)
    msd     = float(cols[N+1])
    origin  = atoms[int(cols[N+2])].copy()  #atom at the matching center
    rotmat  = np.array([float(x) for x in cols[N+3:N+12]]).reshape((3,3))
    print(rotmat)
    #draw matched box
    boxorigin = np.dot(orig, rotmat)
    origin   -= boxorigin #corner of the cell
    box       = np.dot(unitcell, rotmat)
    uatoms    = np.dot(unitatoms, rotmat) + origin
    print(yp.Color(3),end="")
    print(yp.Layer(1),end="")
    drawbox(origin,box)

    print(yp.Layer(2),end="")
    drawbox2(origin,box)

    print(yp.Layer(3),end="")
    print(yp.Size(0.15),end="")
    drawatoms(uatoms)
    for i in range(len(uatoms)):
        print(yp.Line(uatoms[i], atoms[members[i]]),end="")
    #drawbox2(origin,box)
    #print(yp.Color(5), end="")
    #print(yp.Palygon([origin,origin+box[2],origin+box[0]]), end="")
#print(s)

print(yp.Size(0.1),end="")
print(yp.Color(0),end="")
print(yp.Layer(3),end="")
drawatoms(atoms, members=matched)
