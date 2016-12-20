#!/usr/bin/env python
# -*- coding: utf-8 -*-

####Note: xyz are in coord relative to the cell.


import math
import itertools
import numpy as np


def ArrangeAddress(xyz,grid):
    #residents in each grid cell
    residents = dict()
    for i in range(len(xyz)):
        mol = xyz[i]
        mol -= np.floor( mol )
        address = tuple((mol * grid).astype(int))
        if address not in residents:
            residents[address] = set()
        residents[address].add(i)
    return residents



def _pairlist(xyz,grid):
    #print "START Arrange"
    residents = ArrangeAddress(xyz,grid)
    #print "END Arrange"

    pair = set()
    #key-value pairs in the dictionary
    donecellpair = set()
    for address in residents:
        members = residents[address]
        ix,iy,iz = address
        #neighbor cells
        npa = np.array(address)
        k = np.array([(npa + j + grid)%grid for j in range(-1,2)])
        for a2 in itertools.product(k[:,0],k[:,1],k[:,2]):
            if address == a2:
                #print("ISOCELL",address,npa,a2)
                for a,b in itertools.combinations(members,2):
                    pair.add((a,b))
            else:
                if a2 in residents:
                    if not (frozenset((address,a2)) in donecellpair):
                        donecellpair.add(frozenset((address,a2)))
                        #print("HETEROCELL",address,a2,donecellpair)
                        for a in members:
                            for b in residents[a2]:
                                pair.add((a,b))
    #print "PAIRLIST finished"
    return pair


def _pairlist_hetero(xyz,xyz2,grid):
    #print "START Arrange"
    residents  = ArrangeAddress(xyz,grid)
    residents2 = ArrangeAddress(xyz2,grid)
    #print "END Arrange"

    pair = set()
    #key-value pairs in the dictionary
    donecellpair = set()
    for address in residents:
        members = residents[address]
        ix,iy,iz = address
        #neighbor cells
        npa = np.array(address)
        k = np.array([(npa + j + grid)%grid for j in range(-1,2)])
        for a2 in itertools.product(k[:,0],k[:,1],k[:,2]):
            if a2 in residents2:
                if not ((address,a2) in donecellpair):
                        donecellpair.add((address,a2))
                        #print("HETEROCELL",address,a2,donecellpair)
                        for a in members:
                            for b in residents2[a2]:
                                pair.add((a,b))
    #print "PAIRLIST finished"
    return pair


#assume xyz and box are numpy.array
def _pairlist_fine(xyz,rc,cell,grid,distance=True):
    newpairs = []
    for i,j in _pairlist(xyz,grid):
        moli = xyz[i]
        molj = xyz[j]
        d = moli-molj
        d -= np.floor( d + 0.5 )
        d = np.dot(d,cell)
        rr = np.dot(d,d)
            
        if rr < rc**2:
            if distance:
                newpairs.append((i,j,math.sqrt(rr)))
            else:
                newpairs.append((i,j))
    return np.array(newpairs)


#assume xyz and box are numpy.array
def _pairlist_fine_hetero(xyz,xyz2,rc,cell,grid,distance=True):
    newpairs = []
    for i,j in _pairlist_hetero(xyz,xyz2,grid):
        moli = xyz[i]
        molj = xyz2[j]
        d = moli-molj
        d -= np.floor( d + 0.5 )
        d = np.dot(d,cell)
        rr = np.dot(d,d)
            
        if rr < rc**2:
            if distance:
                newpairs.append((i,j,math.sqrt(rr)))
            else:
                newpairs.append((i,j))
    return np.array(newpairs)


#
def determine_grid(cell, radius):
    ct = cell.transpose()
    a = ct[0]
    b = ct[1]
    c = ct[2]
    al = np.linalg.norm(a)   #vector length
    bl = np.linalg.norm(b)
    cl = np.linalg.norm(c)
    ae = a / al              #unit vectors
    be = b / bl
    ce = c / cl
    ad = np.dot(ae,np.cross(be,ce)) #distance to the bc plane
    bd = np.dot(be,np.cross(ce,ae))
    cd = np.dot(ce,np.cross(ae,be))
    ax = radius / ad        # required length of a vector to contain a sphere of radius 
    bx = radius / bd
    cx = radius / cd
    gf = np.array([al/ax, bl/bx, cl/cx])  # required number of grid cells
    #print(cell,radius,gf)
    #import sys
    #sys.exit(1)
    return np.floor(gf).astype(int)


def pairlist(xyz,rc,cell,xyz2=None,distance=True):
    grid = determine_grid(cell, rc)
    if xyz2 is None:
        return _pairlist_fine(xyz,rc,cell,grid,distance)
    else:
        return _pairlist_fine_hetero(xyz,xyz2,rc,cell,grid,distance)

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
    pairs = pairlist(Os, 0.3, cell)
    assert len(pairs) == 864
    #find OH covalent bonds
    pairs = pairlist(Os, 0.12, cell, Hs)
    assert len(pairs) == 864
    #find hydrogen bonds
    pairs = pairlist(Os, 0.245, cell, Hs)
    pairs = pairs[pairs[:,2]> 0.12] #not covalent bond
    assert len(pairs) == 864

if __name__ == "__main__":
    test()
