#!/usr/bin/env python
# -*- coding: utf-8 -*-
# even: stable; odd: develop
__version__ = "0.2.6"

import math
import itertools as it
import numpy as np
import logging
from cpairlist import pairs, pairs2

def Address(pos,grid):
    #residents in each grid cell
    mol = pos % 1 % 1  # avoid cancellation
    return tuple((mol * grid).astype(int))


def ArrangeAddress(xyz,grid):
    #residents in each grid cell
    residents = dict()
    for i in range(len(xyz)):
        address = Address(xyz[i], grid)
        if address not in residents:
            residents[address] = set()
        residents[address].add(i)
    return residents



def pairs_py(xyz,GX,GY,GZ):
    grid = np.array([GX,GY,GZ])
    logger = logging.getLogger()
    logger.debug("START Arrange")
    residents = ArrangeAddress(xyz,grid)
    logger.debug("END Arrange")

    #key-value pairs in the dictionary
    donecellpair = set()
    pairs = []
    for address in residents:
        members = residents[address]
        #neighbor cells
        npa = np.array(address)
        k = np.array([(npa + j + grid)%grid for j in range(-1,2)])
        for a2 in it.product(k[:,0],k[:,1],k[:,2]):
            if address == a2:
                if frozenset((address,a2)) not in donecellpair:
                    donecellpair.add(frozenset((address,a2)))
                    for a,b in it.combinations(members,2):
                        pairs.append((a,b))
            else:
                if a2 in residents:
                    if frozenset((address,a2)) not in donecellpair:
                        donecellpair.add(frozenset((address,a2)))
                        for a in members:
                            for b in residents[a2]:
                                pairs.append((a,b))
    return np.array(pairs)

def pairs2_py(xyz,xyz2,GX,GY,GZ):
    grid = np.array([GX,GY,GZ])
    return _pairs_hetero(xyz,xyz2,grid)

def _pairs_hetero(xyz,xyz2,grid):
    logger = logging.getLogger()
    logger.debug("START Arrange")
    residents  = ArrangeAddress(xyz,grid)
    residents2 = ArrangeAddress(xyz2,grid)
    logger.debug("END Arrange")

    #key-value pairs in the dictionary
    donecellpair = set()
    pairs = []
    for address in residents:
        members = residents[address]
        ix,iy,iz = address
        #neighbor cells
        npa = np.array(address)
        k = np.array([(npa + j + grid)%grid for j in range(-1,2)])
        for a2 in it.product(k[:,0],k[:,1],k[:,2]):
            if a2 in residents2:
                if not ((address,a2) in donecellpair):
                    donecellpair.add((address,a2))
                    for a in members:
                        for b in residents2[a2]:
                            pairs.append((a,b))
    return np.array(pairs)


                                
#assume xyz and box are numpy.array
def pairs_fine_slow(xyz,rc,cell,grid,distance=True):
    logger= logging.getLogger()
    for i,j in pairs(xyz,*grid):
        moli = xyz[i]
        molj = xyz[j]
        d = moli-molj
        d -= np.floor( d + 0.5 )
        d = np.dot(d,cell)
        rr = np.dot(d,d)
        if rr < rc**2:
            if distance:
                yield i,j,math.sqrt(rr)
            else:
                yield i,j


# fully numpy style
def pairs_fine(xyz,rc,cell,grid,distance=True, raw=False, pairs_engine=pairs):
    logger= logging.getLogger()
    #p = np.array(list(pairs_py(xyz, *grid)))
    p = pairs_engine(xyz, *grid)
    idx0 = p[:,0]
    #for i in range(idx0.shape[0]):
    #    if idx0[i] > 1000:
    #        print(i,idx0[i])
    idx1 = p[:,1]
    p0 = xyz[idx0]
    p1 = xyz[idx1]
    d  = p0 - p1
    d -= np.floor(d+0.5)
    a  = np.dot(d,cell)
    L  = np.linalg.norm(a, axis=1)
    cond = L<rc
    # pickup elements satisfying the condition.
    j0 = np.compress(cond, idx0)
    j1 = np.compress(cond, idx1)
    if raw:
        if not distance:
            return j0, j1
        else:
            Ls = np.compress(cond, L)
            return j0,j1,Ls         #no zipping
    else:
        if not distance:
            return np.column_stack((j0, j1))
        else:
            Ls = np.compress(cond, L)
            # return np.column_stack(j0, j1) all the values becomes float...
            return zip(j0,j1,Ls)    #list of tuples
                

def pairs_crude(xyz,rc,cell,distance=True):
    logger = logging.getLogger()
    # logger.debug(xyz)
    logger.debug(rc)
    logger.debug(cell)
    logger.debug(distance)
    for i,j in it.combinations(range(len(xyz)),2):
        moli = xyz[i]
        molj = xyz[j]
        d = moli-molj
        d -= np.floor( d + 0.5 )
        r = np.dot(d,cell)
        rr = np.dot(r,r)
            
        if rr < rc**2:
            # logger.debug((d,r,rr,rc**2))
            if distance:
                yield i,j,math.sqrt(rr)
            else:
                yield i,j




#assume xyz and box are numpy.array
def pairs_fine_hetero_slow(xyz,xyz2,rc,cell,grid,distance=True):
    for i,j in pairs2(xyz,xyz2,*grid):
        moli = xyz[i]
        molj = xyz2[j]
        d = moli-molj
        d -= np.floor( d + 0.5 )
        d = np.dot(d,cell)
        rr = np.dot(d,d)
            
        if rr < rc**2:
            if distance:
                yield i,j,math.sqrt(rr)
            else:
                yield i,j




def pairs_fine_hetero(xyz,xyz2,rc,cell,grid,distance=True, raw=False, pairs_engine=pairs2):
    logger= logging.getLogger()
    p = pairs_engine(xyz, xyz2, *grid)
    idx0 = p[:,0]
    idx1 = p[:,1]
    p0 = xyz[idx0]
    p1 = xyz2[idx1]
    d  = p0 - p1
    d -= np.floor(d+0.5)
    a  = np.dot(d,cell)
    L  = np.linalg.norm(a, axis=1)
    cond = L<rc
    # pickup elements satisfying the condition.
    j0 = np.compress(cond, idx0)
    j1 = np.compress(cond, idx1)
    if raw:
        if not distance:
            return j0, j1
        else:
            Ls = np.compress(cond, L)
            return j0,j1,Ls         #no zipping
    else:
        if not distance:
            return np.column_stack((j0, j1))
        else:
            Ls = np.compress(cond, L)
            # return np.column_stack(j0, j1) all the values becomes float...
            return zip(j0,j1,Ls)    #list of tuples
                
#
def determine_grid(cell, radius):
    logger = logging.getLogger()
    ct = cell  #.transpose()
    #Cell vectors
    a = ct[0]
    b = ct[1]
    c = ct[2]
    logger.debug("cell a {0}".format(a))
    logger.debug("cell b {0}".format(b))
    logger.debug("cell c {0}".format(c))
    #Edge lengths
    al = np.linalg.norm(a)   #vector length
    bl = np.linalg.norm(b)
    cl = np.linalg.norm(c)
    #Unit vectors of the axes.
    ae = a / al              #unit vectors
    be = b / bl
    ce = c / cl
    #Distance between the parallel faces
    an = np.dot(a,np.cross(be,ce)) #distance to the bc plane
    bn = np.dot(b,np.cross(ce,ae))
    cn = np.dot(c,np.cross(ae,be))
    gf = np.array([an/radius, bn/radius, cn/radius])  # required number of grid cells
    if gf[0] < 1.0:
        gf[0] = 1.0
    if gf[1] < 1.0:
        gf[1] = 1.0
    if gf[2] < 1.0:
        gf[2] = 1.0
    #Check the lengths of four diagonals.
    logger.debug("Grid divisions: {0}".format(np.floor(gf)))
    #print(cell,radius,gf)
    #import sys
    #sys.exit(1)
    return np.floor(gf).astype(int)



def main():
    xyz = []
    for x in range(2):
        for y in range(2):
            for z in range(2):
                xyz.append(np.array((x/100.+1,y/100.+1,z/100.+1)) / 4.)
    xyz2 = []
    for x in range(2):
        for y in range(2):
            for z in range(2):
                xyz2.append(np.array((x/100.+2,y/100.+1,z/100.+1)) / 4.)
    xyz = np.array(xyz)
    xyz2 = np.array(xyz2)
    box = np.diag((4,4,4))
    rc = 1.000000001
    grid = determine_grid(box,rc)
    for i,j,l in pairs_fine_hetero(xyz,xyz2,rc,box,grid):
        print(i,j,l)

if __name__ == "__main__":
    main()
