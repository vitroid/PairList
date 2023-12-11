#!/usr/bin/env python
# -*- coding: utf-8 -*-
# even: stable; odd: develop
"""
.. include:: ./README.md
"""

import math
import itertools as it
import numpy as np
from logging import getLogger
from cpairlist import pairs, pairs2
from typing import Generator

__all__ = ["pairs_iter"]


def Address(pos, grid):
    # residents in each grid cell
    mol = pos % 1 % 1  # avoid cancellation
    return tuple((mol * grid).astype(int))


def ArrangeAddress(xyz, grid):
    # residents in each grid cell
    residents = dict()
    for i in range(len(xyz)):
        address = Address(xyz[i], grid)
        if address not in residents:
            residents[address] = set()
        residents[address].add(i)
    return residents


def pairs_py(xyz, GX, GY, GZ):
    grid = np.array([GX, GY, GZ])
    logger = getLogger()
    logger.debug("START Arrange")
    residents = ArrangeAddress(xyz, grid)
    logger.debug("END Arrange")

    # key-value pairs in the dictionary
    donecellpair = set()
    pairs = []
    for address in residents:
        members = residents[address]
        # neighbor cells
        npa = np.array(address)
        k = np.array([(npa + j + grid) % grid for j in range(-1, 2)])
        for a2 in it.product(k[:, 0], k[:, 1], k[:, 2]):
            if address == a2:
                if frozenset((address, a2)) not in donecellpair:
                    donecellpair.add(frozenset((address, a2)))
                    for a, b in it.combinations(members, 2):
                        pairs.append((a, b))
            else:
                if a2 in residents:
                    if frozenset((address, a2)) not in donecellpair:
                        donecellpair.add(frozenset((address, a2)))
                        for a in members:
                            for b in residents[a2]:
                                pairs.append((a, b))
    return np.array(pairs)


def pairs2_py(xyz, xyz2, GX, GY, GZ):
    grid = np.array([GX, GY, GZ])
    return _pairs_hetero(xyz, xyz2, grid)


def _pairs_hetero(xyz, xyz2, grid):
    logger = getLogger()
    logger.debug("START Arrange")
    residents = ArrangeAddress(xyz, grid)
    residents2 = ArrangeAddress(xyz2, grid)
    logger.debug("END Arrange")

    # key-value pairs in the dictionary
    donecellpair = set()
    pairs = []
    for address in residents:
        members = residents[address]
        ix, iy, iz = address
        # neighbor cells
        npa = np.array(address)
        k = np.array([(npa + j + grid) % grid for j in range(-1, 2)])
        for a2 in it.product(k[:, 0], k[:, 1], k[:, 2]):
            if a2 in residents2:
                if not ((address, a2) in donecellpair):
                    donecellpair.add((address, a2))
                    for a in members:
                        for b in residents2[a2]:
                            pairs.append((a, b))
    return np.array(pairs)


# assume xyz and box are numpy.array
def pairs_fine_slow(xyz, rc, cell, grid=None, distance=True):
    logger = getLogger()
    if grid is None:
        grid = determine_grid(cell, rc)
    for i, j in pairs(xyz, *grid):
        moli = xyz[i]
        molj = xyz[j]
        d = moli - molj
        d -= np.floor(d + 0.5)
        d = np.dot(d, cell)
        rr = np.dot(d, d)
        if rr < rc**2:
            if distance:
                yield i, j, math.sqrt(rr)
            else:
                yield i, j


def pairs_nopbc_iter(pos, maxdist=None, pos2=None, distance=False):
    """
    Iterator to find pairs in an open space.
    """
    # assert pos2 is None, "Pairs for hetero group without PBC is not implemented yet."

    def assign_to_bins(atoms):
        # occupants[x] is a list of atom labels in a bin x.
        occupants = dict()
        for i, pos in enumerate(atoms):
            # from position to bin
            ipos = tuple(np.floor(pos / maxdist).astype(int))
            # insert the atom in the bin
            if ipos not in occupants:
                occupants[ipos] = []
            occupants[ipos].append(i)
        return occupants

    def nearby(A, hetero=False):
        """
        Obtain the list of atoms near A.

        A: Label of the target atom
        """
        # bin of A
        p, q, r = np.floor(pos[A] / maxdist).astype(int)

        # adjacent bins
        for ix in range(p - 1, p + 2):
            for iy in range(q - 1, q + 2):
                for iz in range(r - 1, r + 2):
                    if (ix, iy, iz) in occupants:
                        # candidates for the neighbors
                        for B in occupants[ix, iy, iz]:
                            # avoid double counts.
                            if hetero:
                                # relative vector
                                d = pos[A] - pos2[B]
                                d2 = d @ d
                                # if the distance is shorter than the threshold,
                                if d2 < maxdist**2:
                                    yield B, d2**0.5
                            elif A < B:
                                # relative vector
                                d = pos[A] - pos[B]
                                d2 = d @ d
                                # if the distance is shorter than the threshold,
                                if d2 < maxdist**2:
                                    yield B, d2**0.5

    if pos2 is None:
        # homo pairs
        occupants = assign_to_bins(pos)
        for i in range(len(pos)):
            for j, L in nearby(i):
                if distance:
                    yield i, j, L
                else:
                    yield i, j

    else:
        # hetero pairs
        occupants = assign_to_bins(pos2)
        for i in range(len(pos)):
            for j, L in nearby(i, hetero=True):
                if distance:
                    yield i, j, L
                else:
                    yield i, j


# wrapper
def pairs_iter(
    pos: np.ndarray,
    maxdist: float = None,
    cell: np.ndarray[3, 3] = None,
    fractional: bool = True,
    pos2: np.ndarray = None,
    rc: float = None,
    distance: bool = True,
    _raw: bool = False,
    _engine=(pairs, pairs2),
) -> Generator[tuple, None, None]:
    """Yields pairs within the specified distance.

    Args:
        pos (np.ndarray): Positions of the atoms. (in fractional or absolute coordinate)
        maxdist (formerly rc) (float, optional): Upper limit of the pair distance (in the same unit as the cell). Defaults to None.
        cell (np.ndarray[3, 3], optional): Shape of the cell in a 3x3 numpy array. [[ax,ay,az],[bx,by,bz],[cx,cy,cz]]. Defaults to None.
        fractional (bool, optional): If true, pos must be given in fractional coordinate. Defaults to True.
        pos2 (np.ndarray, optional): Specify when you want to find the pairs between two sets of positions. Defaults to None.
        distance (bool, optional): If true, distance between two positions are also returned. Defaults to True.

    Yields:
        Generator[tuple, None, None]: indices of a pair | indices and the distance between the pair (distance=True)
    """
    logger = getLogger()
    if rc is not None:
        logger.warning("rc is deprecated. Use maxdist instead.")
        assert maxdist is None, "rc and maxdist are specified at a time."
        maxdist = rc
    if cell is None:
        return pairs_nopbc_iter(pos, maxdist=maxdist, pos2=pos2, distance=distance)
    grid = None
    if fractional:
        rpos = pos
        rpos2 = pos2
    else:
        celli = np.linalg.inv(cell)
        rpos = pos @ celli
        if pos2 is None:
            rpos2 = None
        else:
            rpos2 = pos2 @ celli

    if rpos2 is None:
        return pairs_fine(
            rpos,
            maxdist,
            cell,
            grid=grid,
            distance=distance,
            raw=_raw,
            engine=_engine[0],
        )
    else:
        return pairs_fine_hetero(
            rpos,
            rpos2,
            maxdist,
            cell,
            grid=grid,
            distance=distance,
            raw=_raw,
            engine=_engine[1],
        )


# fully numpy style
def pairs_fine(xyz, rc, cell, grid=None, distance=True, raw=False, engine=pairs):
    logger = getLogger()
    if grid is None:
        grid = determine_grid(cell, rc)
    p = engine(xyz, *grid)
    idx0 = p[:, 0]
    # for i in range(idx0.shape[0]):
    #    if idx0[i] > 1000:
    #        print(i,idx0[i])
    idx1 = p[:, 1]
    p0 = xyz[idx0]
    p1 = xyz[idx1]
    d = p0 - p1
    d -= np.floor(d + 0.5)
    a = np.dot(d, cell)
    L = np.linalg.norm(a, axis=1)
    cond = L < rc
    # pickup elements satisfying the condition.
    j0 = np.compress(cond, idx0)
    j1 = np.compress(cond, idx1)
    if raw:
        if not distance:
            return j0, j1
        else:
            Ls = np.compress(cond, L)
            return j0, j1, Ls  # no zipping
    else:
        if not distance:
            return np.column_stack((j0, j1))
        else:
            Ls = np.compress(cond, L)
            # return np.column_stack(j0, j1) all the values becomes float...
            return zip(j0, j1, Ls)  # list of tuples


def pairs_crude(xyz, rc, cell, distance=True):
    logger = getLogger()
    # logger.debug(xyz)
    logger.debug(rc)
    logger.debug(cell)
    logger.debug(distance)
    for i, j in it.combinations(range(len(xyz)), 2):
        moli = xyz[i]
        molj = xyz[j]
        d = moli - molj
        d -= np.floor(d + 0.5)
        r = np.dot(d, cell)
        rr = np.dot(r, r)

        if rr < rc**2:
            # logger.debug((d,r,rr,rc**2))
            if distance:
                yield i, j, math.sqrt(rr)
            else:
                yield i, j


# assume xyz and box are numpy.array
def pairs_fine_hetero_slow(xyz, xyz2, rc, cell, grid=None, distance=True):
    if grid is None:
        grid = determine_grid(cell, rc)
    for i, j in pairs2(xyz, xyz2, *grid):
        moli = xyz[i]
        molj = xyz2[j]
        d = moli - molj
        d -= np.floor(d + 0.5)
        d = np.dot(d, cell)
        rr = np.dot(d, d)

        if rr < rc**2:
            if distance:
                yield i, j, math.sqrt(rr)
            else:
                yield i, j


def pairs_fine_hetero(
    xyz, xyz2, rc, cell, grid=None, distance=True, raw=False, engine=pairs2
):
    logger = getLogger()
    if grid is None:
        grid = determine_grid(cell, rc)
    p = engine(xyz, xyz2, *grid)
    idx0 = p[:, 0]
    idx1 = p[:, 1]
    p0 = xyz[idx0]
    p1 = xyz2[idx1]
    d = p0 - p1
    d -= np.floor(d + 0.5)
    a = np.dot(d, cell)
    L = np.linalg.norm(a, axis=1)
    cond = L < rc
    # pickup elements satisfying the condition.
    j0 = np.compress(cond, idx0)
    j1 = np.compress(cond, idx1)
    if raw:
        if not distance:
            return j0, j1
        else:
            Ls = np.compress(cond, L)
            return j0, j1, Ls  # no zipping
    else:
        if not distance:
            return np.column_stack((j0, j1))
        else:
            Ls = np.compress(cond, L)
            # return np.column_stack(j0, j1) all the values becomes float...
            return zip(j0, j1, Ls)  # list of tuples


# @lru_cache(maxsize=None)
def determine_grid(cell, radius):
    """
    Determine grid division based on the cutoff radius.
    """
    logger = getLogger()
    ct = cell  # .transpose()
    # Cell vectors
    a = ct[0]
    b = ct[1]
    c = ct[2]
    logger.debug("cell a {0}".format(a))
    logger.debug("cell b {0}".format(b))
    logger.debug("cell c {0}".format(c))
    # Edge lengths
    al = np.linalg.norm(a)  # vector length
    bl = np.linalg.norm(b)
    cl = np.linalg.norm(c)
    # Unit vectors of the axes.
    ae = a / al  # unit vectors
    be = b / bl
    ce = c / cl
    # Distance between the parallel faces
    an = np.dot(a, np.cross(be, ce))  # distance to the bc plane
    bn = np.dot(b, np.cross(ce, ae))
    cn = np.dot(c, np.cross(ae, be))
    # required number of grid cells
    gf = np.array([an / radius, bn / radius, cn / radius])
    if gf[0] < 1.0:
        gf[0] = 1.0
    if gf[1] < 1.0:
        gf[1] = 1.0
    if gf[2] < 1.0:
        gf[2] = 1.0
    # Check the lengths of four diagonals.
    logger.debug("Grid divisions: {0}".format(np.floor(gf)))
    # print(cell,radius,gf)
    # import sys
    # sys.exit(1)
    return np.floor(gf).astype(int)


def main():
    xyz = []
    for x in range(2):
        for y in range(2):
            for z in range(2):
                xyz.append(
                    np.array((x / 100.0 + 1, y / 100.0 + 1, z / 100.0 + 1)) / 4.0
                )
    xyz2 = []
    for x in range(2):
        for y in range(2):
            for z in range(2):
                xyz2.append(
                    np.array((x / 100.0 + 2, y / 100.0 + 1, z / 100.0 + 1)) / 4.0
                )
    xyz = np.array(xyz)
    xyz2 = np.array(xyz2)
    box = np.diag((4, 4, 4))
    rc = 1.000000001
    for i, j, l in pairs_fine_hetero(xyz, xyz2, rc, box):
        print(i, j, l)


if __name__ == "__main__":
    main()
