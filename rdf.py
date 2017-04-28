#!/usr/bin/env python3

import sys
from math import pi
import numpy as np
from pairlist import pairlist
from collections import defaultdict


def test():
    center = None
    if len(sys.argv) > 1:
        center = sys.argv[1]
    file = sys.stdin
    line = file.readline()
    natom = int(file.readline())
    atoms = defaultdict(list)
    for i in range(natom):
        line = file.readline()
        cols = line[20:].split()
        xyz = np.array([float(x) for x in cols[:3]])
        atom = line[10:15].split()
        #print(atom)
        atoms[atom[0]].append(xyz)
    #the last line is cell shape
    line = file.readline()
    cols = line.split()
    dimen = np.array([float(x) for x in cols])
    celli = 1.0 / dimen
    cell = np.diag(dimen)
    atoms_ordered = sorted(atoms.keys())
    
    #relative
    if center is None:
        #regard them as identical
        allatoms = []
        for atom in atoms_ordered:
            allatoms += atoms[atom]
        atoms = dict()
        atoms["Any"] = np.array(allatoms) * celli
        center = "Any"
    else:
        for atom in atoms_ordered:
            atoms[atom] = np.array(atoms[atom]) * celli
    rdf = dict()
    for atom in atoms_ordered:
        if center == atom:
            pairs = pairlist(atoms[center], 10.0, cell)
        else:
            pairs = pairlist(atoms[center], 10.0, cell, xyz2=atoms[atom])
        rdf[atom] = [0.0] * 10001
        thick = 0.001 #nm
        for i,j,d in pairs:
            bin = int(d/thick+0.5)
            rdf[atom][bin] += 1
    dens = dict()
    for atom in atoms_ordered:
        dens[atom] = len(atoms[atom])/np.product(dimen)
    values = ["# r"] + atoms_ordered
    print(*values)
    for bin,dist in enumerate(rdf[center]):
        r = bin * thick  #in nm
        if r > 0:
            values = []
            for atom in atoms_ordered:
                skin = 4*pi*r**2 * thick*dens[atom]
                values.append(rdf[atom][bin]/skin/(len(atoms[center])-1))
            print(r, *values)



test()

        
