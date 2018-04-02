#!/usr/bin/env python

import random
import cpairlist as cpl
import pairlist as pl
import numpy as np
import time

#1349221 3.978379011154175
#1349221 4.7704668045043945

z = np.array([[random.random() for i in range(3)] for j in range(10000)])
t0 = time.time()
a = set([frozenset((i,j)) for i,j in cpl.pairs(z,10,10,10)])
t1 = time.time()
print(len(a), t1-t0)
t0=t1
b = set([frozenset((i,j)) for i,j in pl.pairs_py(z,10,10,10)])
t1 = time.time()
print(len(b), t1-t0)
t0=t1
print(a-b)
print(b-a)

z2 = np.array([[random.random() for i in range(3)] for j in range(5000)])
a = set([(i,j) for i,j in cpl.pairs2(z,z2,10,10,10)])
print(len(a))
b = set([(i,j) for i,j in pl.pairs2_py(z,z2,10,10,10)])
print(len(b))
print(a-b)
print(b-a)


