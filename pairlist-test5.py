#!/usr/bin/env python

import random
import cpairlist as cpl
import pairlist as pl
import numpy as np
import time


#1352013 0.007421970367431641
#1352013 0.5289030075073242
#209551 16.11849093437195
z = np.array([[random.random() for i in range(3)] for j in range(10000)])
t0 = time.time()
a = cpl.pairs(z,10,10,10)
t1 = time.time()
print(len(a), t1-t0)
t0=t1
b = list(pl.pairs_py(z,10,10,10))
t1 = time.time()
print(len(b), t1-t0)

cell = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
grid = np.array([10,10,10])
rc=0.1
t0 = time.time()
c = pl.pairs_fine(z,rc,cell,grid,distance=True)
t1 = time.time()
print(list(c), t1-t0)
