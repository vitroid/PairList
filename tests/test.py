import time
import numpy as np
import cpairlist as cpl
import pairlist as pl
import random

z = np.array([[random.random() for x in range(3)] for y in range(10000)])
z2 = np.array([[random.random() for x in range(3)] for y in range(10000)])
GX = 10
GY = 10
GZ = 10
t0 = time.time()
print(len(cpl.pairs(z, GX, GY, GZ)))
t1 = time.time()
print(t1 - t0)
print(len(list(pl.pairs(z, GX, GY, GZ))))
t2 = time.time()
print(t2 - t1)
print(len(cpl.pairs2(z, z2, GX, GY, GZ)))
t3 = time.time()
print(t3 - t2)
print(len(list(pl.pairs2(z, z2, GX, GY, GZ))))
t4 = time.time()
print(t4 - t3)
