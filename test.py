import numpy as np
import pairlist as pl
import random

z = np.array([[random.random() for x in range(3)] for y in range(10)])
print(z)
GX=3
GY=4
GZ=5
print(pl.pairs(z,GX,GY,GZ))
print(pl.pairs(z,2,3,4))
