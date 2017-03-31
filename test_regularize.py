#!/usr/bin/env python3


import numpy as np
import itertools as it

x = np.array([1.1, 0.1, -0.2])
y = np.array([0.3, 0.9, 0.2])
z = np.array([-0.2, -0.1, 1.2])

def regularize(x,y,z):
    error = np.dot(x,y)
    print(error)
    x_ort=x-(error/2)*y
    y_ort=y-(error/2)*x
    z_ort=np.cross(x_ort, y_ort)
    print(x_ort, y_ort, z_ort)
    x_new = 1./2*(3-np.dot(x_ort,x_ort))*x_ort
    y_new = 1./2*(3-np.dot(y_ort,y_ort))*y_ort
    z_new = 1./2*(3-np.dot(z_ort,z_ort))*z_ort
    return x_new, y_new, z_new

x,y,z = regularize(x,y,z)
x,y,z = regularize(x,y,z)
print(x,y,z)
for v in (x, y, z):
    print(np.dot(v,v))
for v1, v2 in it.combinations((x, y, z), 2):
    print(np.dot(v1,v2))
