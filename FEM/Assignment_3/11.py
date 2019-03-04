from sympy.solvers import solve
from sympy import Symbol
from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import *
from numpy.linalg import inv
from array import *
from scipy import linalg
x=Symbol('x')
xa=0.3333*2
xb=0.3333*3
psi_1 = (xb-x)/(xb-xa)
psi_2 = (x-xa)/(xb-xa)
# print(psi_1)
# print(psi_2)
f=-x*x
d=integrate(f*psi_1 ,(x, xa, xb))
e=integrate(f*psi_2 ,(x, xa, xb))
print('d',d)
print('e',e)
