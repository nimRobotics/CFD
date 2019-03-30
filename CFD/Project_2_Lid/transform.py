from sympy.solvers import solve
from sympy import Symbol
from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from scipy.sparse import *
from numpy.linalg import inv
from array import *
from scipy import linalg
import scipy.linalg


x = np.arange(0,2.001,0.1)
# y=(x-1)**3+1
beta =8
y = np.arctan(x-1)
print(y)
t = zeros(1,len(x))
print(t)
plt.scatter(y,t)
# plt.subplot(121)
# plt.scatter(y,t)
# plt.subplot(122)
# plt.scatter(yd,t)
plt.show()
