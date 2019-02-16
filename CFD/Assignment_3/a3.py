from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from array import *
from scipy.sparse import *

T=[[]]
# one row for one time  step
#  (r,c) => (n,i)
# T.append([5,2,1])
# T[0].append(6)
r=3
for i in range(20):
    T[0].append(25)
print(T)
print(len(T))
for n in range(10):
    for i in range(0,20):
        T.append([0])
        T[n].append(1)
print(T)
