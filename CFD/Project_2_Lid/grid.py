from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from array import *
from scipy.sparse import *

gama=2
x=[]
y=[]
nx=[]
ny=[]

for i in np.arange(0., 5., 0.2):
    nx.append(i)
for j in nx:
    x.append(1-(np.tanh(gama*(1-(2*j)/len(nx))))/(np.tanh(gama)))

for i in np.arange(0., 5., 0.2):
    ny.append(i)
for i in ny:
    y.append(1-(np.tanh(gama*(1-(2*i)/len(ny))))/(np.tanh(gama)))

for i in range(len(nx)):
    x.append(x[len(nx)+i-1]+x[len(nx)-(i+1)]-x[len(nx)-(i+2)])
    y.append(y[len(ny)+i-1]+y[len(ny)-(i+1)]-y[len(ny)-(i+2)])

# x.append(x[n-1]+x[n-1]-x[n-2])
# x.append(x[n]+x[n-2]-x[n-3])

for i in x:
    for k in y:
        plt.scatter(i,k,color='b',marker='+')
plt.show()
