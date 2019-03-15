"""
MIT License

Copyright (c) 2019 Aakash Yadav

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Reference(s):
http://caefn.com/cfd/hyperbolic-tangent-stretching-grid
"""

from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from array import *
from scipy.sparse import *

def grid(nx,ny, gama):
    # use ty tx to change number of elements
    tx=(2)/((nx+1)+1)
    ty=(2)/((ny+1)+1)
    x=[]
    y=[]
    nx=[]
    ny=[]

    # x elements on left half
    for i in np.arange(0., 1., tx):
        nx.append(i)
    for j in nx:
        x.append(1-(np.tanh(gama*(1-(2*j)/len(nx))))/(np.tanh(gama)))

    # y elements on right half
    for i in np.arange(0., 1., ty):
        ny.append(i)
    for i in ny:
        y.append(1-(np.tanh(gama*(1-(2*i)/len(ny))))/(np.tanh(gama)))

    # mirroring x and y elements for the right half
    for i in range(len(nx)-1):
        x.append(x[len(nx)+i-1]+x[len(nx)-(i+1)]-x[len(nx)-(i+2)])
    for i in range(len(ny)-1):
        y.append(y[len(ny)+i-1]+y[len(ny)-(i+1)]-y[len(ny)-(i+2)])

    xd=[]
    yh=[]
    for i in x:
        xd.append(i/(x[len(x)-1]))
    for i in y:
        yh.append(i/(y[len(y)-1]))
    return(xd,yh)

# number of grid points in x direction
nx=128
# number of grid points un y direction
ny=128
# use gama to change the gaussian distibution
# as gama --> 0 grid becomes uniform
gama=15
c,d=grid(nx,ny,gama)
print(c)
print("\n",d)

for i in c:
    for k in d:
        plt.scatter(i,k,color='b',marker='+')
plt.show()
