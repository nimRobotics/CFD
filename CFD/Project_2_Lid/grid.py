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

def grid(H,D):
    # use gama to change the gaussian distibution
    gama=3
    # use ty tx to change number of elements
    tx=0.1
    ty=0.1
    x=[]
    y=[]
    nx=[]
    ny=[]

    for i in np.arange(0., 3., tx):
        nx.append(i)
    for j in nx:
        x.append(1-(np.tanh(gama*(1-(2*j)/len(nx))))/(np.tanh(gama)))

    for i in np.arange(0., 3., ty):
        ny.append(i)
    for i in ny:
        y.append(1-(np.tanh(gama*(1-(2*i)/len(ny))))/(np.tanh(gama)))

    for i in range(len(nx)-1):
        x.append(x[len(nx)+i-1]+x[len(nx)-(i+1)]-x[len(nx)-(i+2)])
    for i in range(len(ny)-1):
        y.append(y[len(ny)+i-1]+y[len(ny)-(i+1)]-y[len(ny)-(i+2)])
    # D=128
    # H=128
    xd=[]
    yh=[]
    for i in x:
        xd.append((i/(x[len(x)-1]))*D)
    for i in y:
        yh.append((i/(y[len(y)-1]))*H)
    return(xd,yh)
    # NOTE: use below code to print grid
    # for i in x:
    #     for k in y:
    #         plt.scatter(i,k,color='b',marker='+')
    # plt.show()

h=128
d=128
c,d=grid(h,d)
# print(c)
for i in c:
    for k in d:
        plt.scatter(i,k,color='b',marker='+')
plt.show()
