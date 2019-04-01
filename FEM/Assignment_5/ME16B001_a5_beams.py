from numpy import *
import scipy.linalg
import numpy as np
from sympy import Symbol
from sympy import *
from numpy import linalg

x=Symbol('x')

def kfMatrix():
    k = zeros(2*nEle+2,2*nEle+2)
    F = zeros(1,2*nEle+2)
    for i in range(nEle):
        psi = [1-3*(x/h[i])**2+2*(x/h[i])**3,-x*(1-(x/h[i]))**2,3*(x/h[i])**2-2*(x/h[i])**3,-x*((x/h[i])**2-(x/h[i]))]
        for m in range(4):
            F[0,m+i] = F[0,m+i] + integrate(f[i]*psi[m] ,(x, xa[i], xb[i]))
            for n in range(4):
                if i==0:
                    k[m+i,n+i] =k[m+i,n+i]+integrate(a*diff(psi[m],x)*diff(psi[n],x),(x, xa[i], xb[i]))
                else:
                    k[m+i+1,n+i+1] =k[m+i,n+i]+integrate(a*diff(psi[m],x)*diff(psi[n],x),(x, xa[i], xb[i]))
    F[:,:]=F[:,:]+bcs[:,:]
    ff=[]
    for i in range(2*nEle+2):
        ff.append(F[i])
    sl=(linalg.inv(np.array(k).astype(np.float64))).dot(ff)
    return(k,F,sl)

# user inputs, smaple for beam with two nodes
nEle = 2  # number of elements
a=1 # AE
f=[10,0] # load per unit length
bcs=zeros(1,2*nEle+2) # shear1, bending1, shear2, bending2, shear3, bending3
h=[1/2,1/2]  # length of the elements
xa=[0,1/2]   # left end coord of elements
xb=[1/2,1]   # right end coord of elements
K,F,sl=kfMatrix()
print("\n K matrix :\n",K)
print("\n F matrix :\n",F)
print("\n Solution matrix :\n",sl)
