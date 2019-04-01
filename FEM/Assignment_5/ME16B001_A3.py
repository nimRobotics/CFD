from numpy import *
import scipy.linalg
import numpy as np
from sympy import Symbol
from sympy import *
from numpy import linalg

x=Symbol('x')
def nodes():
    xb=[]
    xa=[]
    for i in range(1,nEle+1):
        xb.append(i*(domainLen/nEle))
        xa.append((i-1)*(domainLen/nEle))
    return(xa,xb)

def kfMatrix():
    k = zeros(nEle+1,nEle+1)
    F = zeros(1,nEle+1)
    xa,xb = nodes()
    for i in range(nEle):
        psi = [(xb[i]-x)/(xb[i]-xa[i]),(x-xa[i])/(xb[i]-xa[i])]
        for m in range(2):
            F[0,m+i] = F[0,m+i] + integrate(f*psi[m] ,(x, xa[i], xb[i]))
            for n in range(2):
                k[m+i,n+i] =k[m+i,n+i]+integrate(a*diff(psi[m],x)*diff(psi[n],x)+c*psi[m]*psi[n],(x, xa[i], xb[i]))
    return(k,F)

def calculation():
    k,F = kfMatrix()
    f_modified = F[1:len(F)-1]
    k_modified = k[1:len(F)-1,1:len(F)-1]
    for i in range(0,len(F)-2):
        f_modified[i]=f_modified[i]-k[i+1,0]*bcs[0]-k[i+1,len(F)-1]*bcs[1]
    sol=[]
    sl=(linalg.inv(np.array(k_modified).astype(np.float64))).dot(f_modified)
    sol.append(bcs[0])
    for i in range(len(sl)):
        sol.append(sl[i])
    sol.append(bcs[1])
    return(sol)


# user inputs
nEle = 30  # number of elements
domainLen = 1  # length of the domain
a=1 # from governing equation
c=-1 # from governing equation
f=-x*x # from governing equation
bcs=[5,15] # can be dirichlet (u1,u2...) or neumann (-adu/dx,adu/dx)
print(calculation())
