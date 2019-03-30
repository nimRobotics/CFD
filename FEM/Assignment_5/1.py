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
    if bcType=="dirichlet":
        f_modified = F[1:len(F)-1]
        k_modified = k[1:len(F)-1,1:len(F)-1]
        for i in range(0,len(F)-2):
            f_modified[i]=f_modified[i]-k[i+1,0]*bcs[0]-k[i+1,len(F)-1]*bcs[1]
    elif bcType=="neumann":
        F[0]=F[0]+bcs[0]
        F[len(F)-1]=F[len(F)-1]+bcs[1]

    if bcType=="dirichlet":
        sol=[]
        sl=(linalg.inv(np.array(k_modified).astype(np.float64))).dot(f_modified)
        sol.append(bcs[0])
        for i in range(len(sl)):
            sol.append(sl[i])
        sol.append(bcs[1])
        return(sol)
    elif bcType=="neumann":
        sl=(np.array(F).astype(np.float64)).dot(linalg.inv(np.array(k).astype(np.float64)))
        return(sl)

nEle = 3
domainLen = 1
a=1
c=-1
f=-x*x
bcType = "neumann"
bcs=[5,15]
print(calculation())
