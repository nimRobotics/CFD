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

# function to find k and f for linear elements
def kfmatixLin(nEle,domainLen,a,c,f):
    xb=[]
    xa=[]
    for i in range(1,nEle+1):
        xb.append(i*(domainLen/nEle))
        xa.append((i-1)*(domainLen/nEle))
    k=[]
    Ftemp=[]
    # for ith element  # NOTE: [[element 1 k's],[element 2 k's], ...]
    for i in range(nEle):
        psi_1 = (xb[i]-x)/(xb[i]-xa[i])
        psi_2 = (x-xa[i])/(xb[i]-xa[i])
        # print(psi_1)
        # print(psi_2)
        k.append([])
        Ftemp.append(integrate(f*psi_1 ,(x, xa[i], xb[i])))
        Ftemp.append(integrate(f*psi_2 ,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_1,x)*diff(psi_1,x)+c*psi_1*psi_1,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_1,x)*diff(psi_2,x)+c*psi_1*psi_2,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_2,x)*diff(psi_1,x)+c*psi_2*psi_1,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_2,x)*diff(psi_2,x)+c*psi_2*psi_2,(x, xa[i], xb[i])))

    F=[]
    F.append(Ftemp[0])
    for i in range(0,len(Ftemp)-2,2):
        F.append(Ftemp[i+1]+Ftemp[i+2])
    F.append(Ftemp[len(Ftemp)-1])

    # print('k = ',k)
    # print('Ftemp',Ftemp)
    # print('F = ',F)
    # two diagonals of the tridiagonal matrix
    diagA=[]
    diagB=[]
    # for three element 4*4 k matrix
    for i in range(nEle):
        diagA.append(k[i][1])
    # print('diagA',diagA)
    # NOTE: no need for diagC as it will always be same as diagA

    diagB.append(k[0][0])
    for i in range(nEle-1):
        diagB.append(k[i][3]+k[i+1][0])
    diagB.append(k[nEle-1][3])
    # print('diagB',diagB)
    diagA=np.array(diagA, dtype=np.float64)
    diagB=np.array(diagB, dtype=np.float64)

    K = np.array( diags([diagB,diagA,diagA], [0,-1, 1]).todense() )
    return(K,F)

# function to find Q at the left and right ends
def solQ(k,f,nEle,bcs,method):

    k_modified = k[1:len(k[0])-1,1:len(k[0])-1]
    # print(k_modified)
    f_modified = f[1:len(f)-1]
    # print(f_modified)
    f_modified[0]=f_modified[0]-k[1][0]*bcs[0]-k[1][len(k[0])-1]*bcs[1]
    f_modified[1]=f_modified[1]-k[len(k)-2][0]*bcs[0]-k[len(k)-2][len(k[0])-1]*bcs[1]

    # print(f_modified)
    # print(linalg.inv(k_modified).dot(f_modified))
    sl=linalg.inv(k_modified).dot(f_modified)
    sol = []
    sol.append(bcs[0])
    for i in range(len(sl)):
        sol.append(sl[i])
    sol.append(bcs[1])
    # print(sol)
    return(sol)




# user inputs
domainLen=1
nEle=5
a=1
c=-1
f=-x*x
bcs=[0,10]

k,f=kfmatixLin(nEle,domainLen,a,c,f)
print("k",k)
print("\nf",f)
print("\nSolution",solQ(k,f,nEle,bcs,'linear'))

# for i in range(10):
#     nEle=i+3
#     print(nEle)
#     k,f=kfmatixLin(nEle,domainLen,a,c,f)
#     print("\nSolution",solQ(k,f,nEle,bcs,'linear'))
