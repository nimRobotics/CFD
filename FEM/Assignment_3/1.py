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

# function to find k and f
def kfmatix(nEle,domainLen,a,c,f):
    xb=[]
    xa=[]
    for i in range(1,nEle+1):
        xb.append(i*(domainLen/3))
        xa.append((i-1)*(domainLen/3))
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
        k[i].append(integrate( a*diff(psi_1,x)*diff(psi_1,x)+c*psi_1*psi_2,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_1,x)*diff(psi_2,x)+c*psi_1*psi_2,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_2,x)*diff(psi_1,x)+c*psi_1*psi_2,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_2,x)*diff(psi_2,x)+c*psi_1*psi_2,(x, xa[i], xb[i])))

    F=[]
    F.append(Ftemp[0])
    for i in range(0,len(Ftemp)-2,2):
        F.append(Ftemp[i+1]+Ftemp[i+2])
    F.append(Ftemp[len(Ftemp)-1])

    # print('k = ',k)
    # print('Ftemp',Ftemp)
    # print('F = ',F)
    # three diagonals of the tridiagonal matrix
    diagA=[]
    diagB=[]
    # for three element 4*4 k matrix
    for i in range(nEle):
        diagA.append(k[i][2])
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
def solQ(k,f,nEle,bcs):
    bc1=Symbol('bc1')
    bc2=Symbol('bc2')
    q1=Symbol('q1')
    q2=Symbol('q2')
    c1=0
    c2=0
    for i in range(nEle+1):
        c1=c1+k[0][i]
        c2=c2+k[nEle][i]
    q1_eqn = c1*bc1-f[0]+q1
    q2_eqn = c2*bc2-f[nEle]-q2
    if len(bcs)==0:
        return(solve(q1_eqn,q1),solve(q2_eqn,q2))
    elif len(bcs)==2:
        return(solve(q1_eqn,q1)[0].subs(bc1,bcs[0]),solve(q2_eqn,q2)[0].subs(bc2,bcs[1]))

# length of the domain
domainLen=1
# number of elements
nEle=10
# problem data
a=1
c=-1
f=-x*x

k,f=kfmatix(nEle,domainLen,a,c,f)
print('K matrix : \n',k)
print('f matrix : \n',f)
# boundary conditions
bcs=[2,3]
q1,q2 = solQ(k,f,nEle,bcs)
print(q1,q2)
f[0]=f[0]+q1
f[nEle]=f[nEle]+q2
print('f + Q matrix : \n',f)
print('Solution, U : \n',linalg.inv(k).dot(f))
