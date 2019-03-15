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
# function to find k and f for Quadratic elements
def kfmatixQad(nEle,domainLen,a,c,f):
    xb=[]
    xa=[]
    xc=[]
    for i in range(1,nEle+1):
        xb.append(i*(domainLen/nEle))
        xa.append((i-1)*(domainLen/nEle))
    for i in range(nEle):
        xc.append(0.5*xa[i]+0.5*xb[i])
    # print('xa',xa)
    # print('xb',xb)
    # print('xc',xc)
    k=[]
    Ftemp=[]
    # for ith element  # NOTE: [[element 1 k's],[element 2 k's], ...]
    for i in range(nEle):
        psi_1 = ((x-xc[i])*(x-xb[i]))/((xa[i]-xc[i])*(xa[i]-xb[i]))
        psi_2 = ((x-xa[i])*(x-xb[i]))/((xc[i]-xa[i])*(xc[i]-xb[i]))
        psi_3 = ((x-xa[i])*(x-xc[i]))/((xb[i]-xa[i])*(xb[i]-xc[i]))
        # print(psi_1)
        # print(psi_2)
        k.append([])
        Ftemp.append(integrate(f*psi_1 ,(x, xa[i], xb[i])))
        Ftemp.append(integrate(f*psi_2 ,(x, xa[i], xb[i])))
        Ftemp.append(integrate(f*psi_3 ,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_1,x)*diff(psi_1,x)+c*psi_1*psi_1,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_1,x)*diff(psi_2,x)+c*psi_1*psi_2,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_1,x)*diff(psi_3,x)+c*psi_1*psi_3,(x, xa[i], xb[i])))

        k[i].append(integrate( a*diff(psi_2,x)*diff(psi_1,x)+c*psi_2*psi_1,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_2,x)*diff(psi_2,x)+c*psi_2*psi_2,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_2,x)*diff(psi_3,x)+c*psi_2*psi_3,(x, xa[i], xb[i])))

        k[i].append(integrate( a*diff(psi_3,x)*diff(psi_1,x)+c*psi_3*psi_1,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_3,x)*diff(psi_2,x)+c*psi_3*psi_2,(x, xa[i], xb[i])))
        k[i].append(integrate( a*diff(psi_3,x)*diff(psi_3,x)+c*psi_3*psi_3,(x, xa[i], xb[i])))

    F=[]
    F.append(Ftemp[0])
    F.append(Ftemp[1])
    t=0
    i=0
    while t<nEle-1:
        F.append(Ftemp[i+2]+Ftemp[i+3])
        F.append(Ftemp[i+4])
        i=i+3
        t=t+1
    F.append(Ftemp[len(Ftemp)-1])

    # print('k = ',k)
    # print('Ftemp',Ftemp)
    # print('F = ',F)
    # print('len f',len(F))
    # three diagonals of the quadiagonal matrix
    diagA=[]
    diagB=[]
    diagC=[]
    # for three element 4*4 k matrix
    for i in range(nEle):
        diagA.append(k[i][1])
        diagA.append(k[i][5])
    # print('diagA',diagA)

    for i in range(nEle):
        diagC.append(k[i][2])
        if i!=nEle-1:
            diagC.append(0.0)
    # print('diagC',diagC)

    diagB.append(k[0][0])
    diagB.append(k[0][4])
    for i in range(nEle-1):
        diagB.append(k[i][8]+k[i+1][0])
        diagB.append(k[i+1][4])
    diagB.append(k[nEle-1][8])

    # print('diagB',diagB)
    diagA=np.array(diagA, dtype=np.float64)
    diagB=np.array(diagB, dtype=np.float64)
    diagC=np.array(diagC, dtype=np.float64)

    K = np.array(diags([diagB,diagA,diagA,diagC,diagC], [0,-1, 1,-2, 2]).todense() )
    return(K,F)
# function to find Q at the left and right ends
def solQ(k,f,nEle,bcs,method):
    bc1=Symbol('bc1')
    bc2=Symbol('bc2')
    q1=Symbol('q1')
    q2=Symbol('q2')
    c1=0
    c2=0
    if method=='linear':
        for i in range(nEle+1):
            c1=c1+k[0][i]
            c2=c2+k[nEle][i]
        q1_eqn = c1*bc1-f[0]+q1
        q2_eqn = c2*bc2-f[nEle]-q2
    elif method=='Quadratic':
        for i in range(2*nEle+1):
            c1=c1+k[0][i]
            c2=c2+k[nEle][i]
        q1_eqn = c1*bc1-f[0]+q1
        q2_eqn = c2*bc2-f[2*nEle]-q2
    if len(bcs)==0:
        return(solve(q1_eqn,q1),solve(q2_eqn,q2))
    elif len(bcs)==2:
        return(solve(q1_eqn,q1)[0].subs(bc1,bcs[0]),solve(q2_eqn,q2)[0].subs(bc2,bcs[1]))

# user inputs
domainLen=1
# number of elements
nEle=3
# problem data
a=1
c=-1
f=-x*x
method='Both' # accepts 'Quadratic','linear','Both'
# dirchlet boundary conditions values at left and right end
bcs=[2,3] # assumed boundary values [left end,right end]

if method=='Quadratic':
    kq,fq=kfmatixQad(nEle,domainLen,a,c,f)
    q1q,q2q=solQ(kq,fq,nEle,bcs,'Quadratic')
    fq[0]=fq[0]+q1q
    fq[len(fq)-1]=fq[len(fq)-1]+q2q
    sq=linalg.inv(kq).dot(fq)
    print('Quadratic Solution, U : \n',sq)

elif method=='linear':
    k,f=kfmatixLin(nEle,domainLen,a,c,f)
    q1,q2 = solQ(k,f,nEle,bcs,'linear')
    f[0]=f[0]+q1
    f[nEle]=f[nEle]+q2
    sl=linalg.inv(k).dot(f)
    print('Linear Solution, U : \n',sl)

elif method=='Both':
    kq,fq=kfmatixQad(nEle,domainLen,a,c,f)
    q1q,q2q=solQ(kq,fq,nEle,bcs,'Quadratic')
    fq[0]=fq[0]+q1q
    fq[len(fq)-1]=fq[len(fq)-1]+q2q
    sq=linalg.inv(kq).dot(fq)
    print('Quadratic Solution, U : \n',sq)
    k,f=kfmatixLin(nEle,domainLen,a,c,f)
    q1,q2 = solQ(k,f,nEle,bcs,'linear')
    f[0]=f[0]+q1
    f[nEle]=f[nEle]+q2
    sl=linalg.inv(k).dot(f)
    print('\nLinear Solution, U : \n',sl)
