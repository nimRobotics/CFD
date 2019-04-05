from numpy import *
import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt
from gridGen import grid,gridPlot,spacing

def initialize():
    return(zeros((nx+1,ny+1)),zeros((nx+1,ny+1)))

def bcsApply(wMat,pMat):
    # applying bc on omega
    for j in range(ny+1):
        wMat[j,0] = 2*(pMat[0,j]-pMat[1,j])/(dx[0]**2)    # left wall
        wMat[j,nx] = 2*(pMat[nx,j]-pMat[nx-1,j])/(dx[0]**2) # right wall
    for i in range(nx+1):
        wMat[0,i] = 2*(pMat[i,0]-pMat[i,1])/(dy[0]**2)    # bottom wall
        wMat[ny,i] = 2*(pMat[i,ny]-pMat[i-1,ny])/(dy[0]**2)-(2*U)/dy[0] # top wall

    # applying BC, psi is zero at all the boundaries
    pMat[0,:]=0    # left wall
    pMat[nx,:]=0   # right wall
    pMat[:,0]=0    # bottom wall
    pMat[:,ny]=0   # top wall
    # print(wMat)
    return(np.flipud(wMat),np.flipud(pMat))

def calculation():
    a,b = initialize()
    wMat,pMat = bcsApply(a,b)
    # print("BC applied",bcsApply(a,b))

    nIteration = 0
    while nIteration<5:
        nIteration = nIteration+1
        for i in range(1,nx):
            for j in range(1,ny):
                # stream function equation
                a = (((r**2)*(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))/\
                (((r**2))*(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)*(dx[i]+dx[i-1])+\
                (dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i]+dy[i-1])))

                pMat[i,j] = a*(wMat[i,j]+2*(dx[i-1]*pMat[i+1,j]+dx[i]*pMat[i-1,j])/(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)+\
                (1/(r**2))*2*(dy[i-1]*pMat[i,j+1]+dy[i]*pMat[i,j-1])/(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))

                # vorticty formulation, ignoring temporal derivative
                LHS = (Re/r)*(((pMat[i,j+1]-pMat[i,j-1])/(dy[i]+dy[i-1]))*((wMat[i+1,j]-wMat[i-1,j])/(dx[i]+dx[i-1]))-\
                ((pMat[i+1,j]-pMat[i-1,j])/(dx[i]+dx[i-1]))*((wMat[i,j+1]-wMat[i,j-1])/(dy[i]+dy[i-1])))
                b = (((r**2)*(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))/\
                ((dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i]+dy[i-1])+(r**2)*(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)*(dx[i]+dx[i-1])))

                wMat[i,j] = b*(-LHS+2*(dx[i-1]*wMat[i+1,j]+dx[i]*wMat[i-1,j])/(dx[i-1]*dx[i]**2 + dx[i]*dx[i-1]**2)+\
                (1/(r**2))*(dy[i-1]*wMat[i,j+1]+dy[i]*wMat[i,j-1])/(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))

                # resiudals
                res1 = wMat[i,j] \
                + 2*(dx[i-1]*pMat[i+1,j]+dx[i]*pMat[i-1,j]-(dx[i]+dx[i-1])*pMat[i,j])/(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2) \
                + (1/(r**2))*2*(dy[i-1]*pMat[i,j+1]+dy[i]*pMat[i,j-1]-(dy[i]+dy[i-1])*pMat[i,j])/(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)

                res2 = -LHS + 2*(dx[i-1]*wMat[i+1,j]+dx[i]*wMat[i-1,j]-(dx[i]+dx[i-1])*wMat[i,j])/(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)+\
                   (1/(r**2))*2*(dy[i-1]*pMat[i,j+1]+dy[i]*pMat[i,j-1]-(dy[i]+dy[i-1])*pMat[i,j])/(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)
        # print("wMatrix",wMat)
        # print("psiMatrix",pMat)
        # break if resiudal is close to zero
        if abs(res1)<1 and abs(res2)<1:
            break
            print("Itereation finished")

    print("\nres1\n",res1)
    print("\nres2\n",res2)
    print("\nNo of Itereations\n",nIteration)
    print("\nOmega Matrix\n",wMat)
    print("\nPsi Matrix\n",pMat)


nx=10  # elements in x dir
ny=10  # elements in y dir
U = 1  # mormalized plate velocity
Re = 10 # reynolds number
r = 1 # aspect ratio
x,y=grid(nx,ny,2)   # accepts (nx, ny, stretching param)
# gridPlot(x,y)     # plot the grid
dx=spacing(x)  # grid divisions, symmetric
dy=spacing(y)  # grid divisions, symmetric

calculation()
