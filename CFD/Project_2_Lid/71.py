from numpy import *
import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt
from gridGen import grid,gridPlot,spacing

def initialize():
    return(zeros((nx+1,ny+1)),zeros((nx+1,ny+1)))

def bcsApply(wMatrix,psiMatrix):
    # applying bc on omega
    for j in range(ny+1):
        wMatrix[j,0] = 2*(psiMatrix[0,j]-psiMatrix[1,j])/(dx[0]**2)    # left wall
        wMatrix[j,nx] = 2*(psiMatrix[nx,j]-psiMatrix[nx-1,j])/(dx[0]**2) # right wall
    for i in range(nx+1):
        wMatrix[0,i] = 2*(psiMatrix[i,0]-psiMatrix[i,1])/(dy[0]**2)    # bottom wall
        wMatrix[ny,i] = 2*(psiMatrix[i,ny]-psiMatrix[i-1,ny])/(dy[0]**2)-(2*U)/dy[0] # top wall

    # applying BC, psi is zero at all the boundaries
    psiMatrix[0,:]=0    # left wall
    psiMatrix[nx,:]=0   # right wall
    psiMatrix[:,0]=0    # bottom wall
    psiMatrix[:,ny]=0   # top wall
    # print(wMatrix)
    return(np.flipud(wMatrix),np.flipud(psiMatrix))

def calculation():
    a,b = initialize()
    wMatrix,psiMatrix = bcsApply(a,b)
    # print(bcsApply(a,b))

    nIteration = 0
    while nIteration<200:
        nIteration = nIteration+1
        for i in range(1,nx):
            for j in range(1,ny):
                psiMatrix[i,j] =0.01+ 0.5*(((dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))/((dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)*(dx[i]+dx[i-1])+(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i]+dy[i-1])))*(wMatrix[i,j]+2*(dx[i-1]*psiMatrix[i+1,j]+dx[i]*psiMatrix[i-1,j])/(dx[i-1]*dx[i]**2)+2*(dy[i-1]*psiMatrix[i,j+1]+dy[i]*psiMatrix[i,j-1])/(dy[i-1]*dy[i]**2))


                res1 = wMatrix[i,j] \
                + 2*(dx[i-1]*psiMatrix[i+1,j]+dx[i]*psiMatrix[i-1,j]-(dx[i]+dx[i-1])*psiMatrix[i,j])/(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2) \
                + 2*(dy[i-1]*psiMatrix[i,j+1]+dy[i]*psiMatrix[i,j-1]-(dy[i]+dy[i-1])*psiMatrix[i,j])/(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)
        # break if resiudal is close to zero
        # if rs1<1 and rs2<1:
        #     break

    print("res1",res1)
    print("Iter",nIteration)
    print(wMatrix)
    print(psiMatrix)

nx=10  # elements in x dir
ny=10  # elements in y dir
U = 1  # mormalized plate velocity
x,y=grid(nx,ny,2)   # accepts (nx, ny, stretching param)
# gridPlot(x,y)     # plot the grid
dx=spacing(x)  # grid divisions, symmetric
dy=spacing(y)  # grid divisions, symmetric

calculation()
