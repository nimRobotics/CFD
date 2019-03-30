from numpy import *
import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt

# fn returns delta x delta y
def uGridPoints():
    return(H/nx,D/ny)

# initialize the omega and the psi matrix
def initialize():
    nEleX = nx+1
    nEleY = ny+1
    wMatrix = ones((nEleX,nEleY))
    psiMatrix = ones((nEleX,nEleY))
    # applying BC, psi is zero at all the boundaries
    wMatrix[0,:]=0
    wMatrix[nx,:]=0
    psiMatrix[0,:]=0
    psiMatrix[nx,:]=0
    psiMatrix[:,0]=0
    psiMatrix[:,ny]=0
    return(wMatrix,psiMatrix)

# main function, performs Iteration on the wMatrix and psiMatrix
def calculation():
    wMatrix,psiMatrix = initialize()
    dx,dy = uGridPoints()
    nIteration = 0

    while nIteration<200:
        nIteration = nIteration+1
        for i in range(1,nx):
            for j in range(1,ny):
                psiMatrix[i,j] = 0.5*((r**2*dy**2*dx**2)/(r**2*dy**2+dx**2))*(wMatrix[i,j] + (psiMatrix[i+1,j]+psiMatrix[i-1,j])/(dx**2)+(1/(r**2))*(psiMatrix[i,j+1]+psiMatrix[i,j-1])/(dx**2))
                LHS  = (Re/r)*(((psiMatrix[i,j+1]-psiMatrix[i,j-1])/(2*dy))*((wMatrix[i+1,j]-wMatrix[i-1,j])/(2*dx)) - ((psiMatrix[i+1,j]-psiMatrix[i-1,j])/(2*dx))*((wMatrix[i,j+1]-wMatrix[i,j-1])/(2*dy)))
                wMatrix[i,j] = 0.5*((r**2*dy**2*dx**2)/(r**2*dy**2+dx**2))*(((wMatrix[i+1,j]+wMatrix[i-1,j])/(dx*dx)) + (1/(r*r))*((wMatrix[i,j+1]+wMatrix[i,j-1])/(dy*dy)) - LHS)

                rs1 = wMatrix[i,j]+(psiMatrix[i+1,j]-2*psiMatrix[i,j]+psiMatrix[i-1,j])/(dx**2)+(1/(r**2))*(psiMatrix[i,j+1]-2*psiMatrix[i,j]+psiMatrix[i,j-1])/(dx**2)
                rs2 = LHS - (((wMatrix[i+1,j]-2*wMatrix[i,j]+wMatrix[i-1,j])/(dx**2))+(1/(r**2))*((wMatrix[i,j+1]-2*wMatrix[i,j]+wMatrix[i,j-1])/(dx**2)))

        # break if resiudal is close to zero
        if rs1<1 and rs2<1:
            break

    print(rs1)
    print(rs2)
    print(nIteration)
    print(wMatrix,psiMatrix)
    return(wMatrix,psiMatrix)

def plotContour(matrix):
    cs = plt.contourf(matrix,  extend='both')
    cs.cmap.set_over('red')
    cs.cmap.set_under('blue')
    cs.changed()
    plt.show()

# TODO: check the BCs in initialize fnc
# TODO: use upwind scheme second order

H=1    # length of the cavity
D=1    # depth of the cavity
r=H/D  # aspect ratio
Re=10 # reynolds number
nx=100 # number of elements in x direction
ny=110 # number of elements in y direction

w,psi =calculation()
plotContour(psi)
plotContour(w)
