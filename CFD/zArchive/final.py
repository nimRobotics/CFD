from numpy import *
import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt
from gridGen import grid,gridPlot,spacing,uGrid

def cflAnal(i,j,psiMat,it):
    # find in the velocities from psiMat
    u = zeros((nx+1,ny+1))
    v = zeros((nx+1,ny+1))
    u[0,:]=1  #  top layer with vel = U
    for i in range(1,ny):
        for j in range(1,nx):
            u[i,j] = (psiMat[i,j+1]-psiMat[i,j-1])/(dy[i]+dy[i-1])
            v[i,j] = (psiMat[i+1,j]-psiMat[i-1,j])/(dx[i]+dx[i-1])
    # print("\nu",u)
    # print("\nv",v)

    # applying cfl criteria
    dt1 = (Re/2)*((1/(dx[i]**2))+(1/(r*r*dy[j]**2)))**(-1)
    if it==0:
        return(dt1)
    else:
        dt2 = 1/((u[i,j]/dx[i] + v[i,j]/dy[j])**(-1))
        return(0.9*min(dt1,dt2))

def initialize():
    return(zeros((nx+1,ny+1)),zeros((nx+1,ny+1)))

def bcsApply(wMat,psiMat):
    # applying bc on omega
    for j in range(ny+1):
        wMat[0,j] = 2*(psiMat[0,j]-psiMat[0,j])/(dy[0]**2)-(2*U)/dy[0]
    # TODO: use central differncing for second row
    # applying BC, psi is zero at all the boundaries
    psiMat[0,:]=0    # left wall
    psiMat[nx,:]=0   # right wall
    psiMat[:,0]=0    # bottom wall
    psiMat[:,ny]=0   # top wall
    return(wMat,psiMat)

def calculation():
    a,b = initialize()
    wMat,psiMat = bcsApply(a,b)
    print("\nBC applied\n wmatrix",wMat,"\n psiMat\n",psiMat)

    tItMax=100
    tIt=0
    while tIt<tItMax:
        dt = cflAnal(psiMat,tIt)

        # calculation of wMat at time n+1
        for i in range(1,nx):
            for j in range(1,ny):
                # equations for w(n+1)

                LHS = (1/r)*(((psiMat[i,j+1]-psiMat[i,j-1])/(dy[i]+dy[i-1]))*((wMat[i+1,j]-wMat[i-1,j])/(dx[i]+dx[i-1]))-\
                             ((psiMat[i+1,j]-psiMat[i-1,j])/(dx[i]+dx[i-1]))*((wMat[i,j+1]-wMat[i,j-1])/(dy[i]+dy[i-1])))
                print(LHS)
                ddx = dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2
                ddy = dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2

                wMat[i,j]=wMat[i,j]+ dt*(-LHS + (1/Re)*(2*(dx[i-1]*wMat[i+1,j]+dx[i]*wMat[i-1,j]-(dx[i]+dx[i-1])*wMat[i,j])/ddx+\
                                             (1/(r**2))*2*(dy[i-1]*wMat[i,j+1]+dy[i]*wMat[i,j-1]-(dy[i]+dy[i-1])*wMat[i,j])/ddy))
        # print("\nwMat at time ",it," step ",wMat)

        # calculation of psiMat at time n+1
        # Iteration for psiMat
        pItMax=100
        pIt=0
        while pIt<pItMax:
            for i in range(1,nx):
                for j in range(1,ny):
                    # stream function equation
                    ddx = dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2
                    ddy = dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2

                    a = ((r**2)*ddx*ddy)/((r**2)*ddy*(dx[i]+dx[i-1])+ddx*(dy[i]+dy[i-1]))

                    psiMat[i,j] = a*(wMat[i,j] + 2*(dx[i-1]*psiMat[i+1,j]+dx[i]*psiMat[i-1,j])/ddx + \
                                      (1/(r**2))*2*(dy[i-1]*psiMat[i,j+1]+dy[i]*psiMat[i,j-1])/ddy)
                    # resiudal
                    res1 = wMat[i,j] + 2*(dx[i-1]*psiMat[i+1,j]+dx[i]*psiMat[i-1,j]-(dx[i]+dx[i-1])*psiMat[i,j])/ddx \
                          + (1/(r**2))*2*(dy[i-1]*psiMat[i,j+1]+dy[i]*psiMat[i,j-1]-(dy[i]+dy[i-1])*psiMat[i,j])/ddy

            if abs(res1)<(10**(-5)):
                break
                print("Itereation finished")

            pIt=pIt+1 # counter for psiMat Iterations

        tIt=tIt+1 # counter for time steps
    print("\npsiMat at ",tIt," step ",psiMat)
    print("\nwMat at ",tIt," step ",wMat)




nx=10  # elements in x dir
ny=10  # elements in y dir
U = 1  # mormalized plate velocity
Re = 10 # reynolds number
r = 1 # aspect ratio

x,y=grid(nx,ny,1)   # accepts (nx, ny, stretching param)
dx=spacing(x)  # grid divisions, symmetric
dy=spacing(y)  # grid divisions, symmetric
# print(dx)
# gridPlot(x,y)     # plot the grid
# dx,dy=uGrid(nx,ny)  # for non uniform grid

calculation()
# print(cflAnal(5,4,zeros((nx+1,ny+1)),0))
