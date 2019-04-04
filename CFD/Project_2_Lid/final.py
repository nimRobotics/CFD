from numpy import *
import scipy.linalg
import numpy as np
import matplotlib.pyplot as plt
from gridGen import grid,gridPlot,spacing,uGrid

def plotContour(matrix):
    cs = plt.contourf(matrix,  extend='both')
    cs.cmap.set_over('red')
    cs.cmap.set_under('blue')
    cs.changed()
    plt.show()

def cflAnal(u,v,tIt):
    # applying cfl criteria
    dt1=[]
    for i in range(1,ny):
        for j in range(1,nx):
            dt1.append(abs((Re/2)*((1/(dx[i]**2))+(1/(r*r*dy[j]**2)))**(-1)))
    if tIt==0:
        return(min(dt1))
    else:
        dt2=[]
        for i in range(1,ny):
            for j in range(1,nx):
                dt2.append(abs((u[i,j]/dx[i] + v[i,j]/dy[j])**(-1)))
        return(0.01*min(min(dt1),min(dt2)))
    # return(0.001)

def initialize():
    return(zeros((nx+1,ny+1)),zeros((nx+1,ny+1)),zeros((nx+1,ny+1)),zeros((nx+1,ny+1)))

def bcsApply(wMat,psiMat,u,v):
    # applying bc on omega
    for j in range(ny+1):
        wMat[0,j] = 2*(psiMat[0,j]-psiMat[0,j-1])/(dy[0]**2)-(2*U)/dy[0]
    # use central differncing for second row
    print(u)
    for j in range(ny+1):
        wMat[1,j] = -(u[0,j]-u[2,j])/(dy[0]+dy[1])

    # applying BC, psi is zero at all the boundaries
    psiMat[0,:]=0    # left wall
    psiMat[nx,:]=0   # right wall
    psiMat[:,0]=0    # bottom wall
    psiMat[:,ny]=0   # top wall
    return(wMat,psiMat)

def calculation():
    a,b,u,v = initialize()
    u[0,:]=1  #  top layer with vel = U
    wMat,psiMat = bcsApply(a,b,u,v)
    print("\nBC applied\n wmatrix",wMat,"\n psiMat\n",psiMat)

    tItMax=10
    tIt=0
    while tIt<tItMax:
        # find in the velocities u,v from psiMat
        for i in range(1,ny):
            for j in range(1,nx):
                u[i,j] = (psiMat[i,j+1]-psiMat[i,j-1])/(dy[i]+dy[i-1])
                v[i,j] = -(psiMat[i+1,j]-psiMat[i-1,j])/(dx[i]+dx[i-1])
        # print("\nu",u,"\nv",v)

        dt = cflAnal(u,v,tIt)
        print("\ntime step : ",dt)
        # calculation of wMat at time n+1
        for i in range(1,nx):
            for j in range(1,ny):
                LHS = u[i,j]*((wMat[i+1,j]-wMat[i-1,j])/(dx[i]+dx[i-1])) + v[i,j]*((wMat[i,j+1]-wMat[i,j-1])/(dy[i]+dy[i-1]))
                ddx = dx[i-1]*(dx[i]**2)+dx[i]*(dx[i-1]**2)
                ddy = dy[i-1]*(dy[i]**2)+dy[i]*(dy[i-1]**2)

                wMat[i,j]=wMat[i,j]+ dt*(-LHS + (1/Re)*(2*(dx[i-1]*wMat[i+1,j]+dx[i]*wMat[i-1,j]-(dx[i]+dx[i-1])*wMat[i,j])/ddx+\
                                             (1/(r**2))*2*(dy[i-1]*wMat[i,j+1]+dy[i]*wMat[i,j-1]-(dy[i]+dy[i-1])*wMat[i,j])/ddy))

                print(LHS)
        print("\nwMat at time ",tIt," step ",wMat)

        # calculation of psiMat at time n+1
        # Iteration for psiMat
        pItMax=100
        pIt=0
        while pIt<pItMax:
            for i in range(1,nx):
                for j in range(1,ny):
                    # stream function equation
                    ddx = dx[i-1]*(dx[i]**2)+dx[i]*(dx[i-1]**2)
                    ddy = dy[i-1]*(dy[i]**2)+dy[i]*(dy[i-1]**2)

                    a = 0.5*((r**2)*ddx*ddy)/((r**2)*ddy*(dx[i]+dx[i-1])+ddx*(dy[i]+dy[i-1]))

                    psiMat[i,j] = a*( wMat[i,j] + (2*(dx[i-1]*psiMat[i+1,j]+dx[i]*psiMat[i-1,j])/ddx) + \
                                       ((1/(r**2))*2*(dy[i-1]*psiMat[i,j+1]+dy[i]*psiMat[i,j-1])/ddy))
                    # resiudal
                    res1 = wMat[i,j] + (2*(dx[i-1]*psiMat[i+1,j]+dx[i]*psiMat[i-1,j]-(dx[i]+dx[i-1])*psiMat[i,j])/ddx) \
                          + ((1/(r**2))*2*(dy[i-1]*psiMat[i,j+1]+dy[i]*psiMat[i,j-1]-(dy[i]+dy[i-1])*psiMat[i,j])/ddy)

            if abs(res1)<(10**(-6)):
                break
                print("Iteration finished")

            pIt=pIt+1 # counter for psiMat Iterations
        print("\npsiMat at ",tIt," step ",psiMat)

        tIt=tIt+1 # counter for time steps
    plotContour(psiMat)


nx=10  # elements in x dir
ny=10  # elements in y dir
U = 1  # mormalized plate velocity
Re = 100 # reynolds number
r = 1 # aspect ratio

x,y=grid(nx,ny,3)   # accepts (nx, ny, stretching param)
dx=spacing(x)  # grid divisions, symmetric
dy=spacing(y)  # grid divisions, symmetric
# print(dx)
# gridPlot(x,y)     # plot the grid
# dx,dy=uGrid(nx,ny)  # for non uniform grid

calculation()
# print(cflAnal(5,4,zeros((nx+1,ny+1)),0))
