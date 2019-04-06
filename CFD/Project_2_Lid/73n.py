from numpy import *
import scipy.linalg
from numpy import linalg as LA
import numpy as np
import matplotlib.pyplot as plt
from gridGen import grid,gridPlot,spacing,uGrid

def plotContour(matrix):
    cs = plt.contourf(matrix,  extend='both')
    cs.cmap.set_over('red')
    cs.cmap.set_under('blue')
    cs.changed()
    plt.show()

def initialize():
    return(zeros((nx+1,ny+1)),zeros((nx+1,ny+1)),zeros((nx+1,ny+1)),zeros((nx+1,ny+1)))

def bcsApply(wMat,psiMat,u,v):
    # applying bc on omega
    for j in range(ny+1):
        wMat[0,j] = 2*(psiMat[0,j]-psiMat[0,j-1])/(dy**2)-(2*U)/dy
    # use central differncing for second row
    print(u)
    for j in range(ny+1):
        wMat[1,j] = -(u[0,j]-u[2,j])/(2*dy)

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
    print("\nBC applied\n wMat",wMat,"\n psiMat\n",psiMat)

    pItMax=1000
    pIt=0
    while pIt<pItMax:
        # find in the velocities u,v from psiMat
        for i in range(1,ny):
            for j in range(1,nx):
                u[i,j] = (psiMat[i,j+1]-psiMat[i,j-1])/(2*dy)
                v[i,j] = -(psiMat[i+1,j]-psiMat[i-1,j])/(2*dx)
        # print("\nu",u,"\nv",v)

        # vorticity stream function relation
        for i in range(1,ny):
            for j in range(1,nx):
                wMatOld = wMat
                RHS=(1/Re)*(((wMat[i+1,j]+wMat[i-1,j])/(dx*dx)) + (1/(r*r))*((wMat[i,j+1]+wMat[i,j-1])/(dy*dy)))
                if i==1 or i==ny-1 or j==1 or j==nx-1:
                    if u[i,j]<0 and v[i,j]>0:
                        wMat[i,j]=((-u[i,j]/dx + v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS-(u[i,j]/dx)*wMat[i+1,j]+(v[i,j]/dy)*wMat[i,j-1])
                    if u[i,j]<0 and v[i,j]<0:
                        wMat[i,j]=((-u[i,j]/dx - v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS-(u[i,j]/dx)*wMat[i+1,j]-(v[i,j]/dy)*wMat[i,j+1])
                    if u[i,j]>0 and v[i,j]<0:
                        wMat[i,j]=((u[i,j]/dx - v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS+(u[i,j]/dx)*wMat[i-1,j]-(v[i,j]/dy)*wMat[i,j+1])
                    if u[i,j]>0 and v[i,j]>0:
                        wMat[i,j]=((u[i,j]/dx + v[i,j]/dy + (2/Re)*(1/(dx*dx)+(1/r**2)*(1/(dy*dy))))**-1)*(RHS+(u[i,j]/dx)*wMat[i-1,j]+(v[i,j]/dy)*wMat[i,j-1])
                # TODO: use second order upwind
                else:
                    if u[i,j]<0 and v[i,j]>0:
                    if u[i,j]<0 and v[i,j]<0:
                    if u[i,j]>0 and v[i,j]<0:
                    if u[i,j]>0 and v[i,j]>0:

                rs1 =np.amax(abs(wMatOld-wMat))
                print(rs1)

        # stream function equation
        for i in range(1,nx):
            for j in range(1,ny):
                psiMatOld=psiMat
                psiMat[i,j] = 0.5*((r*r*(dy**2)*(dx**2))/(r*r*(dy**2)+dx**2))*(wMat[i,j] + (psiMat[i+1,j]+psiMat[i-1,j])/(dx**2)+(1/(r**2))*(psiMat[i,j+1]+psiMat[i,j-1])/(dx**2))
                # resiudals
                # rs2 = wMat[i,j]+(psiMat[i+1,j]-2*psiMat[i,j]+psiMat[i-1,j])/(dx**2)+(1/(r**2))*(psiMat[i,j+1]-2*psiMat[i,j]+psiMat[i,j-1])/(dx**2)
                rs2 = np.amax(abs(psiMatOld-psiMat))
                print("rs2",rs2)
        # break if resiudal is close to zero
        if abs(rs1)<1 and abs(rs2)<1:
            break

        pIt=pIt+1 # counter for psiMat Iterations
    print("\npsiMat at ",pIt," Iteration ",psiMat)
    plotContour(psiMat)

nx=10  # elements in x dir
ny=10  # elements in y dir
U = 1  # mormalized plate velocity
Re = 100 # reynolds number
r = 1 # aspect ratio
H=1
D=1
dx=H/nx
dy=D/ny
calculation()
