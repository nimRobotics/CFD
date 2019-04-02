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

def cflAnal(i):
    t1=(1/Re)*0.5*((dx[i]**2)*(dy[i]**2))/((dx[i]**2)+(dy[i]**2))
    return()
    # TODO: inplement two conditions

def initialize():
    return(zeros((nx+1,ny+1)),zeros((nx+1,ny+1)))

def bcsApply(wMat,psiMat):
    # applying bc on omega
    for j in range(ny+1):
        wMat[j,0] = 2*(psiMat[0,j]-psiMat[1,j])/(dx[0]**2)    # left wall
        wMat[j,nx] = 2*(psiMat[nx,j]-psiMat[nx-1,j])/(dx[0]**2) # right wall
    for i in range(nx+1):
        wMat[0,i] = 2*(psiMat[i,0]-psiMat[i,1])/(dy[0]**2)    # bottom wall
        wMat[ny,i] = 2*(psiMat[i,ny]-psiMat[i-1,ny])/(dy[0]**2)-(2*U)/dy[0] # top wall

    # applying BC, psi is zero at all the boundaries
    psiMat[0,:]=0    # left wall
    psiMat[nx,:]=0   # right wall
    psiMat[:,0]=0    # bottom wall
    psiMat[:,ny]=0   # top wall
    # print(wMat)
    return(np.flipud(wMat),np.flipud(psiMat))

def calculation():
    a,b = initialize()
    wMat,psiMat = bcsApply(a,b)
    print("\nBC applied\n",bcsApply(a,b))
    # 
    # # calculation of w1
    # for i in range(1,nx):
    #     for j in range(1,ny):
    #         # equations for w(n+1)
    #         LHS = (1/r)*(((psiMat[i,j+1]-psiMat[i,j-1])/(dy[i]+dy[i-1]))*((wMat[i+1,j]-wMat[i-1,j])/(dx[i]+dx[i-1]))-\
    #                      ((psiMat[i+1,j]-psiMat[i-1,j])/(dx[i]+dx[i-1]))*((wMat[i,j+1]-wMat[i,j-1])/(dy[i]+dy[i-1])))
    #
    #         dt=0.01   # first time step
    #         ddx = dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2
    #         ddy = dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2
    #
    #         wMat[i,j]=wMat[i,j]+ dt*(-LHS + (1/Re)*(2*(dx[i-1]*wMat[i+1,j]+dx[i]*wMat[i-1,j]-(dx[i]+dx[i-1])*wMat[i,j])/ddx+\
    #                                      (1/(r**2))*2*(dy[i-1]*wMat[i,j+1]+dy[i]*wMat[i,j-1]-(dy[i]+dy[i-1])*wMat[i,j])/ddy))
    # # print("w at 1",wMat)
    #
    # # while loop for time steps  psi(n) and w(n+1)
    # nt=0
    # while nt<30:
    #     nt=nt+1
    #
    #     # calculation of psi at new time step
    #     for i in range(1,nx):
    #         for j in range(1,ny):
    #             # stream function equation
    #             ddx = dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2
    #             ddy = dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2
    #
    #             a = ((r**2)*ddx*ddy)/((r**2)*ddy*(dx[i]+dx[i-1])+ddx*(dy[i]+dy[i-1]))
    #
    #             psiMat[i,j] = a*(wMat[i,j] + 2*(dx[i-1]*psiMat[i+1,j]+dx[i]*psiMat[i-1,j])/ddx + \
    #                               (1/(r**2))*2*(dy[i-1]*psiMat[i,j+1]+dy[i]*psiMat[i,j-1])/ddy)
    #
    #     # Iteration of psi at the new time step
    #     nIteration = 0
    #     while nIteration<100:
    #         nIteration = nIteration+1
    #         for i in range(1,nx):
    #             for j in range(1,ny):
    #                 psiMat[i,j]=(psiMat[i-1,j]+psiMat[i+1,j]+psiMat[i,j-1]+psiMat[i,j+1])/4
    #
    #                 ddx = dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2
    #                 ddy = dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2
    #                 # resiudal
    #                 res1 = wMat[i,j] + 2*(dx[i-1]*psiMat[i+1,j]+dx[i]*psiMat[i-1,j]-(dx[i]+dx[i-1])*psiMat[i,j])/ddx \
    #                       + (1/(r**2))*2*(dy[i-1]*psiMat[i,j+1]+dy[i]*psiMat[i,j-1]-(dy[i]+dy[i-1])*psiMat[i,j])/ddy
    #
    #         if abs(res1)<(10**(-5)):
    #             break
    #             print("Itereation finished")
    #
    #     print("res",res1)
    #     print("Interations ",nIteration)
    #     print("psimat",psiMat)
    #
    #     # CFL criteria
    #     it=1
    #     dt = cflAnal(it)
    #     it=it+1
    #
    #     #  omega at next time step
    #     for i in range(1,nx):
    #         for j in range(1,ny):
    #             LHS = (1/r)*(((psiMat[i,j+1]-psiMat[i,j-1])/(dy[i]+dy[i-1]))*((wMat[i+1,j]-wMat[i-1,j])/(dx[i]+dx[i-1]))-\
    #                          ((psiMat[i+1,j]-psiMat[i-1,j])/(dx[i]+dx[i-1]))*((wMat[i,j+1]-wMat[i,j-1])/(dy[i]+dy[i-1])))
    #
    #             ddx = dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2
    #             ddy = dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2
    #
    #             wMat[i,j]=wMat[i,j]+ dt*(-LHS + (1/Re)*(2*(dx[i-1]*wMat[i+1,j]+dx[i]*wMat[i-1,j]-(dx[i]+dx[i-1])*wMat[i,j])/ddx+\
    #                                          (1/(r**2))*2*(dy[i-1]*wMat[i,j+1]+dy[i]*wMat[i,j-1]-(dy[i]+dy[i-1])*wMat[i,j])/ddy))
    #
    #     # TODO: implement a break here for steady state
    #
    # # calculation psi at last time step
    # nIteration = 0
    # while nIteration<100:
    #     nIteration = nIteration+1
    #     for i in range(1,nx):
    #         for j in range(1,ny):
    #             psiMat[i,j]=(psiMat[i-1,j]+psiMat[i+1,j]+psiMat[i,j-1]+psiMat[i,j+1])/4
    #
    #             ddx = dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2
    #             ddy = dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2
    #             # resiudal
    #             res1 = wMat[i,j] + 2*(dx[i-1]*psiMat[i+1,j]+dx[i]*psiMat[i-1,j]-(dx[i]+dx[i-1])*psiMat[i,j])/ddx \
    #                   + (1/(r**2))*2*(dy[i-1]*psiMat[i,j+1]+dy[i]*psiMat[i,j-1]-(dy[i]+dy[i-1])*psiMat[i,j])/ddy
    #
    #     if abs(res1)<(10**(-5)):
    #         break
    #         print("Itereation finished")
    #
    # print("res",res1)
    # print("Interations ",nIteration)
    # print("psimat",psiMat)
    # print("wMat",wMat)

    # plotContour(np.flipud(wMat))
    # plotContour(psiMat)

nx=10  # elements in x dir
ny=10  # elements in y dir
U = 1  # mormalized plate velocity
Re = 10 # reynolds number
r = 1 # aspect ratio

x,y=grid(nx,ny,2)   # accepts (nx, ny, stretching param)
dx=spacing(x)  # grid divisions, symmetric
dy=spacing(y)  # grid divisions, symmetric
# print(dx)
# gridPlot(x,y)     # plot the grid
# dx,dy=uGrid(nx,ny)  # for non uniform grid

calculation()
