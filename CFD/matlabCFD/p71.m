clear all
clear
clc

nx=10  # elements in x dir
ny=10  # elements in y dir
U = 1  # mormalized plate velocity
Re = 0.01 # reynolds number
r = 1 # aspect ratio
%x,y=grid(nx,ny,2)   # accepts (nx, ny, stretching param)
# gridPlot(x,y)     # plot the grid
dx=spacing(x)  # grid divisions, symmetric
dy=spacing(y)  # grid divisions, symmetric

wMat,pMat = bcsApply(zeros((nx+1,ny+1)),zeros((nx+1,ny+1)))
    # print(bcsApply(a,b))
%
%    nIteration = 0
%    while nIteration<200:
%        nIteration = nIteration+1
%        for i in range(1,nx):
%            for j in range(1,ny):
%                # stream function equation
%                a = (((r**2)*(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))/((dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)*(dx[i]+dx[i-1])+(r**2)*(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i]+dy[i-1])))
%                pMat[i,j] = a*(wMat[i,j]+2*(dx[i-1]*pMat[i+1,j]+dx[i]*pMat[i-1,j])/(dx[i-1]*dx[i]**2 + dx[i]*dx[i-1]**2)+2*(dy[i-1]*pMat[i,j+1]+dy[i]*pMat[i,j-1])/(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))
%
%                # vorticty formulation, ignoring temporal derivative
%                LHS = (Re/r)*(((pMat[i,j+1]-pMat[i,j-1])/(dy[i]+dy[i-1]))*((wMat[i+1,j]-wMat[i-1,j])/(dx[i]+dx[i-1]))-((pMat[i+1,j]-pMat[i-1,j])/(dx[i]+dx[i-1]))*((wMat[i,j+1]-wMat[i,j-1])/(dy[i]+dy[i-1])))
%                b = (((r**2)*(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))/((dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2)*(dy[i]+dy[i-1])+(r**2)*(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)*(dx[i]+dx[i-1])))
%                wMat[i,j] = b*(LHS+2*(dx[i-1]*wMat[i+1,j]+dx[i]*wMat[i-1,j])/(dx[i-1]*dx[i]**2 + dx[i]*dx[i-1]**20)+(1/(r**2))*(dy[i-1]*wMat[i,j+1]+dy[i]*wMat[i,j-1])/(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2))
%
%                res1 = wMat[i,j] \
%                + 2*(dx[i-1]*pMat[i+1,j]+dx[i]*pMat[i-1,j]-(dx[i]+dx[i-1])*pMat[i,j])/(dx[i-1]*dx[i]**2+dx[i]*dx[i-1]**2) \
%                + 2*(dy[i-1]*pMat[i,j+1]+dy[i]*pMat[i,j-1]-(dy[i]+dy[i-1])*pMat[i,j])/(dy[i-1]*dy[i]**2+dy[i]*dy[i-1]**2)
%        # break if resiudal is close to zero
%        # if rs1<1 and rs2<1:
%        #     break
%
%    print("res1",res1)
%    print("Iter",nIteration)
%    print(LHS)
%    print(wMat)
%    print(pMat)


