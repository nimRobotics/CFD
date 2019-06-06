from numpy import *
import numpy as np
from itertools import groupby
import matplotlib.pyplot as plt

pi = 3.141592653589793

# returns transformation Matrix
def tMatrix(arg):
    mat=zeros((4,4))  #matrix initialized
    arg=arg*(pi/180)  #degree to radian
    mat[0,0]=np.cos(arg)
    mat[0,1]=-np.sin(arg)
    mat[1,0]=np.sin(arg)
    mat[1,1]=np.cos(arg)
    mat[2,2]=np.cos(arg)
    mat[2,3]=-np.sin(arg)
    mat[3,2]=np.sin(arg)
    mat[3,3]=np.cos(arg)
    return(mat)

# returns angle and length for each truss element
def elementParams():
    theta=[]
    length=[]
    for i in range(len(x1coord)):
        length.append(((y2coord[i]-y1coord[i])**2+(x2coord[i]-x1coord[i])**2)**(0.5))
        try:
            theta.append((180/pi)*np.arctan((y2coord[i]-y1coord[i])/(x2coord[i]-x1coord[i])))
        except ZeroDivisionError:
            theta.append(90)
    return(theta,length)

# global stiffness matrix for each element
def kEle():
    kGloE = zeros((len(E),4,4)) # matrix initialized
    arg,l = elementParams()  # angle and length using elementParams()
    localK = zeros((4,4)) # matrix initialized
    localK[0,0] = 1
    localK[0,2] = -1
    localK[2,0] = -1
    localK[2,2] = 1

    for i in range(len(E)):
        kGloE[i,:,:] = (A[i]*E[i]/l[i])*tMatrix(arg[i])@localK[:,:]@linalg.inv(tMatrix(arg[i]))
        # print("\nGlobal matrix for element "+str(i)+" = \n",kGloE[i,:,:])
    return(kGloE)

# number of unique coordinates/ node coordinates
def uCoord():
    x=x1coord+x2coord
    y=y1coord+y2coord
    keyfunc = lambda p: p[:2]
    mypoints = []
    for k, g in groupby(sorted(zip(x, y), key=keyfunc), keyfunc):
        mypoints.append(list(g)[0])
    # print(mypoints[1][1])
    return(mypoints)

# connectivity Matrix
def conMat():
    uCord=uCoord()
    cMat=zeros((4,len(E)+1))
    cMat[0,0]=1
    cMat[1,0]=2
    cMat[2,0]=3
    cMat[3,0]=4
    for i in range(len(E)):  #element index
        for j in range(len(uCord)):  #node index
            if x1coord[i]==uCord[j][0] and y1coord[i]==uCord[j][1]:
                cMat[0,i+1]=1+j*2
                cMat[1,i+1]=2+j*2
            if x2coord[i]==uCord[j][0] and y2coord[i]==uCord[j][1]:
                cMat[2,i+1]=1+j*2
                cMat[3,i+1]=2+j*2
    # print("connectivity matrix \n",cMat)
    return(cMat)

# assembeled stiffness Matrix
def assMat():
    kGloE=kEle()
    cMat=conMat()
    cMat[:]=cMat[:]-1
    # print(kGloE)
    # print(cMat)
    # print(cMat[1,:].index(1))
    # print(np.where(cMat == 1))
    # a=-1
    # b=-1
    kAss = zeros((len(uCoord())*2,len(uCoord())*2))
    for e in range(len(E)):
        for i in range(len(uCoord())*2):
            for j in range(len(uCoord())*2):
                for k in range(4):
                    if i==cMat[k,e+1]:
                        a=k
                    if j==cMat[k,e+1]:
                        b=k
                # print(a,b)
                # print(int(cMat[b,0]))
                if a>=0 and b>=0:
                    kAss[i,j]=kAss[i,j]+kGloE[e,int(cMat[a,0]),int(cMat[b,0])]
                a=-1
                b=-1
    # print(kAss)
    return(kAss)

# solution
def solution():
    k=assMat()
    delRC =[]
    for i in range(len(uCoord())*2):
        if u[i]==0:
            delRC.append(i)
    k=np.delete(k,delRC, 0)  # deleting rows
    k=np.delete(k,delRC,1)   # deleting cols
    print(k)
    mf=np.delete(f,delRC,0)
    print(mf)
    mf.astype(float64)
    a = np.linalg.inv(k).astype(float)
    b = mf.astype(float)
    sol = a@b
    for i in range(len(uCoord())*2):
        if u[i]=='uk':
            u[i]=sol[0]
            sol=np.delete(sol,0, 0)

    mmf = assMat()@u
    return(u,mmf)

#draw lines
def connectpoints(x1coord, x2coord, y1coord, y2coord,p):
    x1, x2 = x1coord[p], x2coord[p]
    y1, y2 = y1coord[p], y2coord[p]
    plt.scatter([x1,x2],[y1,y2],linewidths=10)
    plt.plot([x1,x2],[y1,y2])
# plot the truss
def trussPlot(x1coord, x2coord, y1coord, y2coord):
    # plot the truss
    for i in range(2):
        connectpoints(x1coord, x2coord, y1coord, y2coord, i)
    plt.show()
#_______________________________________________________________________________
# Users inputs begins
# Note: unknowns are abbrevated as 'uk'
#_____________________________________________________________________________

# force matrix f= [fx1,fy1,fx2,fy2,...]
f = ['uk','uk','uk',-1732,'uk','uk']
# boundary conditions/ node conditions
u = [0,0,0,'uk',0,0]
x1coord=[-5,5]
y1coord=[0,0]
x2coord=[0,0]
y2coord=[-8.66,-8.66]
E = [100000000,100000000]
A = [0.1,0.1]
print(solution())
trussPlot(x1coord, x2coord, y1coord, y2coord)
