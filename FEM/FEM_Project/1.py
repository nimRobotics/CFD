from numpy import *
import numpy as np
from itertools import groupby
# import scipy.linalg

pi = 3.141592653589793

# returns transformation Matrix
def tMatrix(arg):
    mat=zeros((4,4))
    arg=arg*(pi/180)
    mat[0,0]=np.cos(arg)
    mat[0,1]=-np.sin(arg)
    mat[1,0]=np.sin(arg)
    mat[1,1]=np.cos(arg)
    mat[2,2]=np.cos(arg)
    mat[2,3]=-np.sin(arg)
    mat[3,2]=np.sin(arg)
    mat[3,3]=np.cos(arg)
    return(mat)

# returns element angle and length
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
    kGloE = zeros((len(E),4,4))
    arg,l = elementParams()
    localK = zeros((4,4))
    localK[0,0] = 1
    localK[0,2] = -1
    localK[2,0] = -1
    localK[2,2] = 1

    for i in range(len(E)):
        kGloE[i,:,:] = (A[i]*E[i]/l[i])*tMatrix(arg[i])@localK[:,:]@linalg.inv(tMatrix(arg[i]))
        # print("\nGlobal matrix for element "+str(i)+" = \n",kGloE[i,:,:])
    return(kGloE)

# number of unique coordinates
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
    print(kGloE)
    print(cMat)
    # print(cMat[1,:].index(1))
    # print(np.where(cMat == 1))
    # a=-1
    # b=-1
    kAss = zeros((len(uCoord())*2,len(uCoord())*2))
    for e in range(2):
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




x1coord=[-5,0]
y1coord=[0,-8.66]
x2coord=[0,5]
y2coord=[-8.66,0]
# x1coord=[-5,5]
# y1coord=[0,0]
# x2coord=[0,0]
# y2coord=[-8.66,-8.66]
E = [10000,10000]
A = [0.1,0.1]
print(assMat())
