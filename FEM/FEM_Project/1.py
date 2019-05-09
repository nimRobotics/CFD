from numpy import *
import numpy as np
from collections import OrderedDict
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
    tmp = OrderedDict()
    for point in zip(x, y):
        tmp.setdefault(point[:2], point)
    mypoints = tmp.values()
    # print(mypoints)
    return(len(mypoints))

# assembeled stiffness Matrix
def assMat():
    kGloE = kEle()
    print(kGloE)
    kAss = zeros((uCoord()*2,uCoord()*2))




x1coord=[-5,5]
y1coord=[0,0]
x2coord=[0,0]
y2coord=[-8.66,-8.66]
E = [10000,10000]
A = [0.1,0.1]
print(kEle())
