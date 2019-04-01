from numpy import *
import scipy.linalg
import numpy as np
from sympy import Symbol
# from sympy import *
from numpy import linalg
pi = 3.141592653589793

# returns transformation Matrix
def transformMatrix(arg):
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
        length.append(((y1coord[i]-y2coord[i])**2+(x1coord[i]-x2coord[i])**2)**(0.5))
        try:
            theta.append((180/pi)*np.arctan(abs(y1coord[i]-y2coord[i])/abs(x1coord[i]-x2coord[i])))
        except ZeroDivisionError:
            theta.append(90)
    return(theta,length)

# returns kMatrix and solution
def calculation():
    angle,length = elementParams()
    kMatrix = zeros((2*nodes,2*nodes))
    localK = zeros((4,4))
    localK[0,0] = 1
    localK[0,2] = -1
    localK[2,0] = -1
    localK[2,2] = 1
    localK[:,:] = (AE/l)*localK[:,:]
    for i in range(3):
        if i==0:
            kMatrix[0:4,0:4] = (1/length[i])*transformMatrix(angle[i])@localK[:,:]@linalg.inv(transformMatrix(angle[i]))
        elif i==1:
            kMatrix[2:6,2:6] = kMatrix[2:6,2:6]+ (1/length[i])*transformMatrix(angle[i])@localK[:,:]@linalg.inv(transformMatrix(angle[i]))
        elif i==2:
            temp= (1/length[i])*transformMatrix(angle[i])@localK[:,:]@linalg.inv(transformMatrix(angle[i]))
            kMatrix[0:2,0:4]=kMatrix[0:2,0:4]+temp[0:2,0:4]
            kMatrix[4:6,2:6]=kMatrix[4:6,2:6]+temp[2:4,0:4]
    return(kMatrix,linalg.inv(kMatrix)@f)

# user inputs for the problem discussed in class
# assuming all the nodes are placed in unform square grid
x1coord=[1,2,2]
y1coord=[1,1,2]
x2coord=[2,2,1]
y2coord=[1,2,1]
# force matrix f= [f1x,f1y,f2x,f2y,f3x,f3y]
f = [0,0,10,-20,0,0]
AE = 10
nodes=3 #number of nodes
l=2.2 #length of the horizontal or vertical element

a,b=calculation()
print("\nKmatrix\n",a)
print("\nSolutionMatrix\n",b)
