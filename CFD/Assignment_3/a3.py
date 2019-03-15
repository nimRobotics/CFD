from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from array import *
from scipy.sparse import *

#  (r,c) => (n,i) \\ (q,p)  || (time,space)
L=0.24
alpha=0.001
def Temp(scheme,tinter,s,h):
    T=[[]]
    q=tinter/s
    p=L/h
    r=(alpha*s)/(h*h)

    for i in range(int(p)):
        T[0].append(25)
    # explcit scheme
    if scheme==1:
        for n in range(1,int(q)):
            T.append([])
            for i in range(0,int(p-1)):
                if i==0:
                    # asssumption that temperature at node i=-1 is equal to left end
                    T[n].append( r*T[n-1][i] + (1-2*r)*T[n-1][i] + r*T[n-1][i+1] )
                if i==int(p-1):
                    # asssumption that temperature at node i=21 is equal to rightmost end
                    T[n].append(r*T[n-1][i-1] + (1-2*r)*T[n-1][i] + r*T[n-1][i] )
                else:
                    T[n].append(r*T[n-1][i-1] + (1-2*r)*T[n-1][i] + r*T[n-1][i+1] )
            # boundary conditions
            T[n][0]=100
            T[n][int(p-1)]=75
    # if scheme==2:
        # TODO: implemtnt implicit scheme
    # if scheme==3:
        # TODO: implemtnt CN scheme
    return(T)
#  inputs for calling the function (scheme,time interim,time step, space step)
method=1
tin=0.5
s=0.05
h=0.02
print(Temp(method,tin,s,h))
A=Temp(method,tin,s,h)
l=[]
for i in range(int(L/h)):
    l.append(i*h)
for i in range(int(tin/s)):
    plt.plot(l,A[i],label="t="+str(i*s))
    #plt.plot(l,A[i])

plt.legend(loc='center left', bbox_to_anchor=(1, 0.5),ncol=3, fancybox=True, shadow=True,fontsize='small')
plt.xlabel(r'$x (m)$', fontsize=16)
plt.ylabel(r'$ T^{\circ}C $', fontsize=16)
plt.suptitle('Temperature distribution along the plate')
plt.show()
