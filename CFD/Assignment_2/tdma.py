from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from numpy import *
from array import *
from scipy.sparse import *

def TDMA(T,d, n):
	a=[]
	b=[]
	c=[]
	cs=[]
	ds=[]
	u=[]

	for i in range(len(T)):
		for j in range(len(T[0])):
			if (j-i)==1:
				c.append(T[i][j])
			if (i-j==1):
				a.append(T[i][j])
			if i==j:
				b.append(T[i][j])
	# acccomodating n nodes
	tb=b[len(b)-1]
	ta=a[len(a)-1]
	a.pop()
	b.pop()
	for i in range(n-4):
		a.append(a[len(a)-1])
		b.append(b[len(b)-1])
		c.append(c[len(c)-1])
		d.append(d[len(d)-1])
	a.append(ta)
	b.append(tb)
	T = diags([b,a,c], [0,-1, 1]).todense()
	print(T)
	for i in range(len(d)):
		u.append(0)

	for i in range(len(c)):
		if i==0:
			cs.append(c[i]/b[i])
		elif i!=0:
			cs.append(c[i]/(b[i]-a[i-1]*cs[i-1]))

	for i in range(len(d)):
		if i==0:
			ds.append(d[i]/b[i])
		elif i!=0:
			ds.append((d[i]-a[i-1]*ds[i-1])/(b[i]-a[i-1]*cs[i-1]))

	for i in range(len(d)-1,-1,-1):
		if i == len(d)-1:
			u[i]=ds[i]
		elif i!= len(d)-1:
			u[i]=ds[i]-cs[i]*u[i+1]
	return(u)
# Matrix equation AX=B, n is the number of nodes
A = [[10, 5, 0,0,0],[5,15,5,0,0],[0,5,15,5,0],[0,0,5,15,5],[0,0,0,5,10]]
B=[1100,100,100,100,100]
n=10
print(TDMA(A,B,n))
