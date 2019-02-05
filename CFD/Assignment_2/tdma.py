from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from array import *

def TDMA(T,d):
	a=[]
	b=[]
	c=[]
	cs=[]
	ds=[]
	u=[]

	for i in range(len(d)):
		u.append(0)

	# ith row and jth col
	for i in range(len(T)):
		for j in range(len(T[0])):
			print(T[i][j],  end ="   ")
		print()

	for i in range(len(T)):
		for j in range(len(T[0])):
			if (j-i)==1:
				c.append(T[i][j])
			if (i-j==1):
				a.append(T[i][j])
			if i==j:
				b.append(T[i][j])
	# print("a",a)
	# print("b",b)
	# print("c",c)

	for i in range(len(c)):
		if i==0:
			cs.append(c[i]/b[i])
		elif i!=0:
			cs.append(c[i]/(b[i]-a[i-1]*cs[i-1]))
	# print("cs",cs)

	for i in range(len(d)):
		if i==0:
			ds.append(d[i]/b[i])
		elif i!=0:
			ds.append((d[i]-a[i-1]*ds[i-1])/(b[i]-a[i-1]*cs[i-1]))
	# print("ds",ds)

	for i in range(len(d)-1,-1,-1):
		if i == len(d)-1:
			u[i]=ds[i]
		elif i!= len(d)-1:
			u[i]=ds[i]-cs[i]*u[i+1]
	return(u)

T = [[10, 5, 0,0,0],[5,15,5,0,0],[0,5,15,5,0],[0,0,5,15,5],[0,0,0,5,10]]
B=[1100,100,100,100,100]
lr=[]
# print(TDMA(A,B))


a=[]
b=[]
c=[]
# ith row and jth col
for i in range(len(T)):
	for j in range(len(T[0])):
		print(T[i][j],  end ="   ")
	print()

for i in range(len(T)):
	for j in range(len(T[0])):
		if (j-i)==1:
			c.append(T[i][j])
		if (i-j==1):
			a.append(T[i][j])
		if i==j:
			b.append(T[i][j])
print("a",a)
print("b",b)
print("c",c)

# 10   5   0   0   0
# 5   15   5   0   0
# 0   5   15   5   0
# 0   0   5   15   5
# 0   0   0   5   10
# a [5, 5, 5, 5]
# b [10, 15, 15, 15, 10]
# c [5, 5, 5, 5]

# a [5, 5, 5]
# b [10, 15, 15, 15]
# c [5, 5, 5, 5]

n=5
tb=b[len(b)-1]
ta=a[len(a)-1]
a.pop()
b.pop()
for i in range(n-4):
	a.append(a[len(a)-1])
	b.append(b[len(b)-1])
	c.append(c[len(c)-1])


print("a",a)
print("b",b)
print("c",c)


# for j in range(len(A[0])):
# 	lr.append(A[len(A)-1][j])
# print(lr)
# A.remove(lr)
# # TODO: build the mmatrix
# n=5
# for i in range(n-4):
# 	print(i)
# 	for j in range(len(A[0])-1):
# 		A[j].append(0)
# 	print()
#
# for i in range(n-3):
# 	print(i)
# 	A.append([])
# 	print(len(A))
# 	for j in range(len(A)-1,):
# 		A[j].append(0)
# 	print()
#
# for i in range(len(A)):
# 	for j in range(len(A[0])):
# 		print(A[i][j],  end ="   ")
# 	print()
