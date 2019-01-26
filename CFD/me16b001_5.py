from astropy.table import Table, Column
import matplotlib.pyplot as plt
import numpy as np 

a=0
b=2
ival=[]
hval=[]
no = []
exact = []
error =[]
def calculate(n):
	h=(b-a)/n
	hval.append(h)
	A=0
	for x in range(n):
		A=A+0.5*h*(1/(1+x*h*x*h) + 1/(1+(x*h+h)*(x*h+h)))
	return(A)

for x in range(1,11):
	ival.append(calculate(10*x))


for x in range(1,11):
	no.append(x)
	exact.append(1.1071487177943273)
# exact solution   1.1071487177943273

for x in range(10):
	error.append(abs(ival[x]-exact[x]))

t = Table([no , hval, ival, exact , error], names=('No','h', 'IntegralApprox' , ' IntegralExact', 'Error'))
print(t)

slope, intercept = np.polyfit(np.log(hval), np.log(error), 1)
print("Slope of the log-log curve is: ",slope)

plt.scatter(np.log(hval), np.log(error))
plt.ylabel('log(h)', fontsize=16)
plt.xlabel('log(error)', fontsize=16)
# function to show the plot
plt.show()