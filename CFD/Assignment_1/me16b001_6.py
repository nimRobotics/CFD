import numpy as np 
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import math 

x=6
y=np.exp(0.3*x)*np.log(x)-2-x
lx = []
ly = []
error = []
exact = []
i = 0
itera = []
while y > 0.00001:
	y=np.exp(0.3*x)*np.log(x)-2-x
	dy=0.3*np.exp(0.3*x)*np.log(x)+(1/x)*np.exp(0.3*x) - 1
	xn = x - (y/dy)
	lx.append(x)
	x = xn
	ly.append(y)
	i = i+1
print("x = ",x)

# exact solution 4.89389362525
for x in range(i):
	exact.append(4.89389362525)
	itera.append(x)
    
for x in range(i):
	error.append(abs(lx[x]-exact[x]))
# potting the points
plt.scatter(itera, error, marker='o')
plt.ylabel('Absolute error', fontsize=16)
plt.xlabel('Iteration number', fontsize=16)

t = Table([itera, error], names=('Iteration number', 'Absolute error'))
print(t)
plt.show()
