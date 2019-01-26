import numpy as np 
import matplotlib.pyplot as plt
from astropy.table import Table, Column
import math 
# print (math.log(14))
# print(np.log(14)) 
# print(np.exp(1))
x=60
y=np.exp(0.3*x)*np.log(x)-2-x
lx = []
ly = []
while y > 0.00001:
	y=np.exp(0.3*x)*np.log(x)-2-x
	dy=0.3*np.exp(0.3*x)*np.log(x)+(1/x)*np.exp(0.3*x) - 1
	xn = x - (y/dy)
	print("x = ",x)
	lx.append(x)
	x = xn
	print("y = ",y)
	ly.append(y)
# potting the points
plt.scatter(lx, ly)
plt.ylabel('y', fontsize=16)
plt.xlabel('x', fontsize=16)
# function to show the plot
plt.show()
t = Table([ls, ly], names=('a', 'b'))
print(t)