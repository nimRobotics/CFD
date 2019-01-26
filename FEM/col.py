import matplotlib.pyplot as plt
import numpy as np 
from sympy import *

a = Symbol('a')
b = Symbol('b')
val = [4/3 ,5/3]
eqn =[]

def function(x):
	for x in val:
		# calculating residuals 
		eqn.append(-0.25+4*(x-1)*a+3*(3*(x**2)-4)*b-(2/(x**2)))
	return(eqn)

# print(solve(function(val), [a,b]))
# print(eqn)

z = solve(function(val), [a,b])
# {a: 2.09925000000000, b: -0.356000000000000}

x = np.arange(0.1, 5., 0.01)
u = 2 - 0.25*(x-1)+(x-1)*(x-3)*z.get(a)+(x-1)*(x**2+x-11)*z.get(b)
plt.plot(u, x)
plt.xlabel('u', fontsize=16)
plt.ylabel('x', fontsize=16)
plt.show()