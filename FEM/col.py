import matplotlib.pyplot as plt
import numpy as np 
from sympy import *

a_col = Symbol('a')
b_col = Symbol('b')
val_col = [4/3 ,5/3]
eqn_col =[]

def function(x):
	for x in val_col:
		# calculating residuals 
		eqn_col.append(-0.25+4*(x-1)*a_col+3*(3*(x**2)-4)*b_col-(2/(x**2)))
	return(eqn_col)

# print(solve(function(val_col), [a_col,b_col]))
# print("eqn",eqn_col)

z_col = solve(function(val_col), [a_col,b_col])
# {a: 2.09925000000000, b: -0.356000000000000}
# print("a ",z.get(a_col))
# print("b ",z.get(b_col))

x = np.arange(1, 2, 0.01)
u_col = 2. - (0.25)*(x-1)+(x-1)*(x-3)*z_col.get(a_col)+(x-1)*(x**2+x-11)*z_col.get(b_col)
plt.plot(u_col, x)
plt.xlabel('u', fontsize=16)
plt.ylabel('x', fontsize=16)
plt.show()