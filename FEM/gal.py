from sympy import *
import matplotlib.pyplot as plt
import numpy as np 
import matplotlib.patches as mpatches

a_gal = Symbol('a')
b_gal = Symbol('b')
x_gal = Symbol('x')
eqn_gal =[]
r_gal = -0.25+4*(x_gal-1)*a_gal+3*(3*(x_gal**2)-4)*b_gal-(2/(x_gal**2))
print(diff(r_gal,a_gal))
eqn_gal.append(integrate(r_gal*x_gal, (x_gal, 1, 2)))
eqn_gal.append(integrate(r_gal*x_gal*x_gal, (x_gal, 1, 2)))
print(eqn_gal)
z_gal = solve(eqn_gal, [a_gal,b_gal])

x = np.arange(1, 2., 0.01)
u_gal = 2 - 0.25*(x-1)+(x-1)*(x-3)*z_gal.get(a_gal)+(x-1)*((x**2)+x-11)*z_gal.get(b_gal)
plt.plot(x, u_gal , color='b')
plt.xlabel('x', fontsize=16)
plt.ylabel('u', fontsize=16)
plt.show()