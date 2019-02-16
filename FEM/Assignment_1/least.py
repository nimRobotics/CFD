from sympy import *
import matplotlib.pyplot as plt
import numpy as np 
import matplotlib.patches as mpatches

a_lea = Symbol('a')
b_lea = Symbol('b')
x_lea = Symbol('x')
eqn_lea =[]
r_lea = -0.25+4*(x_lea-1)*a_lea+3*(3*(x_lea**2)-4)*b_lea-(2/(x_lea**2))
# print(diff(r_lea,a_lea))
eqn_lea.append(integrate(r_lea*diff(r_lea,a_lea), (x_lea, 1, 2)))
eqn_lea.append(integrate(r_lea*diff(r_lea,b_lea), (x_lea, 1, 2)))
# print(eqn_lea)
z_lea = solve(eqn_lea, [a_lea,b_lea])
# print(z.get(a))
# print(z.get(b))
x = np.arange(1, 2., 0.01)
u_lea = 2 - 0.25*(x-1)+(x-1)*(x-3)*z_lea.get(a_lea)+(x-1)*((x**2)+x-11)*z_lea.get(b_lea)
plt.plot(x, u_lea , color='b')
plt.xlabel('u', fontsize=16)
plt.ylabel('x', fontsize=16)
plt.show()