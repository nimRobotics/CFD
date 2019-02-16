from sympy import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

a_sub = Symbol('a')
b_sub = Symbol('b')
x_sub = Symbol('x_sub')
eqn_sub =[]
eqn_sub.append(integrate(-0.25+4*(-1)*a_sub+3*(3*(x_sub**2)-4)*b_sub-(2/(x_sub**2)), (x_sub, 1, 1.5)))
eqn_sub.append(integrate(-0.25+4*(x_sub-1)*a_sub+3*(3*(x_sub**2)-4)*b_sub-(2/(x_sub**2)), (x_sub, 1.5, 2)))
z_sub = solve(eqn_sub, [a_sub,b_sub])
# print(z_sub.get(a_sub))
# print(z_sub.get(b_sub))
x_sub = np.arange(1, 2., 0.01)
u_sub = 2 - 0.25*(x_sub-1)+(x_sub-1)*(x_sub-3)*z_sub.get(a_sub)+(x_sub-1)*(x_sub**2+x_sub-11)*z_sub.get(b_sub)
plt.plot(x_sub, u_sub , color='b', label="Subdomain")

plt.xlabel('x_sub', fontsize=16)
plt.ylabel('u', fontsize=16)
plt.show()
