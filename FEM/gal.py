from sympy import *
import matplotlib.pyplot as plt
import numpy as np 
import matplotlib.patches as mpatches

a = Symbol('a')
b = Symbol('b')
x = Symbol('x')
eqn2 =[]
r = -0.25+4*(x-1)*a+3*(3*(x**2)-4)*b-(2/(x**2))
print(diff(r,a))
eqn2.append(integrate(r*x, (x, 1, 2)))
eqn2.append(integrate(r*x*x, (x, 1, 2)))
print(eqn2)
z = solve(eqn2, [a,b])
# print(z.get(a))
# print(z.get(b))
x = np.arange(0.1, 5., 0.01)
u = 2 - 0.25*(x-1)+(x-1)*(x-3)*z.get(a)+(x-1)*((x**2)+x-11)*z.get(b)
plt.plot(u, x , color='b')
plt.xlabel('u', fontsize=16)
plt.ylabel('x', fontsize=16)
plt.show()