from sympy import *
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.patches as mpatches

x = Symbol('x')
x=Symbol('x')
a1=Symbol('a1')
a2=Symbol('a2')
a3=Symbol('a3')
a4=Symbol('a4')

eqn_gal =[]
# trial function
u = a1 + a2*x + a3*x*x + a4*x*x*x
# applying BCs
z=u.subs(a1,solve(u.subs(x,1)-2 , a1)[0])
z=z.subs(a2,solve((x*diff(u,x)+0.5).subs(x,2),a2)[0])
# finding residual
r_gal=diff(x*(diff(z,x)),x)-(2/(x*x))
# galerkian integral wrt weighting functions
eqn_gal.append(integrate(r_gal*x, (x, 1, 2)))
eqn_gal.append(integrate(r_gal*x*x, (x, 1, 2)))
z_gal = solve(eqn_gal, [a3,a4])
print(z_gal)
# plotting the curve
x = np.arange(1, 2., 0.01)
u_gal = 2 - 0.25*(x-1)+(x-1)*(x-3)*z_gal.get(a3)+(x-1)*((x**2)+x-11)*z_gal.get(a4)
plt.plot(x, u_gal , color='b')
plt.xlabel('x', fontsize=16)
plt.ylabel('u', fontsize=16)
plt.show()
