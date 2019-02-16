import matplotlib.pyplot as plt
import numpy as np
from sympy import *
############################# col
a_col = Symbol('a')
b_col = Symbol('b')
val_col = [4/3 ,5/3]
eqn_col =[]

def function(x):
	for x in val_col:
		# calculating residuals
		eqn_col.append(-0.25+4*(x-1)*a_col+3*(3*(x**2)-4)*b_col-(2/(x**2)))
	return(eqn_col)

z_col = solve(function(val_col), [a_col,b_col])
############################# sub
a_sub = Symbol('a')
b_sub = Symbol('b')
x_sub = Symbol('x_sub')
eqn_sub =[]
eqn_sub.append(integrate(-0.25+4*(-1)*a_sub+3*(3*(x_sub**2)-4)*b_sub-(2/(x_sub**2)), (x_sub, 1, 1.5)))
eqn_sub.append(integrate(-0.25+4*(x_sub-1)*a_sub+3*(3*(x_sub**2)-4)*b_sub-(2/(x_sub**2)), (x_sub, 1.5, 2)))
z_sub = solve(eqn_sub, [a_sub,b_sub])

############################# least square
a_lea = Symbol('a')
b_lea = Symbol('b')
x_lea = Symbol('x')
eqn_lea =[]
r_lea = -0.25+4*(x_lea-1)*a_lea+3*(3*(x_lea**2)-4)*b_lea-(2/(x_lea**2))
eqn_lea.append(integrate(r_lea*diff(r_lea,a_lea), (x_lea, 1, 2)))
eqn_lea.append(integrate(r_lea*diff(r_lea,b_lea), (x_lea, 1, 2)))
z_lea = solve(eqn_lea, [a_lea,b_lea])
############################# galerkian method
a_gal = Symbol('a')
b_gal = Symbol('b')
x_gal = Symbol('x')
eqn_gal =[]
r_gal = -0.25+4*(x_gal-1)*a_gal+3*(3*(x_gal**2)-4)*b_gal-(2/(x_gal**2))
eqn_gal.append(integrate(r_gal*x_gal, (x_gal, 1, 2)))
eqn_gal.append(integrate(r_gal*x_gal*x_gal, (x_gal, 1, 2)))
z_gal = solve(eqn_gal, [a_gal,b_gal])
# _____________________________________________

x = np.arange(1, 2, 0.01)

u_col = 2. - (0.25)*(x-1)+(x-1)*(x-3)*z_col.get(a_col)+(x-1)*(x**2+x-11)*z_col.get(b_col)
u_exact = (2/x)+0.5*np.log(x)
u_sub = 2 - 0.25*(x-1)+(x-1)*(x-3)*z_sub.get(a_sub)+(x-1)*(x**2+x-11)*z_sub.get(b_sub)
u_lea = 2 - 0.25*(x-1)+(x-1)*(x-3)*z_lea.get(a_lea)+(x-1)*((x**2)+x-11)*z_lea.get(b_lea)
u_gal = 2 - 0.25*(x-1)+(x-1)*(x-3)*z_gal.get(a_gal)+(x-1)*((x**2)+x-11)*z_gal.get(b_gal)

plt.plot(x, u_exact, label="Exact")
plt.plot(x, u_col, label="Collocation")
plt.plot(x, u_sub, label="Subdomain")
plt.plot(x, u_lea, label="Least Squares", color="black")
plt.plot(x, u_gal, label="Galerian", color ="m")

plt.legend()
plt.xlabel('x', fontsize=16)
plt.ylabel('Solution, u', fontsize=16)
plt.show()
