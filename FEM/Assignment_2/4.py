from sympy import *
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
import numpy as np

x=Symbol('x')
a1=Symbol('a1')
a2=Symbol('a2')
a3=Symbol('a3')
a4=Symbol('a4')
a5=Symbol('a5')
eqn_gal2 =[]
eqn_gal3 =[]
eqn_gal4 =[]


u2 = a1 + a2*x + a3*x*x
z2=u2.subs(a1,solve(u2.subs(x,1)-2 , a1)[0])
z2=z2.subs(a2,solve((x*diff(u2,x)+0.5).subs(x,2),a2)[0])
# finding residual
r_gal2=diff(x*(diff(z2,x)),x)-(2/(x*x))
eqn_gal2.append(integrate(r_gal2*x*x, (x, 1, 2)))
z_gal2 = solve(eqn_gal2, a3)
print("Equation for Quadratic",eqn_gal2)
print(z_gal2)

u3 = a1 + a2*x + a3*x*x + a4*x*x*x
z3=u3.subs(a1,solve(u3.subs(x,1)-2 , a1)[0])
z3=z3.subs(a2,solve((x*diff(u3,x)+0.5).subs(x,2),a2)[0])
# finding residual
r_gal3=diff(x*(diff(z3,x)),x)-(2/(x*x))
# galerkian integral wrt weighting functions
eqn_gal3.append(integrate(r_gal3*x*x, (x, 1, 2)))
eqn_gal3.append(integrate(r_gal3*x*x*x, (x, 1, 2)))
z_gal3 = solve(eqn_gal3, [a3,a4])
print("Equation for Qubic",eqn_gal3)
print(z_gal3)

u4 = a1 + a2*x + a3*x*x + a4*x*x*x + a5*x*x*x*x
z4=u4.subs(a1,solve(u4.subs(x,1)-2 , a1)[0])
z4=z4.subs(a2,solve((x*diff(u4,x)+0.5).subs(x,2),a2)[0])
# finding residual
r_gal4=diff(x*(diff(z4,x)),x)-(2/(x*x))
# galerkian integral wrt weighting functions
eqn_gal4.append(integrate(r_gal4*x*x, (x, 1, 2)))
eqn_gal4.append(integrate(r_gal4*x*x*x, (x, 1, 2)))
eqn_gal4.append(integrate(r_gal4*x*x*x*x, (x, 1, 2)))
z_gal4 = solve(eqn_gal4, [a3,a4,a5])
print("Equation for Biquadratic",eqn_gal4)
print(z_gal4)

# plotting the curve
u_gal2=[]
u_gal3=[]
u_gal4=[]
u_exact=[]
y = []
i=1
while i<=2:
    y.append(i)
    u_gal2.append((z2.subs(a3,z_gal2.get(a3))).subs(x,i))
    u_gal3.append(((z3.subs(a3,z_gal3.get(a3))).subs(a4,z_gal3.get(a4))).subs(x,i))
    u_gal4.append((((z4.subs(a3,z_gal4.get(a3))).subs(a4,z_gal4.get(a4))).subs(a5,z_gal4.get(a5))).subs(x,i))
    u_exact.append((2/i)+0.5*np.log(i))
    i+=0.1

plt.plot(y, u_gal2, label="Quadratic")
plt.plot(y, u_gal3, label="Qubic")
plt.plot(y, u_gal4, label="Biquadratic")
plt.plot(y, u_exact, label="Exact")
plt.legend()
plt.xlabel('x', fontsize=16)
plt.ylabel('u', fontsize=16)
plt.suptitle('Solutions with different Trial functions')
print("It can be observed from the graph that Biquadratic trial function \n has the least error and almost coincides with the exact soluton")
plt.show()
