from sympy.solvers import solve
from sympy import Symbol
from sympy import *
import matplotlib.pyplot as plt
import numpy as np

x=Symbol('x')
a1=Symbol('a1')
a2=Symbol('a2')
a3=Symbol('a3')
a4=Symbol('a4')
eqn_gal3 =[]
eqn_gal2=[]
eqn_gal1=[]
u_gal3=[]
u_gal2=[]
u_gal1=[]

# trial function
u2 = a1 + a2*x + a3*x*x
# imposing BCs
z2=u2.subs(a1,solve(u2.subs(x,1)-2 , a1)[0])
z2=z2.subs(a2,solve(((diff(u2,x)+2*u2-5).subs(x,2)).subs(a1,solve(u2.subs(x,1)-2 , a1)[0]),a2)[0])
# finding residual
r_gal2=x*x*diff(diff(z2,x),x)+2*x*diff(z2,x)+x-1
# galerkian integral wrt weighting functions
eqn_gal2.append(integrate(r_gal2*x, (x, 1, 2)))
z_gal2 = solve(eqn_gal2, a3)
print(z_gal2)
