from sympy import *
from sympy.solvers import solve
from sympy import Symbol
import matplotlib.pyplot as plt
import numpy as np

x=Symbol('x')
a1=Symbol('a1')
a2=Symbol('a2')
a3=Symbol('a3')
wf1=Symbol('wf1')
wf2=Symbol('wf2')
eqn=[]

u=2+(x-1)*a1+(x*x-1)*a2+a3*(x**3-1)
up=diff(u,x)
n=[x-1,x*x-1,x**3-1]
print(n)
for i in range(3):
    print(n[i])
    w=n[i]*(x**2)*(5-2*u)-integrate(diff(n[i],x)*(x**2)*up,x)+integrate((x-1)*n[i],x)
    wf1=w.subs(x,1)
    wf2=w.subs(x,2)
    eqn.append(wf2-wf1)
print(w)
print(eqn)
sol=solve(eqn,[a1, a2, a3])
print(sol)

u_ext=[]
u_gal=[]
y = []
i=1
while i<=2:
    y.append(i)
    u_ext.append(((0.1)*(29.9096-(5*x)-(4.9096/x))+log(x)).subs(x,i))
    u_gal.append((((u.subs(a1,sol.get(a1))).subs(a2,sol.get(a2))).subs(a3,sol.get(a3))).subs(x,i))
    i+=0.01
plt.plot(y, u_ext, label="Exact")
plt.plot(y, u_gal, label="Galerkin")
plt.legend()
plt.xlabel('x', fontsize=16)
plt.ylabel('u', fontsize=16)
plt.suptitle('Exact soluton vs Galerkin weak form')
plt.show()
