import numpy as np
import matplotlib.pyplot as plt
import math
# explicit scheme
u=[]
t=0.01
for i in range(200):
    if i==0: # case t=0
        u.append(1-0.01)
    else:
        u.append(u[i-1]+(0.01)*(2*t*np.exp(-t)-u[i-1]))
    t=t+0.01
# implicit scheme
ub=[]
t=0.01
for i in range(200):
    if i==0: # case t=0
        ub.append((1+ (0.01)*2*(0.01)*np.exp(-0.01))/(1+0.01))
    else:
        ub.append((ub[i-1]+(0.01)*2*(t+0.01)*np.exp(-(t+0.01)))/(1+0.01))
    t=t+0.01
# CN scheme
uc=[]
t=0.01
for i in range(200):
    if i==0:
        uc.append((1+(0.5)*(0.01)*(2*t*np.exp(-t)-1+2*(t+0.01)*np.exp(-(t+0.01))))/(1+0.005))
    else:
        uc.append((uc[i-1]+(0.5)*(0.01)*(2*t*np.exp(-t)-uc[i-1]+2*(t+0.01)*np.exp(-(t+0.01))))/(1+0.005))
    t=t+0.01
# CN scheme end
t = np.arange(0, 2, 0.01)

plt.plot(t, u, label="explicit")
plt.plot(t, ub, label="implicit")
plt.plot(t, uc, label="CN")
plt.legend()
plt.xlabel('t (s)', fontsize=16)
plt.ylabel('Solution, u', fontsize=16)
plt.suptitle('u vs time')
plt.show()
