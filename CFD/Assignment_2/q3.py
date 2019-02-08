import numpy as np
import matplotlib.pyplot as plt
import math
# explicit scheme
u=[]
t=0.01
for i in range(200):
    if i==0:
        u.append(-1)
    else:
        u.append(u[i-1]+(0.01)*t*np.exp(-t))
    t=t+0.01
# implicit scheme
ub=[]
t=0.01
for i in range(200):
    if i==0:
        ub.append(-1+ (0.01)*(t+0.01)*np.exp(-(t+0.01)))
    else:
        ub.append(ub[i-1]+(0.01)*(t+0.01)*np.exp(-(t+0.01)))
    t=t+0.01
# CN scheme
uc=[]
t=0.01
for i in range(200):
    if i==0:
        uc.append(-1+(0.5)*(0.01)*(t*np.exp(-t)+(t+0.01)*np.exp(-(t+0.01))))
    else:
        uc.append(uc[i-1]+(0.5)*(0.01)*(t*np.exp(-t)+(t+0.01)*np.exp(-(t+0.01))))
    t=t+0.01
# CN scheme end
t = np.arange(0, 2, 0.01)
u_exact = -1*(t+1)*np.exp(-t)

plt.plot(t, u_exact, label="Exact")
plt.plot(t, u, label="explicit")
plt.plot(t, ub, label="implicit")
plt.plot(t, uc, label="CN")
plt.legend()
plt.xlabel('Time,t', fontsize=16)
plt.ylabel('Solution, u', fontsize=16)
plt.show()
