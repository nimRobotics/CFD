import matplotlib.pyplot as plt
import numpy as np 

x = np.arange(0.1, 5., 0.01)
u = (2/x)+0.5*np.log(x)-2

plt.plot(u, x)
plt.xlabel('u', fontsize=16)
plt.ylabel('x', fontsize=16)
plt.show()