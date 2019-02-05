import matplotlib.pyplot as plt
import numpy as np 

x = np.arange(1, 2, 0.01)
u_exact = (2/x)+0.5*np.log(x)

plt.plot(x, u_exact)
plt.xlabel('x', fontsize=16)
plt.ylabel('u', fontsize=16)
plt.show()