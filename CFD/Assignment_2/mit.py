import numpy as np
from numpy import *
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt
from matplotlib import animation


N = 100
T = 100
alpha = 2.0
u = zeros(N)
uu = zeros([T, N])
u[N/2] = 1.0
k = arange(N)
uu[0,:] = u


def tridiag(a, b, c, k1=-1, k2=0, k3=1):
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

## Tri Diagonal Matrix Algorithm(a.k.a Thomas algorithm) solver
def TDMAsolver(a, b, c, d):
    '''
    TDMA solver, a b c d can be NumPy array type or Python list type.
    refer to http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    '''
    nf = len(a)     # number of equations
    ac, bc, cc, dc = map(np.array, (a, b, c, d))     # copy the array
    for it in xrange(1, nf):
        mc = ac[it]/bc[it-1]
        bc[it] = bc[it] - mc*cc[it-1] 
        dc[it] = dc[it] - mc*dc[it-1]
 
    xc = ac
    xc[-1] = dc[-1]/bc[-1]
 
    for il in xrange(nf-2, -1, -1):
        xc[il] = (dc[il]-cc[il]*xc[il+1])/bc[il]
 
    del bc, cc, dc  # delete variables from memory
 
    return xc

def buildMatrix():
	global alpha, N
	a = -alpha*ones(N-1); a[N-2] = 0 #a[0] = 0;
	b = (1+2*alpha)*ones(N); b[0] = 1; b[N-1] = 1
	c = -alpha*ones(N-1); c[0] = 0 #c[N-2] = 0
	A = tridiag(a,b,c)
	return A


for n in range(T):
	a = -alpha*ones(N); a[0] = 0; a[N-1] = 0 #a[0] = 0;
	b = (1+2*alpha)*ones(N); b[0] = 1; b[N-1] = 1
	c = -alpha*ones(N); c[0] = 0; c[N-1] = 0
	d = u
	u = TDMAsolver(a,b,c,d)
	# print(u)
	uu[n,:] = u

# GAUSS ELIMINATION SOLVER
# A = buildMatrix()
# # print(A)
# # print(u)
# for n in range(T):
# 	A = buildMatrix()
# 	for i in range(1,N):
# 		a = A[i-1,i-1]; b = A[i,i-1]
# 		A[i,:] = np.add(A[i,:],np.multiply(A[i-1,:],-b/a))
# 		u[i] += u[i-1]*-b/a
# 		# print(-b/a)
# 		b = A[i,i]
# 		A[i,:] = np.divide(A[i,:],b)
# 		u[i] /= b
# 	# 	print(b)
# 	# 	print(A)
# 	# 	print(u)
# 	# print("Down")
# 	for i in reversed(range(1,N)):
# 		c = A[i-1,i]; b = A[i,i]
# 		A[i-1,:] = np.add(A[i-1,:],np.multiply(A[i,:],-c/b))
# 		u[i-1] += u[i]*-c/b
# 	# 	print(A)
# 	# 	print(u)
# 	# print("back up")
# 	uu[n,:] = u
		

# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure()
ax = fig.add_subplot(111,xlim=(0, N), ylim=(-0.1, 1))
ax.axis("off")
point, = ax.plot([], [],'b')

# initialization function: plot the background of each frame
def init_ani():
    point.set_data([], [])
    return point,
 
# animation function.  This is called sequentially
def animate(i):
	global uu, k
	x = uu[i,:]
	point.set_data(k,x)
	return point,

# # call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init_ani,
                               frames=T, interval=20, blit=True)
 
# this is how you save your animation to file:
anim.save('implicit_diffusion_2.gif', writer='imagemagick', fps=10)
 
plt.show()