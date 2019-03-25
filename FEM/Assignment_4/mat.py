from numpy import *
import scipy.linalg
# https://docs.scipy.org/doc/numpy/user/numpy-for-matlab-users.html
# http://mathesaurus.sourceforge.net/matlab-numpy.html

a =array([[1.,2.,3.], [4.,5.,6.]])
print(a[0,2]*0.4)
a[:]=55
print(a)
