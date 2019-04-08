from numpy import *
import scipy.linalg
from numpy import linalg as LA
import numpy as np

wMat = ones((2,2))
wMatOld = 2*ones((2,2))
wMatOld[1,1]=3
rs1 =np.amax(abs(-wMatOld+wMat))
print("rs1 ",rs1)
