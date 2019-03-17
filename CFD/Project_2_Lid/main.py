from sympy import *
import matplotlib.pyplot as plt
import numpy as np
from array import *
from scipy.sparse import *
from gridGen import grid,gridPlot

x,y=grid(10,10,3)   # accepts (nx, ny, game)
gridPlot(grid(10,10,3)[0],grid(10,10,3)[1])
