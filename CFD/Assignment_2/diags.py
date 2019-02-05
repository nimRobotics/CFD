from scipy import *
from scipy.sparse import *

diagonals = [[1,2,3,4], [1,2,3], [1,2]]
print(diags(diagonals, [0, -1, 2]).todense())
