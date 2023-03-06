"""
@author: Mad Dreamer
"""
import numpy as np
from functions import *

# quadrilateral
#x = np.array([[0,0],
#              [1,0],
#              [1,2],
#              [0,2]])
#GPE = 4 # (gauss per element ?!)

# triangle
x = np.array([[0,0],
              [2,0],
              [1,1]])
GPE = 1 # (gauss per element ?!)

K = stiffness(x, GPE)
print(K)