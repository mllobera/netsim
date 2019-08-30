''' 
This module defines constants and coefficients used to calculate distance transforms
'''


# imports
import numpy as np

#: chamfer window size
wsize = 5

# neighborhood side size
nsize = wsize // 2

#: chamfer num coefficients
ncoef = 8

#: chamfer coefficient
a1, a2, a3 = 2.2062, 1.4141, 0.9866

# local distance metric
ldm = np.array([[a1,a1,a1,a2,a3,a2,a1,a3],
               [a3,a1,a2,a3,a2,a1,a1,a1]])

# chamfer window offsets
dx  = np.array([[-2,-2,-1,-1,-1,-1,-1, 0],
                [0, 1, 1, 1, 1, 1, 2, 2]])
dy  = np.array([[-1, 1,-2,-1, 0, 1, 2,-1],
                [1,-2,-1, 0, 1, 2,-1, 1]])

#: threshold for distance transform convergence
threshold = 0.5

# dtype
DTYPE = np.float64
