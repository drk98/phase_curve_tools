cimport numpy as np
import numpy as np
cimport cython

ctypedef fused valid_arg:
    np.ndarray
    cython.float
    cython.double