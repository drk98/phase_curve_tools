cimport numpy as np
import numpy as np
cimport cython

#cython: language_level=3

ctypedef fused valid_arg:
    np.ndarray
    cython.float
    cython.double

cdef _calcReducedMag(valid_arg mag, valid_arg helioDist, valid_arg obsDist):
    return mag - 5*np.log10(helioDist*obsDist)

def calcReducedMag(valid_arg mag, valid_arg helioDist, valid_arg obsDist):
    return _calcReducedMag(mag, helioDist, obsDist)