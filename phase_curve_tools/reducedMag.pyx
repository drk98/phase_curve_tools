cimport numpy as np
import numpy as np
cimport cython

ctypedef fused valid_arg:
    np.ndarray
    cython.float
    cython.double

cdef _calcReducedMag(valid_arg mag, valid_arg helioDist, valid_arg obsDist):
    return mag - 5*np.log10(helioDist*obsDist)

def calcReducedMag(valid_arg mag, valid_arg helioDist, valid_arg obsDist):
    r"""Calculate teh reduced magnitude
    ..math:

        V(\alpha) = m - 5\log10(\Delta r)
    
    :param mag: The magnitudes of the observations
    :type mag: np.ndarray
    :param helioDist: The heliocentric distances
    :type helioDist: np.ndarray
    :param obsDist: The observer/Earth distances
    :type obsDist: np.ndarray
    :return: The reduced magnitudes
    :rtype: np.ndarray
    """
    return _calcReducedMag(mag, helioDist, obsDist)