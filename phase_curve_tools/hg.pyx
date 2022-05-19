cimport numpy as np
import numpy as np
cimport cython

#cython: language_level=3

ctypedef fused valid_arg:
    np.ndarray
    cython.float
    cython.double
    

cdef extern from "constants.h":
    cdef cython.float HG_SIMPLE_A1
    cdef cython.float HG_SIMPLE_A2 
    cdef cython.float HG_SIMPLE_B1 
    cdef cython.float HG_SIMPLE_B2 
    cdef cython.float HG_A1
    cdef cython.float HG_A2
    cdef cython.float HG_B1
    cdef cython.float HG_B2
    cdef cython.float HG_C1
    cdef cython.float HG_C2

cdef _phaseIntegral_simple(valid_arg phaseAngle, cython.float G):
    cdef valid_arg phi1, phi2, tan_val

    tan_val = np.tan(.5 * (phaseAngle * (np.pi/180)))
    phi1 = np.exp(-HG_SIMPLE_A1 * np.power(tan_val, HG_SIMPLE_B1))
    phi2 = np.exp(-HG_SIMPLE_A2 * np.power(tan_val, HG_SIMPLE_B2))

    return (1-G) * phi1 + G * phi2

cdef _phaseIntegral(valid_arg phaseAngle, cython.float G):
    cdef valid_arg W, phi1, phi2, phiS, phiL, sin_pa, tan_val

    sin_pa = np.sin(phaseAngle* (np.pi/180))
    tan_val = np.tan(.5 * (phaseAngle * (np.pi/180)))

    W = np.exp(-90.56 * tan_val*tan_val)

    phiS = 1 - ((HG_C1*sin_pa)/((.119+1.341*sin_pa) - (.754*sin_pa*sin_pa)))
    phiL = np.exp(-HG_A1 * np.power(tan_val, HG_B1))
    phi1 = W*phiS + (1-W)*phiL

    phiS = 1 - ((HG_C2*sin_pa)/((.119+1.341*sin_pa) - (.754*sin_pa*sin_pa)))
    phiL = np.exp(-HG_A2 * np.power(tan_val, HG_B2))
    phi2 = W*phiS + (1-W)*phiL

    return (1-G) * phi1 + G * phi2

cdef _HG(valid_arg apmag, valid_arg helioDist, valid_arg obsDist, valid_arg phaseAngle ,cython.float G, const cython.int simple):
    cdef valid_arg phaseInt
    if simple:
        phaseInt = _phaseIntegral_simple(phaseAngle, G)
    else:
        phaseInt = _phaseIntegral(phaseAngle, G)

    return apmag - (5 * np.log10(helioDist*obsDist) - 2.5*np.log10(phaseInt))

def HG(valid_arg apmag, valid_arg helioDist, valid_arg obsDist, valid_arg phaseAngle ,const cython.float G, const cython.int simple = True):
    """Calculate the absolute mags for an object from the Bowell HG system 

    :param apmag: The apparent magnitudes to use
    :type apmag: np.ndarray
    :param helioDist: The distance from the object to the Sun (AU)
    :type helioDist: np.ndarray
    :param obsDist: The distance from the object to the Observer/Earth (AU)
    :type obsDist: np.ndarray
    :param phaseAngle: The phase angle, the sun - object - observer angle (degrees)
    :type phaseAngle: np.ndarray
    :param G: The Bowell G value
    :type G: np.ndarray
    :return: The absolute magnitudes for the observations.
    :rtype: np.ndarray
    """
    return _HG(apmag, helioDist, obsDist, phaseAngle ,G, simple)
