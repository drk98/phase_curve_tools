cimport numpy as np
import numpy as np
cimport cython

from phase_curve_tools.valid_arg cimport valid_arg


cdef extern from "constants.h":
    cdef cython.float M_1_PI 
    cdef cython.float M_1_5PI 


cdef _calcPhis( valid_arg phaseAngle):
    cdef valid_arg phi1, phi2, phi3, pa_rad
    pa_rad = phaseAngle * np.pi/180
    phi1 = 1- (6 * M_1_PI * pa_rad)
    phi2 = 1- (9 * M_1_5PI * pa_rad)

    phi3 = np.exp(-4*np.pi*np.power(np.tan(.5 * pa_rad), 2/3))

    return phi1, phi2, phi3

cdef _calcPhaseInt( valid_arg phaseAngle,  float G1,  float G2):
    cdef valid_arg phi1, phi2, phi3

    phi1, phi2, phi3 = _calcPhis(phaseAngle)

    return G1*phi1 + G2*phi2 + (1-G1-G2)*phi3

    
cdef _HG1G2( valid_arg apmag,  valid_arg helioDist,  valid_arg obsDist,  valid_arg phaseAngle , cython.float G1,  cython.float G2):
    cdef valid_arg phaseInt

    phaseInt = _calcPhaseInt(phaseAngle, G1, G2)

    return apmag - (5 * np.log10(helioDist*obsDist) - 2.5*np.log10(phaseInt))

def calcHG1G2( valid_arg apmag,  valid_arg helioDist,  valid_arg obsDist,  valid_arg phaseAngle , cython.float G1,  cython.float G2):
    r"""Calculate the absolute mags for an object from the HG1G2 system\ :footcite:p:`MUINONEN2010542`
    
    

    :param apmag: The apparent magnitudes to use
    :type apmag: np.ndarray
    :param helioDist: The distance from the object to the Sun (AU)
    :type helioDist: np.ndarray
    :param obsDist: The distance from the object to the Observer/Earth (AU)
    :type obsDist: np.ndarray
    :param phaseAngle: The phase angle, the sun - object - observer angle (degrees)
    :type phaseAngle: np.ndarray
    :param G1: The G1 value
    :type G1: float
    :param G2: The G2 value
    :type G2: float
    :return: The absolute magnitudes for the observations.
    :rtype: np.ndarray
    """

    return _HG1G2(apmag, helioDist, obsDist, phaseAngle, G1, G2)