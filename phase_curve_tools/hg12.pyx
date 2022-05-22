cimport numpy as np
import numpy as np
cimport cython

cimport phase_curve_tools.hg1g2 as hg1g2
import phase_curve_tools.hg1g2 as hg1g2

from phase_curve_tools.valid_arg cimport valid_arg


cdef extern from "constants.h":
    cdef cython.float G1_LT_SCAL 
    cdef cython.float G1_LT_CONST
    cdef cython.float G1_GTE_SCAL
    cdef cython.float G1_GTE_CONST

    cdef cython.float G2_LT_SCAL
    cdef cython.float G2_LT_CONST
    cdef cython.float G2_GTE_SCAL
    cdef cython.float G2_GTE_CONST


cdef _HG12(valid_arg apmag,  valid_arg helioDist,  valid_arg obsDist,  valid_arg phaseAngle , cython.float G12):
    cdef valid_arg phaseInt
    cdef float G1, G2

    if G12 >= .2:
        G1 = G1_GTE_SCAL*G12 + G1_GTE_CONST
        G2 = G2_GTE_SCAL*G12 + G2_GTE_CONST
    else:
        G1 = G1_LT_SCAL*G12 + G1_LT_CONST
        G2 = G2_LT_SCAL*G12 + G2_LT_CONST

    phaseInt = hg1g2._calcPhaseInt(phaseAngle, G1, G2)

    return apmag - (5 * np.log10(helioDist*obsDist) - 2.5*np.log10(phaseInt))

def calcHG12(valid_arg apmag,  valid_arg helioDist,  valid_arg obsDist,  valid_arg phaseAngle , cython.float G12):
    r"""Calculate the absolute mags for an object from the HG12 system\ :footcite:p:`MUINONEN2010542`
    
    

    :param apmag: The apparent magnitudes to use
    :type apmag: np.ndarray
    :param helioDist: The distance from the object to the Sun (AU)
    :type helioDist: np.ndarray
    :param obsDist: The distance from the object to the Observer/Earth (AU)
    :type obsDist: np.ndarray
    :param phaseAngle: The phase angle, the sun - object - observer angle (degrees)
    :type phaseAngle: np.ndarray
    :param G12: The G12 value
    :type G12: float
    :return: The absolute magnitudes for the observations.
    :rtype: np.ndarray
    """

    return _HG12(apmag, helioDist, obsDist, phaseAngle, G12)