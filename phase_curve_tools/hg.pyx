cimport numpy as np
import numpy as np
cimport cython


from phase_curve_tools.valid_arg cimport valid_arg
    

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

cdef _calcPhis_simple(valid_arg phaseAngle, cython.float G):
    cdef valid_arg phi1, phi2, tan_val

    tan_val = np.tan(.5 * (phaseAngle * (np.pi/180)))
    phi1 = np.exp(-HG_SIMPLE_A1 * np.power(tan_val, HG_SIMPLE_B1))
    phi2 = np.exp(-HG_SIMPLE_A2 * np.power(tan_val, HG_SIMPLE_B2))
    return phi1, phi2

cdef _phaseIntegral_simple(valid_arg phaseAngle, cython.float G):
    cdef valid_arg phi1, phi2
    phi1, phi2 = _calcPhis_simple(phaseAngle, G)
    return (1-G) * phi1 + G * phi2

cdef _calcPhis(valid_arg phaseAngle, cython.float G):
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

    return phi1, phi2

cdef _phaseIntegral(valid_arg phaseAngle, cython.float G):
    
    cdef valid_arg phi1, phi2

    phi1, phi2 = _calcPhis(phaseAngle, G)
    return (1-G) * phi1 + G * phi2

cdef _HG(valid_arg apmag, valid_arg helioDist, valid_arg obsDist, valid_arg phaseAngle ,cython.float G, const cython.int simple):
    cdef valid_arg phaseInt
    if simple:
        phaseInt = _phaseIntegral_simple(phaseAngle, G)
    else:
        phaseInt = _phaseIntegral(phaseAngle, G)

    return apmag - (5 * np.log10(helioDist*obsDist) - 2.5*np.log10(phaseInt))


cdef _calcPhaseCurve(float H, float G, valid_arg phaseAngle):
    cdef float a1, a2
    cdef valid_arg phi1, phi2

    a1 = np.power(10, -H/2.5) / (1+(G/(1-G)))
    a2 = G*a1/(1-G)

    phi1, phi2 = _calcPhis(phaseAngle, G)
    return -2.5 * np.log10(a1*phi1 + a2*phi2)

def calcPhaseCurve(float H, float G, valid_arg phaseAngle):
    r"""Calculate the reduced magnitudes, the Bowell HG system\ :footcite:p:`bowellHG` for the given phase angles.

    Will use the non-simple :math:`\phi_1,\phi_2` calculations

    :param H: The Bowell H value
    :type H: float
    :param G: The Bowell G value 
    :type G: float
    :param phaseAngle: The phase angle, the sun - object - observer angle (degrees)
    :type phaseAngle: np.ndarray
    :return: The reduced magnitudes
    :rtype: valid_arg
    """

    return _calcPhaseCurve(H, G, phaseAngle)

def bowellCalcAbsMag(valid_arg apmag, valid_arg helioDist, valid_arg obsDist, valid_arg phaseAngle ,const cython.float G, const cython.int simple = True):
    r"""Calculate the absolute mags for an object from the Bowell HG system\ :footcite:p:`bowellHG`
    
    

    :param apmag: The apparent magnitudes to use
    :type apmag: np.ndarray
    :param helioDist: The distance from the object to the Sun (AU)
    :type helioDist: np.ndarray
    :param obsDist: The distance from the object to the Observer/Earth (AU)
    :type obsDist: np.ndarray
    :param phaseAngle: The phase angle, the sun - object - observer angle (degrees)
    :type phaseAngle: np.ndarray
    :param G: The Bowell G value
    :type G: float
    :param simple: If to use the simple or other version of :math:`\phi` calculations. The usual version used is the simple version. defaults to True
    :type simple: bool
    :return: The absolute magnitudes for the observations.
    :rtype: np.ndarray
    """
    return _HG(apmag, helioDist, obsDist, phaseAngle ,G, simple)

