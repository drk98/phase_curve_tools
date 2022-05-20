from dataclasses import dataclass
from scipy.optimize import minimize
import numpy as np
from typing import Tuple, Union
from astroquery.jplhorizons import Horizons
from .hg import calcPhaseCurve

@dataclass(frozen=True)
class HGResult:
    """Class for hoding the reult of HG minimumzation 

    :param H: The Bowell H value
    :type H: float
    :param G: The Bowell G value 
    :type G: float
    :param sigH: The error in the Bowell H value
    :type sigH: float
    :param sigG: The error in the Bowell G value 
    :type sigG: float
    :param success: If the minimization succedded (according to scipy)
    :type success: bool
    :param atEdge: If the derived G value is at the edge of the porvided G bounds 
    :type atEdge: bool
    """
    H:float
    G:float
    sigH:float = 0
    sigG:float = 0
    success:bool = False
    atEdge:bool = False

    @property
    def fullSuccess(self)->bool:
        """fullSuccess If the minimzation succeded and the G value is not at an edge

        :getter: Returns if there was a full succedd
        :type: bool
        """
        return self.success and self.atEdge

class BestFit:
    """ Class to find the best fit in the Bowell HG system.

    If the helioDist nor obsDist is provided, then Horizons will be queried, so the object ID and times are needed. They don't have to be provided if helioDist and obsDist are provided.
    
    If no phase angles are provided, then they will be computed, assuming the distance from the observer/Earth to the Sun is always 1 AU.

    :param mags: The magnitudes of the observations
    :type mags: np.ndarray
    :param helioDist: The distance from the Sun to the object (AU), defaults to None
    :type helioDist: np.ndarray, optional
    :param obsDist: The distance from the Observer/Earth to the object (AU), defaults to None
    :type obsDist: np.ndarray, optional
    :param phaseAngle: The phase angle of the observations (deg), defaults to None
    :type phaseAngle: np.ndarray, optional
    :param times: The times of the observations (JD), defaults to None
    :type times: np.ndarray, optional
    :param objectID: The object ID, needs to be compatable with JPL Horizons, defaults to None
    :type objectID: Union[str, int], optional
    :raises ValueError: Needs the object ID and the times of the observations if helioDist or obsDist is not provided
    """

    def __init__(self, mags:np.ndarray, helioDist:np.ndarray=None, obsDist:np.ndarray=None, phaseAngle:np.ndarray=None, times:np.ndarray=None, objectID:Union[str, int]=None) -> None:
        """Constructor
        """
        
        self.helioDist:np.ndarray = helioDist
        self.obsDist:np.ndarray = obsDist
        self.phaseAngle:np.ndarray = phaseAngle

        if helioDist is None or obsDist is None:

            # Need an object ID and times to query horizions
            if objectID is None:
                raise ValueError("No heliocentric distances orobserver-centric distances provided and no object id provided. Unable to procede.")
            if times is None:
                raise ValueError("No heliocentric distances orobserver-centric distances provided and no times provided. Unable to procede.")

            eph = Horizons(id=objectID, epochs=times).ephemerides()
            self.helioDist = self.eph["r"].data
            self.obsDist = self.eph["delta"].data
            self.phaseAngle = self.eph["alpha"].data
        elif phaseAngle is None:
            self.phaseAngle = BestFit.calcPhaseAngle(self.helioDist, self.obsDist)
        
        self.mags:np.ndarray = mags

        self.reducedMagHG:np.ndarray = None

    @staticmethod
    def calcPhaseAngle(helioDist:np.ndarray, obsDist:np.ndarray, sunObsDist:Union[np.ndarray, float]=1)->np.ndarray:
        r"""calcPhaseAngle Calculate the phase angles.
        To calculate the s-t-o phase angle, using the law of cosines:
        
        .. math:: \alpha=\cos^{-1}\left(\dfrac{\Delta^2+r^2-d^2}{2\Delta{r}}\right) 
        
        where :math:`\Delta` is the heliocentric distance, :math:`r` is the observer/Earth distance, and :math:`d` is the distance from the observer/Earth to the Sun.

        :param helioDist: The heliocentric distance (AU)
        :type helioDist: np.ndarray
        :param obsDist: The observer distance (AU)
        :type obsDist: np.ndarray
        :param sunObsDist: The distance from the observer/Earth to the Sun (AU), defaults to 1
        :type sunObsDist: Union[np.ndarray, float], optional
        :return: The phase angles (deg)
        :rtype: np.ndarray
        """
        return np.arccos((helioDist**2 + obsDist**2 + sunObsDist**2)/(2*helioDist*obsDist))

            
    def fitHG(self, H_0:float=None, G_0:float=.15, HRange:Tuple=(4.69236933828342, 29.276432904419803), GRange:Tuple=(-0.301744075798055, 0.9073912874142601))->HGResult:
        r"""fitHG Fit the data to the Bowell HG system\ :footcite:p:`bowellHG`

        :param H_0: The initial H guess, if None, then will use the reduced magnitude at the minimum phase angle, defaults to None
        :type H_0: float, optional
        :param G_0: The inital G guess, defaults to .15
        :type G_0: float, optional
        :param HRange: The valid range for H values, defaults to (4.69236933828342, 29.276432904419803), which are :math:`\pm10\%` the values from \ :footcite:t:`verdes`
        :type HRange: Tuple, optional
        :param GRange: The valid range for G values, defaults to (-0.301744075798055, 0.9073912874142601), which are :math:`\pm10\%` the values from \ :footcite:t:`verdes`
        :type GRange: Tuple, optional
        :return: The result of the fit
        :rtype: HGResult
        """
        
        if self.reducedMag is None:
            self.reducedMagHG = self.mags - 5*np.log10(self.helioDist*self.obsDist)

        if H_0 is None:
            H_0 = self.reducedMag[np.argmin(self.phaseAngle)]
        
        res = minimize(self._HGminimize_target, x0=[H_0,G_0], bounds=(HRange, GRange))

        return HGResult(res.x[0], res.x[1], 0, 0, res.success, res.x[1] in GRange)
        
    def _HGminimize_target(self, _hg:Tuple[float,float])->float:
        
        adj_mag = calcPhaseCurve(_hg[0], _hg[1], self.phaseAngle)
        return np.sqrt(np.sum(np.power(self.reducedMagHG - adj_mag, 2)) / (len(self.reducedMagHG) - 1))
