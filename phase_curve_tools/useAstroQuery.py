from ast import Call

import numpy as np
from astroquery.jplhorizons import HorizonsClass
from astropy import units as u
from typing import Callable
from .hg import bowellCalcAbsMag

class HorAbsMags(HorizonsClass):
    """Wrapper for astroquries Horizons class. can calculate the absolute magnitude.

    :param mags: The mags of observations at the provided ephocs
    :type mags: np.ndarray
    """

    def __init__(self, mags:np.ndarray, *args, **kwargs) -> None:
        """Constructor
        """
        self.mags:np.ndarray = mags
        eph = self.ephemerides()
        super().__init__(*args, **kwargs)

    def calcAbsMag(self, absMagCalcFunc:Callable = bowellCalcAbsMag, G:float=None) -> np.ndarray:
        """calcAbsMag Function to calculate the absolute magnitudes. Will grab the data from horizons

        :param absMagCalcFunc: The function that calcualtes the absolute magnitude, defaults to HG
        :type absMagCalcFunc: Callable, optional
        :param G: The G value to use. Use None if using the value from horizons, defaults to None
        :type G: float, optional
        :return: The absolute magnitudes for the observations.
        :rtype: np.ndarray
        """
        if G is None:
            G = self.eph["G"][0]
        helioDist = self.eph["r"]
        obsDist = self.eph["delta"]
        phaseAngle = self.eph["alpha"]

        return absMagCalcFunc(self.mags, helioDist, obsDist, phaseAngle, G)
    
    def getLightTimeCorrection(self, outputUnit:u = u.day)-> np.ndarray:
        """getLightTimeCorrection Calculate the light time corrction for the observations

        :param outputUnit: The output unit to use, defaults to u.day
        :type outputUnit: u, optional
        :return: The light time correction for each ephoch
        :rtype: np.ndarray
        """
        return self.eph["lighttime"].to(outputUnit).data