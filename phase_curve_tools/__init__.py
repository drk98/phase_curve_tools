"""Top-level package for Phase Curve Tools."""

__author__ = """Daniel Kramer"""
__email__ = 'drk98@nau.edu'
__version__ = '0.4.4'


from .hg import bowellCalcAbsMag, calcPhaseCurve
from .reducedMag import calcReducedMag
from .hg1g2 import calcHG1G2
from .hg12 import calcHG12

from .useAstroQuery import HorAbsMags
from .bestFit import BestFit, HGResult