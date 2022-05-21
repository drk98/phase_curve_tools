"""Top-level package for Phase Curve Tools."""

__author__ = """Daniel Kramer"""
__email__ = 'drk98@nau.edu'
__version__ = '0.2.1'


from .hg import bowellCalcAbsMag, calcPhaseCurve
from .reducedMag import calcReducedMag
from .useAstroQuery import HorAbsMags
from .bestFit import BestFit