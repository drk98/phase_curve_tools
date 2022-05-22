#!/usr/bin/env python

"""Tests for `phase_curve_tools` package."""

import pytest


from phase_curve_tools import bowellCalcAbsMag, calcHG1G2, calcHG12



def test_hg_simple():
    assert bowellCalcAbsMag(18.168, 1.997, 1.187, 22.9769, .15) == pytest.approx(15.201,.05)

def test_hg():
    assert bowellCalcAbsMag(18.168, 1.997, 1.187, 22.9769, .15, simple=False) == pytest.approx(15.201,.05)


def test_hg1g2():

    assert calcHG1G2(18.168, 1.997, 1.187, 22.9769, .1, .2) == pytest.approx(15.201,.05)

def test_hg12():

    assert calcHG12(18.168, 1.997, 1.187, 22.9769, .15) == pytest.approx(15.201,.05)