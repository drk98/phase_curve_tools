#!/usr/bin/env python

"""Tests for `phase_curve_tools` package."""

import pytest


from phase_curve_tools import phase_curve_tools, hg



def test_hg_simple():
    assert hg.HG(18.168, 1.997, 1.187, 22.9769, .15) == pytest.approx(15.201,.05)

def test_hg():
    assert hg.HG(18.168, 1.997, 1.187, 22.9769, .15, simple=False) == pytest.approx(15.201,.05)