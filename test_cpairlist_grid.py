#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for ngrid lower bound and fractional-coordinate grid clamping."""

import numpy as np
import pytest

from cpairlist import pairs, pairs2


def test_pairs_rejects_ngrid_zero():
    rpos = np.zeros((2, 3), dtype=np.float64)
    with pytest.raises(ValueError, match=r"ngrid\[0\] must be >= 1"):
        pairs(rpos, 0, 2, 2)


def test_pairs_rejects_ngrid_negative():
    rpos = np.zeros((2, 3), dtype=np.float64)
    with pytest.raises(ValueError, match=r"ngrid\[2\] must be >= 1"):
        pairs(rpos, 2, 2, -1)


def test_pairs_accepts_ngrid_one():
    rpos = np.array([[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]], dtype=np.float64)
    result = pairs(rpos, 1, 1, 1)
    # single cell → one unordered pair
    assert result.shape == (1, 2)


def test_pairs2_ngrid_one_no_duplicate():
    rpos0 = np.array([[0.1, 0.1, 0.1]], dtype=np.float64)
    rpos1 = np.array([[0.5, 0.5, 0.5]], dtype=np.float64)
    result = pairs2(rpos0, rpos1, 1, 1, 1)
    # one atom each in the only cell → exactly one pair (not 27×)
    assert result.shape == (1, 2)


def test_pairs_clamps_product_rounding_to_ngrid():
    # x in [0,1) but x*ngrid may round to ngrid → must clamp, not OOB
    x = np.nextafter(1.0, 0.0)
    rpos = np.array([[x, x, x], [0.0, 0.0, 0.0]], dtype=np.float64)
    result = pairs(rpos, 8, 8, 8)
    assert result.ndim == 2
    assert result.shape[1] == 2
