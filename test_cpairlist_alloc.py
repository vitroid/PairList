#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for empty results, overflow rejection, and allocation-safe paths."""

import numpy as np
import pytest

from cpairlist import pairs, pairs2


def test_pairs_empty_input():
    rpos = np.zeros((0, 3), dtype=np.float64)
    result = pairs(rpos, 2, 2, 2)
    assert result.shape == (0, 2)
    assert result.dtype == np.dtype(np.intc) or result.dtype.kind == "i"


def test_pairs2_empty_inputs():
    r0 = np.zeros((0, 3), dtype=np.float64)
    r1 = np.zeros((0, 3), dtype=np.float64)
    result = pairs2(r0, r1, 2, 2, 2)
    assert result.shape == (0, 2)


def test_pairs_no_neighbors_still_ok():
    # single atom → zero pairs
    rpos = np.array([[0.1, 0.2, 0.3]], dtype=np.float64)
    result = pairs(rpos, 4, 4, 4)
    assert result.shape == (0, 2)


def test_pairs_rejects_huge_ngrid_product():
    rpos = np.zeros((2, 3), dtype=np.float64)
    # product exceeds 32-bit int
    with pytest.raises(ValueError, match="too large"):
        pairs(rpos, 100000, 100000, 100000)
