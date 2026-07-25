#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Tests for cpairlist input array validation (contiguous / shape / dtype)."""

import numpy as np
import pytest

from cpairlist import pairs, pairs2


def _sample_rpos(n=4, dtype=np.float64, order="C"):
    base = np.array(
        [
            [0.1, 0.1, 0.1],
            [0.1, 0.1, 0.6],
            [0.6, 0.1, 0.1],
            [0.6, 0.6, 0.6],
        ],
        dtype=np.float64,
    )[:n]
    return np.array(base, dtype=dtype, order=order)


def test_pairs_accepts_c_contiguous_float64():
    rpos = _sample_rpos()
    assert rpos.flags.c_contiguous
    result = pairs(rpos, 2, 2, 2)
    assert result.ndim == 2
    assert result.shape[1] == 2


def test_pairs_accepts_float32_via_cast():
    rpos = _sample_rpos(dtype=np.float32)
    result = pairs(rpos, 2, 2, 2)
    assert result.ndim == 2


def test_pairs_accepts_fortran_order_via_copy():
    rpos = _sample_rpos(order="F")
    assert not rpos.flags.c_contiguous
    result = pairs(rpos, 2, 2, 2)
    assert result.ndim == 2


def test_pairs_accepts_noncontiguous_view_via_copy():
    wide = np.zeros((4, 6), dtype=np.float64)
    wide[:, 0:3] = _sample_rpos()
    rpos = wide[:, 0:3]  # non-contiguous view
    assert not rpos.flags.c_contiguous
    result = pairs(rpos, 2, 2, 2)
    assert result.ndim == 2


def test_pairs_rejects_wrong_column_count():
    rpos = np.zeros((4, 2), dtype=np.float64)
    with pytest.raises(ValueError, match=r"shape \(N, 3\)"):
        pairs(rpos, 2, 2, 2)


def test_pairs_rejects_1d():
    rpos = np.zeros(12, dtype=np.float64)
    with pytest.raises(ValueError, match="2-dimensional"):
        pairs(rpos, 2, 2, 2)


def test_pairs2_rejects_wrong_shape_on_second_arg():
    rpos0 = _sample_rpos()
    rpos1 = np.zeros((3, 4), dtype=np.float64)
    with pytest.raises(ValueError, match="rpos1"):
        pairs2(rpos0, rpos1, 2, 2, 2)
