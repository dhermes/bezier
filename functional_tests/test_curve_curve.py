# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

import struct

import numpy as np
import six

import bezier
from bezier import _intersection_helpers


PACK_DOUBLE = struct.Struct('>d').pack
CURVE1 = bezier.Curve(np.array([
    [0.0, 0.0],
    [0.5, 1.0],
    [1.0, 0.0],
]))
CURVE2 = bezier.Curve(np.array([
    [1.125, 0.5],
    [0.625, -0.5],
    [0.125, 0.5],
]))


def to_bits(byte_):
    binary_rep = bin(ord(byte_))
    assert binary_rep[:2] == '0b'
    binary_rep = binary_rep[2:].zfill(8)
    assert len(binary_rep) == 8
    return binary_rep


def binary_representation(val):
    as_bytes = PACK_DOUBLE(val)
    as_bits = ''.join(map(to_bits, as_bytes))

    sign = as_bits[0]
    exponent = as_bits[1:12]
    mantissa = as_bits[12:]
    return sign, exponent, mantissa


def assert_close(approximated, exact):
    sign_a, exponent_a, mantissa_a = binary_representation(approximated)
    sign_e, exponent_e, mantissa_e = binary_representation(exact)
    assert sign_a == sign_e
    assert exponent_a == exponent_e

    mantissa_err = abs(int(mantissa_a, 2) - int(mantissa_e, 2))
    # Make sure the error is isolated to the last 3 bits.
    assert mantissa_err < 8


def curve_curve_check(curve1, curve2, s_vals, t_vals, points):
    assert len(s_vals) == len(t_vals)
    assert len(s_vals) == len(points)

    intersections = _intersection_helpers.all_intersections(
        [(curve1, curve2)])
    assert len(intersections) == len(s_vals)

    info = six.moves.zip(intersections, s_vals, t_vals, points)
    for intersection, s_val, t_val, point in info:
        assert intersection.left is curve1
        assert intersection.right is curve2

        assert_close(intersection._s_val, s_val)
        assert_close(intersection._t_val, t_val)

        assert_close(intersection.point[0], point[0])
        assert_close(intersection.point[1], point[1])


def test_curves1_and_2():
    sq31 = np.sqrt(31.0)
    s_val0 = 0.0625 * (9.0 - sq31)
    s_val1 = 0.0625 * (9.0 + sq31)

    s_vals = np.array([s_val0, s_val1])
    t_vals = np.array([s_val1, s_val0])
    points = np.array([
        [s_val0, (16.0 + sq31) / 64.0],
        [s_val1, (16.0 - sq31) / 64.0],
    ])
    curve_curve_check(CURVE1, CURVE2, s_vals, t_vals, points)
