# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Pure Python helper methods for :mod:`bezier.triangle`.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import functools
import operator

import numpy as np

from bezier.hazmat import curve_helpers
from bezier.hazmat import helpers as _py_helpers
from bezier.hazmat import intersection_helpers


_MAX_POLY_SUBDIVISIONS = 5
_SIGN = np.sign  # pylint: disable=no-member
_FLOAT64 = np.float64  # pylint: disable=no-member
_SAME_CURVATURE = "Tangent curves have same curvature."
_WRONG_CURVE = "Start and end node not defined on same curve"
CLASSIFICATION_T = intersection_helpers.IntersectionClassification
# NOTE: The ``SUBDIVIDE`` matrices are public since used in
#       the ``triangle`` module.
LINEAR_SUBDIVIDE_A = (
    np.asfortranarray([[2, 1, 1], [0, 1, 0], [0, 0, 1]], dtype=_FLOAT64) / 2.0
)
LINEAR_SUBDIVIDE_B = (
    np.asfortranarray([[0, 1, 1], [1, 0, 1], [1, 1, 0]], dtype=_FLOAT64) / 2.0
)
LINEAR_SUBDIVIDE_C = (
    np.asfortranarray([[1, 0, 0], [1, 2, 1], [0, 0, 1]], dtype=_FLOAT64) / 2.0
)
LINEAR_SUBDIVIDE_D = (
    np.asfortranarray([[1, 0, 0], [0, 1, 0], [1, 1, 2]], dtype=_FLOAT64) / 2.0
)
QUADRATIC_SUBDIVIDE_A = (
    np.asfortranarray(
        [
            [4, 2, 1, 2, 1, 1],
            [0, 2, 2, 0, 1, 0],
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 2, 1, 2],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ],
        dtype=_FLOAT64,
    )
    / 4.0
)
QUADRATIC_SUBDIVIDE_B = (
    np.asfortranarray(
        [
            [0, 0, 1, 0, 1, 1],
            [0, 1, 0, 1, 1, 2],
            [1, 0, 0, 1, 0, 1],
            [0, 1, 2, 1, 1, 0],
            [2, 1, 0, 1, 1, 0],
            [1, 1, 1, 0, 0, 0],
        ],
        dtype=_FLOAT64,
    )
    / 4.0
)
QUADRATIC_SUBDIVIDE_C = (
    np.asfortranarray(
        [
            [1, 0, 0, 0, 0, 0],
            [2, 2, 0, 1, 0, 0],
            [1, 2, 4, 1, 2, 1],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 1, 2, 2],
            [0, 0, 0, 0, 0, 1],
        ],
        dtype=_FLOAT64,
    )
    / 4.0
)
QUADRATIC_SUBDIVIDE_D = (
    np.asfortranarray(
        [
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0],
            [2, 1, 0, 2, 0, 0],
            [0, 1, 2, 0, 2, 0],
            [1, 1, 1, 2, 2, 4],
        ],
        dtype=_FLOAT64,
    )
    / 4.0
)
CUBIC_SUBDIVIDE_A = (
    np.asfortranarray(
        [
            [8, 4, 2, 1, 4, 2, 1, 2, 1, 1],
            [0, 4, 4, 3, 0, 2, 2, 0, 1, 0],
            [0, 0, 2, 3, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 4, 2, 1, 4, 2, 3],
            [0, 0, 0, 0, 0, 2, 2, 0, 2, 0],
            [0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 2, 1, 3],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        ],
        dtype=_FLOAT64,
    )
    / 8.0
)
CUBIC_SUBDIVIDE_B = (
    np.asfortranarray(
        [
            [0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
            [0, 0, 1, 0, 0, 1, 1, 1, 2, 3],
            [0, 1, 0, 0, 1, 1, 0, 2, 1, 3],
            [1, 0, 0, 0, 1, 0, 0, 1, 0, 1],
            [0, 0, 1, 3, 0, 1, 2, 1, 1, 0],
            [0, 2, 2, 0, 2, 2, 2, 2, 2, 0],
            [3, 1, 0, 0, 2, 1, 0, 1, 1, 0],
            [0, 1, 2, 3, 1, 1, 1, 0, 0, 0],
            [3, 2, 1, 0, 1, 1, 1, 0, 0, 0],
            [1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
        ],
        dtype=_FLOAT64,
    )
    / 8.0
)
CUBIC_SUBDIVIDE_C = (
    np.asfortranarray(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [3, 2, 0, 0, 1, 0, 0, 0, 0, 0],
            [3, 4, 4, 0, 2, 2, 0, 1, 0, 0],
            [1, 2, 4, 8, 1, 2, 4, 1, 2, 1],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 2, 2, 0, 2, 0, 0],
            [0, 0, 0, 0, 1, 2, 4, 2, 4, 3],
            [0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 1, 2, 3],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        ],
        dtype=_FLOAT64,
    )
    / 8.0
)
CUBIC_SUBDIVIDE_D = (
    np.asfortranarray(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [3, 1, 0, 0, 2, 0, 0, 0, 0, 0],
            [0, 2, 2, 0, 0, 2, 0, 0, 0, 0],
            [0, 0, 1, 3, 0, 0, 2, 0, 0, 0],
            [3, 2, 1, 0, 4, 2, 0, 4, 0, 0],
            [0, 1, 2, 3, 0, 2, 4, 0, 4, 0],
            [1, 1, 1, 1, 2, 2, 2, 4, 4, 8],
        ],
        dtype=_FLOAT64,
    )
    / 8.0
)
QUARTIC_SUBDIVIDE_A = (
    np.asfortranarray(
        [
            [16, 8, 4, 2, 1, 8, 4, 2, 1, 4, 2, 1, 2, 1, 1],
            [0, 8, 8, 6, 4, 0, 4, 4, 3, 0, 2, 2, 0, 1, 0],
            [0, 0, 4, 6, 6, 0, 0, 2, 3, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 2, 4, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 8, 4, 2, 1, 8, 4, 2, 6, 3, 4],
            [0, 0, 0, 0, 0, 0, 4, 4, 3, 0, 4, 4, 0, 3, 0],
            [0, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 2, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 2, 1, 6, 3, 6],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 3, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 4],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        ],
        dtype=_FLOAT64,
    )
    / 16.0
)
QUARTIC_SUBDIVIDE_B = (
    np.asfortranarray(
        [
            [0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1],
            [0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 2, 1, 3, 4],
            [0, 0, 1, 0, 0, 0, 1, 1, 0, 1, 2, 1, 3, 3, 6],
            [0, 1, 0, 0, 0, 1, 1, 0, 0, 2, 1, 0, 3, 1, 4],
            [1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1],
            [0, 0, 0, 1, 4, 0, 0, 1, 3, 0, 1, 2, 1, 1, 0],
            [0, 0, 2, 3, 0, 0, 2, 3, 3, 2, 3, 4, 3, 3, 0],
            [0, 3, 2, 0, 0, 3, 3, 2, 0, 4, 3, 2, 3, 3, 0],
            [4, 1, 0, 0, 0, 3, 1, 0, 0, 2, 1, 0, 1, 1, 0],
            [0, 0, 1, 3, 6, 0, 1, 2, 3, 1, 1, 1, 0, 0, 0],
            [0, 3, 4, 3, 0, 3, 3, 3, 3, 2, 2, 2, 0, 0, 0],
            [6, 3, 1, 0, 0, 3, 2, 1, 0, 1, 1, 1, 0, 0, 0],
            [0, 1, 2, 3, 4, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
            [4, 3, 2, 1, 0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
            [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        ],
        dtype=_FLOAT64,
    )
    / 16.0
)
QUARTIC_SUBDIVIDE_C = (
    np.asfortranarray(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [4, 2, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [6, 6, 4, 0, 0, 3, 2, 0, 0, 1, 0, 0, 0, 0, 0],
            [4, 6, 8, 8, 0, 3, 4, 4, 0, 2, 2, 0, 1, 0, 0],
            [1, 2, 4, 8, 16, 1, 2, 4, 8, 1, 2, 4, 1, 2, 1],
            [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 3, 2, 0, 0, 2, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 3, 4, 4, 0, 4, 4, 0, 3, 0, 0],
            [0, 0, 0, 0, 0, 1, 2, 4, 8, 2, 4, 8, 3, 6, 4],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 2, 0, 3, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4, 3, 6, 6],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 4],
            [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
        ],
        dtype=_FLOAT64,
    )
    / 16.0
)
QUARTIC_SUBDIVIDE_D = (
    np.asfortranarray(
        [
            [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [4, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 3, 2, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 2, 3, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 4, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0],
            [6, 3, 1, 0, 0, 6, 2, 0, 0, 4, 0, 0, 0, 0, 0],
            [0, 3, 4, 3, 0, 0, 4, 4, 0, 0, 4, 0, 0, 0, 0],
            [0, 0, 1, 3, 6, 0, 0, 2, 6, 0, 0, 4, 0, 0, 0],
            [4, 3, 2, 1, 0, 6, 4, 2, 0, 8, 4, 0, 8, 0, 0],
            [0, 1, 2, 3, 4, 0, 2, 4, 6, 0, 4, 8, 0, 8, 0],
            [1, 1, 1, 1, 1, 2, 2, 2, 2, 4, 4, 4, 8, 8, 16],
        ],
        dtype=_FLOAT64,
    )
    / 16.0
)
_WEIGHTS_SUBDIVIDE0 = np.asfortranarray([1.0, 0.0, 0.0])
_WEIGHTS_SUBDIVIDE1 = np.asfortranarray([0.5, 0.5, 0.0])
_WEIGHTS_SUBDIVIDE2 = np.asfortranarray([0.5, 0.0, 0.5])
_WEIGHTS_SUBDIVIDE3 = np.asfortranarray([0.0, 0.5, 0.5])
_WEIGHTS_SUBDIVIDE4 = np.asfortranarray([0.0, 1.0, 0.0])
_WEIGHTS_SUBDIVIDE5 = np.asfortranarray([0.0, 0.0, 1.0])
# The Jacobian of a quadratric (in any dimension) as given by
# dB/ds = [-2L1, 2(L1 - L2), 2L2, -2L3, 2L3, 0] * nodes
# dB/dt = [-2L1, -2L2, 0, 2(L1 - L3), 2L2, 2L3] * nodes
# We evaluate this at each of the 6 points in the quadratic
# triangle and then stack them (2 columns * 6 = 12 columns)
_QUADRATIC_JACOBIAN_HELPER = np.asfortranarray(
    [
        [-2, -2, -1, -1, 0, 0, -1, -1, 0, 0, 0, 0],
        [2, 0, 0, -1, -2, -2, 1, 0, -1, -1, 0, 0],
        [0, 0, 1, 0, 2, 0, 0, 0, 1, 0, 0, 0],
        [0, 2, 0, 1, 0, 0, -1, 0, -1, -1, -2, -2],
        [0, 0, 0, 1, 0, 2, 1, 0, 1, 1, 2, 0],
        [0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2],
    ],
    dtype=_FLOAT64,
)
_QUADRATIC_TO_BERNSTEIN = (
    np.asfortranarray(
        [
            [2, -1, 0, -1, 0, 0],
            [0, 4, 0, 0, 0, 0],
            [0, -1, 2, 0, -1, 0],
            [0, 0, 0, 4, 0, 0],
            [0, 0, 0, 0, 4, 0],
            [0, 0, 0, -1, -1, 2],
        ],
        dtype=_FLOAT64,
    )
    / 2.0
)
# The Jacobian of a cubic (in any dimension) as given by
# dB/ds = [-3 L1^2, 3 L1(L1 - 2 L2), 3 L2(2 L1 - L2), 3 L2^2, -6 L1 L3,
#          6 L3(L1 - L2), 6 L2 L3, -3 L3^2, 3 L3^2, 0] * nodes
# dB/dt = [-3 L1^2, -6 L1 L2, -3 L2^2, 0, 3 L1(L1 - 2 L3), 6 L2 (L1 - L3),
#          3 L2^2, 3 L3(2 L1 - L3), 6 L2 L3, 3 L3^2] * nodes
# We evaluate this at each of the 15 points in the quartic
# triangle and then stack them (2 columns * 15 = 30 columns)
_CUBIC_JACOBIAN_HELPER = (
    np.asfortranarray(
        [
            [
                -48,
                -48,
                -27,
                -27,
                -12,
                -12,
                -3,
                -3,
                0,
                0,
                -27,
                -27,
                -12,
                -12,
                -3,
                -3,
                0,
                0,
                -12,
                -12,
                -3,
                -3,
                0,
                0,
                -3,
                -3,
                0,
                0,
                0,
                0,
            ],
            [
                48,
                0,
                9,
                -18,
                -12,
                -24,
                -15,
                -18,
                0,
                0,
                27,
                0,
                0,
                -12,
                -9,
                -12,
                0,
                0,
                12,
                0,
                -3,
                -6,
                0,
                0,
                3,
                0,
                0,
                0,
                0,
                0,
            ],
            [
                0,
                0,
                15,
                -3,
                12,
                -12,
                -9,
                -27,
                -48,
                -48,
                0,
                0,
                9,
                -3,
                0,
                -12,
                -27,
                -27,
                0,
                0,
                3,
                -3,
                -12,
                -12,
                0,
                0,
                -3,
                -3,
                0,
                0,
            ],
            [
                0,
                0,
                3,
                0,
                12,
                0,
                27,
                0,
                48,
                0,
                0,
                0,
                3,
                0,
                12,
                0,
                27,
                0,
                0,
                0,
                3,
                0,
                12,
                0,
                0,
                0,
                3,
                0,
                0,
                0,
            ],
            [
                0,
                48,
                0,
                27,
                0,
                12,
                0,
                3,
                0,
                0,
                -18,
                9,
                -12,
                0,
                -6,
                -3,
                0,
                0,
                -24,
                -12,
                -12,
                -9,
                0,
                0,
                -18,
                -15,
                0,
                0,
                0,
                0,
            ],
            [
                0,
                0,
                0,
                18,
                0,
                24,
                0,
                18,
                0,
                0,
                18,
                0,
                6,
                6,
                -6,
                0,
                -18,
                -18,
                24,
                0,
                0,
                -6,
                -24,
                -24,
                18,
                0,
                -18,
                -18,
                0,
                0,
            ],
            [
                0,
                0,
                0,
                3,
                0,
                12,
                0,
                27,
                0,
                48,
                0,
                0,
                6,
                3,
                12,
                12,
                18,
                27,
                0,
                0,
                12,
                3,
                24,
                12,
                0,
                0,
                18,
                3,
                0,
                0,
            ],
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                -3,
                15,
                -3,
                9,
                -3,
                3,
                -3,
                -3,
                -12,
                12,
                -12,
                0,
                -12,
                -12,
                -27,
                -9,
                -27,
                -27,
                -48,
                -48,
            ],
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                3,
                0,
                3,
                6,
                3,
                12,
                3,
                18,
                12,
                0,
                12,
                12,
                12,
                24,
                27,
                0,
                27,
                18,
                48,
                0,
            ],
            [
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                0,
                3,
                0,
                3,
                0,
                3,
                0,
                3,
                0,
                12,
                0,
                12,
                0,
                12,
                0,
                27,
                0,
                27,
                0,
                48,
            ],
        ],
        dtype=_FLOAT64,
    )
    / 16.0
)
_QUARTIC_TO_BERNSTEIN = np.asfortranarray(
    [
        [36, -39, 26, -9, 0, -39, 26, -9, 0, 26, -9, 0, -9, 0, 0],
        [0, 144, -128, 48, 0, 0, -64, 32, 0, 0, 16, 0, 0, 0, 0],
        [0, -108, 240, -108, 0, 0, -24, -24, 0, 0, 12, 0, 0, 0, 0],
        [0, 48, -128, 144, 0, 0, 32, -64, 0, 0, 16, 0, 0, 0, 0],
        [0, -9, 26, -39, 36, 0, -9, 26, -39, 0, -9, 26, 0, -9, 0],
        [0, 0, 0, 0, 0, 144, -64, 16, 0, -128, 32, 0, 48, 0, 0],
        [0, 0, 0, 0, 0, 0, 288, -96, 0, 0, -96, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, -96, 288, 0, 0, -96, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 16, -64, 144, 0, 32, -128, 0, 48, 0],
        [0, 0, 0, 0, 0, -108, -24, 12, 0, 240, -24, 0, -108, 0, 0],
        [0, 0, 0, 0, 0, 0, -96, -96, 0, 0, 288, 0, 0, 0, 0],
        [0, 0, 0, 0, 0, 0, 12, -24, -108, 0, -24, 240, 0, -108, 0],
        [0, 0, 0, 0, 0, 48, 32, 16, 0, -128, -64, 0, 144, 0, 0],
        [0, 0, 0, 0, 0, 0, 16, 32, 48, 0, -64, -128, 0, 144, 0],
        [0, 0, 0, 0, 0, -9, -9, -9, -9, 26, 26, 26, -39, -39, 36],
    ],
    dtype=_FLOAT64,
)
# NOTE: We avoid round-off until after ``_QUARTIC_TO_BERNSTEIN``
#       has been applied.
_QUARTIC_BERNSTEIN_FACTOR = 36.0
# List of constants for ``basic_interior_combine()``. In each constant, each
# row is a return value of ``ends_to_curve()``. The second and third constant
# are just obtained from the first by rotating the rows.
FIRST_TRIANGLE_INFO = (
    ((0, 0.0, 1.0), (1, 0.0, 1.0), (2, 0.0, 1.0)),
    ((1, 0.0, 1.0), (2, 0.0, 1.0), (0, 0.0, 1.0)),
    ((2, 0.0, 1.0), (0, 0.0, 1.0), (1, 0.0, 1.0)),
)
SECOND_TRIANGLE_INFO = (
    ((3, 0.0, 1.0), (4, 0.0, 1.0), (5, 0.0, 1.0)),
    ((4, 0.0, 1.0), (5, 0.0, 1.0), (3, 0.0, 1.0)),
    ((5, 0.0, 1.0), (3, 0.0, 1.0), (4, 0.0, 1.0)),
)
# Threshold where a vector cross-product (u x v) is considered
# to be "zero". This is a "hack", since it doesn't take ||u||
# or ||v|| into account.
ALMOST_TANGENT = 0.5**50
# Hardcoded "line integral" helpers for ``shoelace_for_area()``.
SHOELACE_LINEAR = ((1, 0, 1),)
SHOELACE_QUADRATIC = ((2, 0, 1), (1, 0, 2), (2, 1, 2))
SHOELACE_CUBIC = (
    (6, 0, 1),
    (3, 0, 2),
    (1, 0, 3),
    (3, 1, 2),
    (3, 1, 3),
    (6, 2, 3),
)
SHOELACE_QUARTIC = (
    (20, 0, 1),
    (10, 0, 2),
    (4, 0, 3),
    (1, 0, 4),
    (8, 1, 2),
    (8, 1, 3),
    (4, 1, 4),
    (8, 2, 3),
    (10, 2, 4),
    (20, 3, 4),
)


def polynomial_sign(poly_triangle, degree):
    r"""Determine the "sign" of a polynomial on the reference triangle.

    .. note::

       This is used **only** by ``Triangle._compute_valid()`` (which is
       in turn used to compute / cache the :attr:`.Triangle.is_valid`
       property).

    Checks if a polynomial :math:`p(s, t)` is positive, negative
    or mixed sign on the reference triangle.

    Does this by utilizing the B |eacute| zier form of :math:`p`: it is a
    convex combination of the Bernstein basis (real numbers) hence
    if the Bernstein basis is all positive, the polynomial must be.

    If the values are mixed, then we can recursively subdivide
    until we are in a region where the coefficients are all one
    sign.

    Args:
        poly_triangle (numpy.ndarray): 2D array (with 1 row) of control
            points for a "triangle", i.e. a bivariate polynomial.
        degree (int): The degree of the triangle / polynomial given by
            ``poly_triangle``.

    Returns:
        int: The sign of the polynomial. Will be one of ``-1``, ``1``
        or ``0``. A value of ``0`` indicates a mixed sign or the
        zero polynomial.

    Raises:
        ValueError: If no conclusion is reached after the maximum
            number of subdivisions.
    """
    # The indices where the corner nodes in a triangle are.
    corner_indices = (0, degree, -1)
    sub_polys = [poly_triangle]
    signs = set()
    for _ in range(_MAX_POLY_SUBDIVISIONS):
        undecided = []
        for poly in sub_polys:
            # First add all the signs of the corner nodes.
            signs.update(_SIGN(poly[0, corner_indices]).astype(int))
            # Then check if the ``poly`` nodes are **uniformly** one sign.
            if np.all(poly == 0.0):
                signs.add(0)
            elif np.all(poly > 0.0):
                signs.add(1)
            elif np.all(poly < 0.0):
                signs.add(-1)
            else:
                undecided.append(poly)
            if len(signs) > 1:
                return 0

        sub_polys = functools.reduce(
            operator.add,
            [subdivide_nodes(poly, degree) for poly in undecided],
            (),
        )
        if not sub_polys:
            break

    if sub_polys:
        raise ValueError(
            "Did not reach a conclusion after max subdivisions",
            _MAX_POLY_SUBDIVISIONS,
        )

    # NOTE: We are guaranteed that ``len(signs) <= 1``.
    return signs.pop()


def two_by_two_det(mat):
    r"""Compute the determinant of a 2x2 matrix.

    .. note::

       This is used **only** by :func:`quadratic_jacobian_polynomial` and
       :func:`cubic_jacobian_polynomial`.

    This is "needed" because :func:`numpy.linalg.det` uses a more generic
    determinant implementation which can introduce rounding even when the
    simple :math:`a d - b c` will suffice. For example:

    Args:
        mat (numpy.ndarray): A 2x2 matrix.

    Returns:
        float: The determinant of ``mat``.
    """
    return mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]


def quadratic_jacobian_polynomial(nodes):
    r"""Compute the Jacobian determinant of a quadratic triangle.

    .. note::

       This is used **only** by ``Triangle._compute_valid()`` (which is
       in turn used to compute / cache the :attr:`.Triangle.is_valid`
       property).

    Converts :math:`\det(J(s, t))` to a polynomial on the reference
    triangle and represents it as a triangle object.

    .. note::

       This assumes that ``nodes`` is ``2 x 6`` but doesn't verify this.
       (However, the right multiplication by ``_QUADRATIC_JACOBIAN_HELPER``
       would fail if ``nodes`` wasn't ``R x 6`` and then the ensuing
       determinants would fail if there weren't 2 rows.)

    Args:
        nodes (numpy.ndarray): A 2 x 6 array of nodes in a triangle.

    Returns:
        numpy.ndarray: 1 x 6 array, coefficients in Bernstein basis.
    """
    # First evaluate the Jacobian at each of the 6 nodes.
    jac_parts = _py_helpers.matrix_product(nodes, _QUADRATIC_JACOBIAN_HELPER)
    jac_at_nodes = np.empty((1, 6), order="F")
    # pylint: disable=unsubscriptable-object
    jac_at_nodes[0, 0] = two_by_two_det(jac_parts[:, :2])
    jac_at_nodes[0, 1] = two_by_two_det(jac_parts[:, 2:4])
    jac_at_nodes[0, 2] = two_by_two_det(jac_parts[:, 4:6])
    jac_at_nodes[0, 3] = two_by_two_det(jac_parts[:, 6:8])
    jac_at_nodes[0, 4] = two_by_two_det(jac_parts[:, 8:10])
    jac_at_nodes[0, 5] = two_by_two_det(jac_parts[:, 10:])
    # pylint: enable=unsubscriptable-object
    # Convert the nodal values to the Bernstein basis...
    bernstein = _py_helpers.matrix_product(
        jac_at_nodes, _QUADRATIC_TO_BERNSTEIN
    )
    return bernstein


def cubic_jacobian_polynomial(nodes):
    r"""Compute the Jacobian determinant of a cubic triangle.

    .. note::

       This is used **only** by ``Triangle._compute_valid()`` (which is
       in turn used to compute / cache the :attr:`.Triangle.is_valid`
       property).

    Converts :math:`\det(J(s, t))` to a polynomial on the reference
    triangle and represents it as a triangle object.

    .. note::

       This assumes that ``nodes`` is ``2 x 10`` but doesn't verify this.
       (However, the right multiplication by ``_CUBIC_JACOBIAN_HELPER``
       would fail if ``nodes`` wasn't ``R x 10`` and then the ensuing
       determinants would fail if there weren't 2 rows.)

    Args:
        nodes (numpy.ndarray): A 2 x 10 array of nodes in a triangle.

    Returns:
        numpy.ndarray: 1 x 15 array, coefficients in Bernstein basis.
    """
    # First evaluate the Jacobian at each of the 15 nodes
    # in the quartic triangle.
    jac_parts = _py_helpers.matrix_product(nodes, _CUBIC_JACOBIAN_HELPER)
    jac_at_nodes = np.empty((1, 15), order="F")
    # pylint: disable=unsubscriptable-object
    jac_at_nodes[0, 0] = two_by_two_det(jac_parts[:, :2])
    jac_at_nodes[0, 1] = two_by_two_det(jac_parts[:, 2:4])
    jac_at_nodes[0, 2] = two_by_two_det(jac_parts[:, 4:6])
    jac_at_nodes[0, 3] = two_by_two_det(jac_parts[:, 6:8])
    jac_at_nodes[0, 4] = two_by_two_det(jac_parts[:, 8:10])
    jac_at_nodes[0, 5] = two_by_two_det(jac_parts[:, 10:12])
    jac_at_nodes[0, 6] = two_by_two_det(jac_parts[:, 12:14])
    jac_at_nodes[0, 7] = two_by_two_det(jac_parts[:, 14:16])
    jac_at_nodes[0, 8] = two_by_two_det(jac_parts[:, 16:18])
    jac_at_nodes[0, 9] = two_by_two_det(jac_parts[:, 18:20])
    jac_at_nodes[0, 10] = two_by_two_det(jac_parts[:, 20:22])
    jac_at_nodes[0, 11] = two_by_two_det(jac_parts[:, 22:24])
    jac_at_nodes[0, 12] = two_by_two_det(jac_parts[:, 24:26])
    jac_at_nodes[0, 13] = two_by_two_det(jac_parts[:, 26:28])
    jac_at_nodes[0, 14] = two_by_two_det(jac_parts[:, 28:])
    # pylint: enable=unsubscriptable-object
    # Convert the nodal values to the Bernstein basis...
    bernstein = _py_helpers.matrix_product(jac_at_nodes, _QUARTIC_TO_BERNSTEIN)
    bernstein /= _QUARTIC_BERNSTEIN_FACTOR
    return bernstein


def de_casteljau_one_round(nodes, degree, lambda1, lambda2, lambda3):
    r"""Performs one "round" of the de Casteljau algorithm for triangles.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    .. note::

       This is a helper function, used by :func:`make_transform` and
       :func:`specialize_triangle` (and :func:`make_transform` is **only**
       used by :func:`specialize_triangle`).

    Converts the ``nodes`` into a basis for a triangle one degree smaller
    by using the barycentric weights:

    .. math::

       q_{i, j, k} = \lambda_1 \cdot p_{i + 1, j, k} +
           \lambda_2 \cdot p_{i, j + 1, k} + \lambda_2 \cdot p_{i, j, k + 1}

    .. note:

       For degree :math:`d`, the number of nodes should be
       :math:`(d + 1)(d + 2)/2`, but we don't verify this.

    Args:
        nodes (numpy.ndarray): The nodes to reduce.
        degree (int): The degree of the triangle.
        lambda1 (float): Parameter along the reference triangle.
        lambda2 (float): Parameter along the reference triangle.
        lambda3 (float): Parameter along the reference triangle.

    Returns:
        numpy.ndarray: The converted nodes.
    """
    dimension, num_nodes = nodes.shape
    num_new_nodes = num_nodes - degree - 1
    new_nodes = np.empty((dimension, num_new_nodes), order="F")
    index = 0
    # parent_i1 = index + k
    # parent_i2 = index + k + 1
    # parent_i3 = index + degree + 1
    parent_i1 = 0
    parent_i2 = 1
    parent_i3 = degree + 1
    for k in range(degree):
        for unused_j in range(degree - k):
            # NOTE: i = (degree - 1) - j - k
            new_nodes[:, index] = (
                lambda1 * nodes[:, parent_i1]
                + lambda2 * nodes[:, parent_i2]
                + lambda3 * nodes[:, parent_i3]
            )
            # Update all the indices.
            parent_i1 += 1
            parent_i2 += 1
            parent_i3 += 1
            index += 1
        # Update the indices that depend on k.
        parent_i1 += 1
        parent_i2 += 1
    return new_nodes


def make_transform(degree, weights_a, weights_b, weights_c):
    """Compute matrices corresponding to the de Casteljau algorithm.

    .. note::

       This is a helper used only by :func:`specialize_triangle`.

    Applies the de Casteljau algorithm to the identity matrix, thus
    effectively caching the algorithm in a transformation matrix.

    .. note::

       This is premature optimization. It's unclear if the time
       saved from "caching" one round of de Casteljau is cancelled
       out by the extra storage required for the 3 matrices.

    Args:
        degree (int): The degree of a candidate triangle.
        weights_a (numpy.ndarray): Triple (1D array) of barycentric weights
            for a point in the reference triangle
        weights_b (numpy.ndarray): Triple (1D array) of barycentric weights
            for a point in the reference triangle
        weights_c (numpy.ndarray): Triple (1D array) of barycentric weights
            for a point in the reference triangle

    Returns:
        Mapping[int, numpy.ndarray]: Mapping from keys to the de Casteljau
        transformation mappings. The keys are ``0`` corresponding to
        ``weights_a``, ``1`` to ``weights_b`` and ``2`` to ``weights_c``.
    """
    num_nodes = ((degree + 1) * (degree + 2)) // 2
    id_mat = np.eye(num_nodes, order="F")
    # Pre-compute the matrices that do the reduction so we don't
    # have to **actually** perform the de Casteljau algorithm
    # every time.
    transform = {
        0: de_casteljau_one_round(id_mat, degree, *weights_a),
        1: de_casteljau_one_round(id_mat, degree, *weights_b),
        2: de_casteljau_one_round(id_mat, degree, *weights_c),
    }
    return transform


def reduced_to_matrix(shape, degree, vals_by_weight):
    r"""Converts a reduced values dictionary into a matrix.

    .. note::

       This is a helper used only by :func:`specialize_triangle`.

    The ``vals_by_weight`` mapping has keys of the form:
    ``(0, ..., 1, ..., 2, ...)`` where the ``0`` corresponds
    to the number of times the first set of barycentric
    weights was used in the reduction process, and similarly
    for ``1`` and ``2``.

    These points correspond to barycentric weights in their
    own right. For example ``(0, 0, 0, 1, 2, 2)`` corresponds to
    the barycentric weight
    :math:`\left(\frac{3}{6}, \frac{1}{6}, \frac{2}{6}\right)`.

    Once the keys in ``vals_by_weight`` have been converted
    to barycentric coordinates, we order them according to
    our rule (bottom to top, left to right) and then return
    them in a single matrix.

    Args:
        shape (tuple): The shape of the result matrix.
        degree (int): The degree of the triangle.
        vals_by_weight (Mapping[tuple, numpy.ndarray]): Dictionary
            of reduced nodes according to blending of each of the
            three sets of weights in a reduction.

    Returns:
        numpy.ndarray: The newly created reduced control points.
    """
    result = np.empty(shape, order="F")
    index = 0
    for k in range(degree + 1):
        for j in range(degree + 1 - k):
            i = degree - j - k
            key = (0,) * i + (1,) * j + (2,) * k
            result[:, index] = vals_by_weight[key][:, 0]
            index += 1
    return result


def specialize_triangle(nodes, degree, weights_a, weights_b, weights_c):
    """Specialize a triangle to a reparameterization.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Does so by taking three points (in barycentric form) within the
    reference triangle and then reparameterizing the triangle onto
    the triangle formed by those three points.

    .. note::

       This assumes the triangle is degree 1 or greater but doesn't check.

    .. note::

       This is used **only** as a helper for :func:`subdivide_nodes`, however
       it may be worth adding this to :class:`.Triangle` as an analogue to
       :meth:`.Curve.specialize`.

    Args:
        nodes (numpy.ndarray): Control points for a triangle.
        degree (int): The degree of the triangle.
        weights_a (numpy.ndarray): Triple (1D array) of barycentric weights
            for a point in the reference triangle
        weights_b (numpy.ndarray): Triple (1D array) of barycentric weights
            for a point in the reference triangle
        weights_c (numpy.ndarray): Triple (1D array) of barycentric weights
            for a point in the reference triangle

    Returns:
        numpy.ndarray: The control points for the specialized triangle.
    """
    # Uses A-->0, B-->1, C-->2 to represent the specialization used.
    partial_vals = {
        (0,): de_casteljau_one_round(nodes, degree, *weights_a),
        (1,): de_casteljau_one_round(nodes, degree, *weights_b),
        (2,): de_casteljau_one_round(nodes, degree, *weights_c),
    }
    for reduced_deg in range(degree - 1, 0, -1):
        new_partial = {}
        transform = make_transform(
            reduced_deg, weights_a, weights_b, weights_c
        )
        for key, sub_nodes in partial_vals.items():
            # Our keys are ascending so we increment from the last value.
            for next_id in range(key[-1], 2 + 1):
                new_key = key + (next_id,)
                new_partial[new_key] = _py_helpers.matrix_product(
                    sub_nodes, transform[next_id]
                )
        partial_vals = new_partial
    return reduced_to_matrix(nodes.shape, degree, partial_vals)


def subdivide_nodes(nodes, degree):
    """Subdivide a triangle into four sub-triangles.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Does so by taking the unit triangle (i.e. the domain of the triangle) and
    splitting it into four sub-triangles by connecting the midpoints of each
    side.

    Args:
        nodes (numpy.ndarray): Control points for a triangle.
        degree (int): The degree of the triangle.

    Returns:
        Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray, numpy.ndarray]: The
        nodes for the four sub-triangles.
    """
    if degree == 1:
        nodes_a = _py_helpers.matrix_product(nodes, LINEAR_SUBDIVIDE_A)
        nodes_b = _py_helpers.matrix_product(nodes, LINEAR_SUBDIVIDE_B)
        nodes_c = _py_helpers.matrix_product(nodes, LINEAR_SUBDIVIDE_C)
        nodes_d = _py_helpers.matrix_product(nodes, LINEAR_SUBDIVIDE_D)
    elif degree == 2:
        nodes_a = _py_helpers.matrix_product(nodes, QUADRATIC_SUBDIVIDE_A)
        nodes_b = _py_helpers.matrix_product(nodes, QUADRATIC_SUBDIVIDE_B)
        nodes_c = _py_helpers.matrix_product(nodes, QUADRATIC_SUBDIVIDE_C)
        nodes_d = _py_helpers.matrix_product(nodes, QUADRATIC_SUBDIVIDE_D)
    elif degree == 3:
        nodes_a = _py_helpers.matrix_product(nodes, CUBIC_SUBDIVIDE_A)
        nodes_b = _py_helpers.matrix_product(nodes, CUBIC_SUBDIVIDE_B)
        nodes_c = _py_helpers.matrix_product(nodes, CUBIC_SUBDIVIDE_C)
        nodes_d = _py_helpers.matrix_product(nodes, CUBIC_SUBDIVIDE_D)
    elif degree == 4:
        nodes_a = _py_helpers.matrix_product(nodes, QUARTIC_SUBDIVIDE_A)
        nodes_b = _py_helpers.matrix_product(nodes, QUARTIC_SUBDIVIDE_B)
        nodes_c = _py_helpers.matrix_product(nodes, QUARTIC_SUBDIVIDE_C)
        nodes_d = _py_helpers.matrix_product(nodes, QUARTIC_SUBDIVIDE_D)
    else:
        nodes_a = specialize_triangle(
            nodes,
            degree,
            _WEIGHTS_SUBDIVIDE0,
            _WEIGHTS_SUBDIVIDE1,
            _WEIGHTS_SUBDIVIDE2,
        )
        nodes_b = specialize_triangle(
            nodes,
            degree,
            _WEIGHTS_SUBDIVIDE3,
            _WEIGHTS_SUBDIVIDE2,
            _WEIGHTS_SUBDIVIDE1,
        )
        nodes_c = specialize_triangle(
            nodes,
            degree,
            _WEIGHTS_SUBDIVIDE1,
            _WEIGHTS_SUBDIVIDE4,
            _WEIGHTS_SUBDIVIDE3,
        )
        nodes_d = specialize_triangle(
            nodes,
            degree,
            _WEIGHTS_SUBDIVIDE2,
            _WEIGHTS_SUBDIVIDE3,
            _WEIGHTS_SUBDIVIDE5,
        )
    return nodes_a, nodes_b, nodes_c, nodes_d


def jacobian_s(nodes, degree, dimension):
    r"""Compute :math:`\frac{\partial B}{\partial s}`.

    .. note::

       This is a helper for :func:`jacobian_both`, which has an
       equivalent Fortran implementation.

    Args:
        nodes (numpy.ndarray): Array of nodes in a triangle.
        degree (int): The degree of the triangle.
        dimension (int): The dimension the triangle lives in.

    Returns:
        numpy.ndarray: Nodes of the Jacobian triangle in
        B |eacute| zier form.
    """
    num_nodes = (degree * (degree + 1)) // 2
    result = np.empty((dimension, num_nodes), order="F")
    index = 0
    i = 0
    for num_vals in range(degree, 0, -1):
        for _ in range(num_vals):
            result[:, index] = nodes[:, i + 1] - nodes[:, i]
            # Update the indices
            index += 1
            i += 1
        # In between each row, the index gains an extra value.
        i += 1
    return float(degree) * result


def jacobian_t(nodes, degree, dimension):
    r"""Compute :math:`\frac{\partial B}{\partial t}`.

    .. note::

       This is a helper for :func:`jacobian_both`, which has an
       equivalent Fortran implementation.

    Args:
        nodes (numpy.ndarray): Array of nodes in a triangle.
        degree (int): The degree of the triangle.
        dimension (int): The dimension the triangle lives in.

    Returns:
        numpy.ndarray: Nodes of the Jacobian triangle in
        B |eacute| zier form.
    """
    num_nodes = (degree * (degree + 1)) // 2
    result = np.empty((dimension, num_nodes), order="F")
    index = 0
    i = 0
    j = degree + 1
    for num_vals in range(degree, 0, -1):
        for _ in range(num_vals):
            result[:, index] = nodes[:, j] - nodes[:, i]
            # Update the indices
            index += 1
            i += 1
            j += 1
        # In between each row, the index gains an extra value.
        i += 1
    return float(degree) * result


def jacobian_both(nodes, degree, dimension):
    r"""Compute :math:`s` and :math:`t` partial of :math:`B`.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Array of nodes in a triangle.
        degree (int): The degree of the triangle.
        dimension (int): The dimension the triangle lives in.

    Returns:
        numpy.ndarray: Nodes of the Jacobian triangles in
        B |eacute| zier form.
    """
    _, num_nodes = nodes.shape
    result = np.empty((2 * dimension, num_nodes - degree - 1), order="F")
    result[:dimension, :] = jacobian_s(nodes, degree, dimension)
    result[dimension:, :] = jacobian_t(nodes, degree, dimension)
    return result


def jacobian_det(nodes, degree, st_vals):
    r"""Compute :math:`\det(D B)` at a set of values.

    This requires that :math:`B \in \mathbf{R}^2`.

    .. note::

       This assumes but does not check that each ``(s, t)``
       in ``st_vals`` is inside the reference triangle.

    .. warning::

       This relies on helpers in :mod:`bezier` for computing the
       Jacobian of the triangle. However, these helpers are not
       part of the public API and may change or be removed.

    .. testsetup:: jacobian-det

       import numpy as np

       import bezier
       from bezier.hazmat.triangle_helpers import jacobian_det

    .. doctest:: jacobian-det
       :options: +NORMALIZE_WHITESPACE

       >>> import bezier
       >>> import numpy as np
       >>> nodes = np.asfortranarray([
       ...     [0.0, 1.0, 2.0, 0.0, 1.5, 0.0],
       ...     [0.0, 0.0, 0.0, 1.0, 1.5, 2.0],
       ... ])
       >>> triangle = bezier.Triangle(nodes, degree=2)
       >>> st_vals = np.asfortranarray([
       ...     [0.25, 0.0  ],
       ...     [0.75, 0.125],
       ...     [0.5 , 0.5  ],
       ... ])
       >>> s_vals, t_vals = st_vals.T
       >>> triangle.evaluate_cartesian_multi(st_vals)
       array([[0.5 , 1.59375, 1.25 ],
              [0.  , 0.34375, 1.25 ]])
       >>> # B(s, t) = [s(t + 2), t(s + 2)]
       >>> s_vals * (t_vals + 2)
       array([0.5 , 1.59375, 1.25 ])
       >>> t_vals * (s_vals + 2)
       array([0. , 0.34375, 1.25 ])
       >>> jacobian_det(nodes, 2, st_vals)
       array([4.5 , 5.75, 6. ])
       >>> # det(DB) = 2(s + t + 2)
       >>> 2 * (s_vals + t_vals + 2)
       array([4.5 , 5.75, 6. ])

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Nodes defining a B |eacute| zier
            triangle :math:`B(s, t)`.
        degree (int): The degree of the triangle :math:`B`.
        st_vals (numpy.ndarray): ``N x 2`` array of Cartesian
            inputs to B |eacute| zier triangles defined by
            :math:`B_s` and :math:`B_t`.

    Returns:
        numpy.ndarray: Array of all determinant values, one
        for each row in ``st_vals``.
    """
    jac_nodes = jacobian_both(nodes, degree, 2)
    if degree == 1:
        num_vals, _ = st_vals.shape
        bs_bt_vals = np.repeat(jac_nodes, num_vals, axis=1)
    else:
        bs_bt_vals = evaluate_cartesian_multi(
            jac_nodes, degree - 1, st_vals, 4
        )
    # Take the determinant for each (s, t).
    return (
        bs_bt_vals[0, :] * bs_bt_vals[3, :]
        - bs_bt_vals[1, :] * bs_bt_vals[2, :]
    )


def classify_tangent_intersection(
    intersection, nodes1, tangent1, nodes2, tangent2
):
    """Helper for :func:`classify_intersection` at tangencies.

    .. note::

       This is a helper used only by :func:`classify_intersection`.

    Args:
        intersection (Intersection): An intersection object.
        nodes1 (numpy.ndarray): Control points for the first curve at
            the intersection.
        tangent1 (numpy.ndarray): The tangent vector to the first curve
            at the intersection (``2 x 1`` array).
        nodes2 (numpy.ndarray): Control points for the second curve at
            the intersection.
        tangent2 (numpy.ndarray): The tangent vector to the second curve
            at the intersection (``2 x 1`` array).

    Returns:
        IntersectionClassification: The "inside" curve type, based on
        the classification enum. Will either be ``opposed`` or one
        of the ``tangent`` values.

    Raises:
        NotImplementedError: If the curves are tangent at the intersection
            and have the same curvature.
    """
    # Each array is 2 x 1 (i.e. a column vector), we want the vector
    # dot product.
    dot_prod = np.vdot(tangent1[:, 0], tangent2[:, 0])
    # NOTE: When computing curvatures we assume that we don't have lines
    #       here, because lines that are tangent at an intersection are
    #       parallel and we don't handle that case.
    curvature1 = curve_helpers.get_curvature(nodes1, tangent1, intersection.s)
    curvature2 = curve_helpers.get_curvature(nodes2, tangent2, intersection.t)
    if dot_prod < 0:
        # If the tangent vectors are pointing in the opposite direction,
        # then the curves are facing opposite directions.
        # pylint: disable=assignment-from-no-return
        sign1, sign2 = _SIGN([curvature1, curvature2])
        # pylint: enable=assignment-from-no-return
        if sign1 == sign2:
            # If both curvatures are positive, since the curves are
            # moving in opposite directions, the tangency isn't part of
            # the triangle intersection.
            if sign1 == 1.0:
                return CLASSIFICATION_T.OPPOSED

            else:
                return CLASSIFICATION_T.TANGENT_BOTH

        else:
            delta_c = abs(curvature1) - abs(curvature2)
            if delta_c == 0.0:
                raise NotImplementedError(_SAME_CURVATURE)

            if sign1 == _SIGN(delta_c):
                return CLASSIFICATION_T.OPPOSED

            return CLASSIFICATION_T.TANGENT_BOTH

    else:
        if curvature1 > curvature2:
            return CLASSIFICATION_T.TANGENT_FIRST

        elif curvature1 < curvature2:
            return CLASSIFICATION_T.TANGENT_SECOND

        else:
            raise NotImplementedError(_SAME_CURVATURE)


def ignored_edge_corner(edge_tangent, corner_tangent, corner_previous_edge):
    """Check ignored when a corner lies **inside** another edge.

    .. note::

       This is a helper used only by :func:`ignored_corner`, which in turn is
       only used by :func:`classify_intersection`.

    Helper for :func:`ignored_corner` where one of ``s`` and
    ``t`` are ``0``, but **not both**.

    Args:
        edge_tangent (numpy.ndarray): Tangent vector (``2 x 1`` array) along
            the edge that the intersection occurs in the middle of.
        corner_tangent (numpy.ndarray): Tangent vector (``2 x 1`` array) at
            the corner where intersection occurs (at the beginning of edge).
        corner_previous_edge (numpy.ndarray): Edge that ends at the corner
            intersection (whereas ``corner_tangent`` comes from the edge
            that **begins** at the corner intersection). This is a ``2 x N``
            array where ``N`` is the number of nodes on the edge.

    Returns:
        bool: Indicates if the corner intersection should be ignored.
    """
    cross_prod = _py_helpers.cross_product(
        edge_tangent.ravel(order="F"), corner_tangent.ravel(order="F")
    )
    # A negative cross product indicates that ``edge_tangent`` is
    # "inside" / "to the left" of ``corner_tangent`` (due to right-hand rule).
    if cross_prod > 0.0:
        return False

    # Do the same for the **other** tangent at the corner.
    alt_corner_tangent = curve_helpers.evaluate_hodograph(
        1.0, corner_previous_edge
    )
    # Change the direction of the "in" tangent so that it points "out".
    alt_corner_tangent *= -1.0
    cross_prod = _py_helpers.cross_product(
        edge_tangent.ravel(order="F"), alt_corner_tangent.ravel(order="F")
    )
    return cross_prod <= 0.0


def ignored_double_corner(
    intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2
):
    """Check if an intersection is an "ignored" double corner.

    .. note::

       This is a helper used only by :func:`ignored_corner`, which in turn is
       only used by :func:`classify_intersection`.

    Helper for :func:`ignored_corner` where both ``s`` and
    ``t`` are ``0``.

    Does so by checking if either edge through the ``t`` corner goes
    through the interior of the other triangle. An interior check
    is done by checking that a few cross products are positive.

    Args:
        intersection (Intersection): An intersection to "diagnose".
        tangent_s (numpy.ndarray): The tangent vector (``2 x 1`` array) to
            the first curve at the intersection.
        tangent_t (numpy.ndarray): The tangent vector (``2 x 1`` array) to
            the second curve at the intersection.
        edge_nodes1 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the first triangle being intersected.
        edge_nodes2 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the second triangle being intersected.

    Returns:
        bool: Indicates if the corner is to be ignored.
    """
    # Compute the other edge for the ``s`` triangle.
    prev_index = (intersection.index_first - 1) % 3
    prev_edge = edge_nodes1[prev_index]
    alt_tangent_s = curve_helpers.evaluate_hodograph(1.0, prev_edge)
    # First check if ``tangent_t`` is interior to the ``s`` triangle.
    cross_prod1 = _py_helpers.cross_product(
        tangent_s.ravel(order="F"), tangent_t.ravel(order="F")
    )
    # A positive cross product indicates that ``tangent_t`` is
    # interior to ``tangent_s``. Similar for ``alt_tangent_s``.
    # If ``tangent_t`` is interior to both, then the triangles
    # do more than just "kiss" at the corner, so the corner should
    # not be ignored.
    if cross_prod1 >= 0.0:
        # Only compute ``cross_prod2`` if we need to.
        cross_prod2 = _py_helpers.cross_product(
            alt_tangent_s.ravel(order="F"), tangent_t.ravel(order="F")
        )
        if cross_prod2 >= 0.0:
            return False

    # If ``tangent_t`` is not interior, we check the other ``t``
    # edge that ends at the corner.
    prev_index = (intersection.index_second - 1) % 3
    prev_edge = edge_nodes2[prev_index]
    alt_tangent_t = curve_helpers.evaluate_hodograph(1.0, prev_edge)
    # Change the direction of the "in" tangent so that it points "out".
    alt_tangent_t *= -1.0
    cross_prod3 = _py_helpers.cross_product(
        tangent_s.ravel(order="F"), alt_tangent_t.ravel(order="F")
    )
    if cross_prod3 >= 0.0:
        # Only compute ``cross_prod4`` if we need to.
        cross_prod4 = _py_helpers.cross_product(
            alt_tangent_s.ravel(order="F"), alt_tangent_t.ravel(order="F")
        )
        if cross_prod4 >= 0.0:
            return False

    # If neither of ``tangent_t`` or ``alt_tangent_t`` are interior
    # to the ``s`` triangle, one of two things is true. Either
    # the two triangles have no interior intersection (1) or the
    # ``s`` triangle is bounded by both edges of the ``t`` triangle
    # at the corner intersection (2). To detect (2), we only need
    # check if ``tangent_s`` is interior to both ``tangent_t``
    # and ``alt_tangent_t``. ``cross_prod1`` contains
    # (tangent_s) x (tangent_t), so it's negative will tell if
    # ``tangent_s`` is interior. Similarly, ``cross_prod3``
    # contains (tangent_s) x (alt_tangent_t), but we also reversed
    # the sign on ``alt_tangent_t`` so switching the sign back
    # and reversing the arguments in the cross product cancel out.
    return cross_prod1 > 0.0 or cross_prod3 < 0.0


def ignored_corner(
    intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2
):
    """Check if an intersection is an "ignored" corner.

    .. note::

       This is a helper used only by :func:`classify_intersection`.

    An "ignored" corner is one where the triangles just "kiss" at
    the point of intersection but their interiors do not meet.

    We can determine this by comparing the tangent lines from
    the point of intersection.

    .. note::

       This assumes the ``intersection`` has been shifted to the
       beginning of a curve so only checks if ``s == 0.0`` or ``t == 0.0``
       (rather than also checking for ``1.0``).

    Args:
        intersection (Intersection): An intersection to "diagnose".
        tangent_s (numpy.ndarray): The tangent vector (``2 x 1`` array) to
            the first curve at the intersection.
        tangent_t (numpy.ndarray): The tangent vector (``2 x 1`` array) to
            the second curve at the intersection.
        edge_nodes1 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the first triangle being intersected.
        edge_nodes2 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the second triangle being intersected.

    Returns:
        bool: Indicates if the corner is to be ignored.
    """
    if intersection.s == 0.0:
        if intersection.t == 0.0:
            # Double corner.
            return ignored_double_corner(
                intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2
            )

        else:
            # s-only corner.
            prev_index = (intersection.index_first - 1) % 3
            prev_edge = edge_nodes1[prev_index]
            return ignored_edge_corner(tangent_t, tangent_s, prev_edge)

    elif intersection.t == 0.0:
        # t-only corner.
        prev_index = (intersection.index_second - 1) % 3
        prev_edge = edge_nodes2[prev_index]
        return ignored_edge_corner(tangent_s, tangent_t, prev_edge)

    else:
        # Not a corner.
        return False


def classify_intersection(intersection, edge_nodes1, edge_nodes2):
    r"""Determine which curve is on the "inside of the intersection".

    .. note::

       This is a helper used only by :meth:`.Triangle.intersect`.

    This is intended to be a helper for forming a :class:`.CurvedPolygon`
    from the edge intersections of two :class:`.Triangle`-s. In order
    to move from one intersection to another (or to the end of an edge),
    the interior edge must be determined at the point of intersection.

    The "typical" case is on the interior of both edges:

    .. image:: ../../images/classify_intersection1.png
       :align: center

    .. testsetup:: classify-intersection1, classify-intersection2,
                   classify-intersection3, classify-intersection4,
                   classify-intersection5, classify-intersection6,
                   classify-intersection7, classify-intersection8,
                   classify-intersection9

       import numpy as np
       import bezier
       from bezier.hazmat import curve_helpers
       from bezier.hazmat.intersection_helpers import Intersection
       from bezier.hazmat.triangle_helpers import classify_intersection

       def hodograph(curve, s):
           return curve_helpers.evaluate_hodograph(
               s, curve._nodes)

       def curvature(curve, s):
           nodes = curve._nodes
           tangent = curve_helpers.evaluate_hodograph(
               s, nodes)
           return curve_helpers.get_curvature(
               nodes, tangent, s)

    .. doctest:: classify-intersection1
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.asfortranarray([
       ...     [1.0, 1.75, 2.0],
       ...     [0.0, 0.25, 1.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [0.0, 1.6875, 2.0],
       ...     [0.0, 0.0625, 0.5],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.25, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> tangent1 = hodograph(curve1, s)
       >>> tangent1
       array([[1.25],
              [0.75]])
       >>> tangent2 = hodograph(curve2, t)
       >>> tangent2
       array([[2. ],
              [0.5]])
       >>> intersection = Intersection(0, s, 0, t)
       >>> edge_nodes1 = (nodes1, None, None)
       >>> edge_nodes2 = (nodes2, None, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       <IntersectionClassification.FIRST: 0>

    .. testcleanup:: classify-intersection1

       import make_images
       make_images.classify_intersection1(
           s, curve1, tangent1, curve2, tangent2)

    We determine the interior (i.e. left) one by using the `right-hand rule`_:
    by embedding the tangent vectors in :math:`\mathbf{R}^3`, we
    compute

    .. _right-hand rule: https://en.wikipedia.org/wiki/Right-hand_rule

    .. math::

       \left[\begin{array}{c}
           x_1'(s) \\ y_1'(s) \\ 0 \end{array}\right] \times
       \left[\begin{array}{c}
           x_2'(t) \\ y_2'(t) \\ 0 \end{array}\right] =
       \left[\begin{array}{c}
           0 \\ 0 \\ x_1'(s) y_2'(t) - x_2'(t) y_1'(s) \end{array}\right].

    If the cross product quantity
    :math:`B_1'(s) \times B_2'(t) = x_1'(s) y_2'(t) - x_2'(t) y_1'(s)`
    is positive, then the first curve is "outside" / "to the right", i.e.
    the second curve is interior. If the cross product is negative, the
    first curve is interior.

    When :math:`B_1'(s) \times B_2'(t) = 0`, the tangent
    vectors are parallel, i.e. the intersection is a point of tangency:

    .. image:: ../../images/classify_intersection2.png
       :align: center

    .. doctest:: classify-intersection2
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.asfortranarray([
       ...     [1.0, 1.5, 2.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [0.0, 1.5, 3.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> intersection = Intersection(0, s, 0, t)
       >>> edge_nodes1 = (nodes1, None, None)
       >>> edge_nodes2 = (nodes2, None, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       <IntersectionClassification.TANGENT_SECOND: 4>

    .. testcleanup:: classify-intersection2

       import make_images
       make_images.classify_intersection2(s, curve1, curve2)

    Depending on the direction of the parameterizations, the interior
    curve may change, but we can use the (signed) `curvature`_ of each
    curve at that point to determine which is on the interior:

    .. _curvature: https://en.wikipedia.org/wiki/Curvature

    .. image:: ../../images/classify_intersection3.png
       :align: center

    .. doctest:: classify-intersection3
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.asfortranarray([
       ...     [2.0, 1.5, 1.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [3.0, 1.5, 0.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> intersection = Intersection(0, s, 0, t)
       >>> edge_nodes1 = (nodes1, None, None)
       >>> edge_nodes2 = (nodes2, None, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       <IntersectionClassification.TANGENT_FIRST: 3>

    .. testcleanup:: classify-intersection3

       import make_images
       make_images.classify_intersection3(s, curve1, curve2)

    When the curves are moving in opposite directions at a point
    of tangency, there is no side to choose. Either the point of tangency
    is not part of any :class:`.CurvedPolygon` intersection

    .. image:: ../../images/classify_intersection4.png
       :align: center

    .. doctest:: classify-intersection4
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.asfortranarray([
       ...     [2.0, 1.5, 1.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [0.0, 1.5, 3.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> intersection = Intersection(0, s, 0, t)
       >>> edge_nodes1 = (nodes1, None, None)
       >>> edge_nodes2 = (nodes2, None, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       <IntersectionClassification.OPPOSED: 2>

    .. testcleanup:: classify-intersection4

       import make_images
       make_images.classify_intersection4(s, curve1, curve2)

    or the point of tangency is a "degenerate" part of two
    :class:`.CurvedPolygon` intersections. It is "degenerate"
    because from one direction, the point should be classified as
    :attr:`~.IntersectionClassification.FIRST` and from another as
    :attr:`~.IntersectionClassification.SECOND`.

    .. image:: ../../images/classify_intersection5.png
       :align: center

    .. doctest:: classify-intersection5
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.asfortranarray([
       ...     [1.0, 1.5, 2.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [3.0, 1.5, 0.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> intersection = Intersection(0, s, 0, t)
       >>> edge_nodes1 = (nodes1, None, None)
       >>> edge_nodes2 = (nodes2, None, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       <IntersectionClassification.TANGENT_BOTH: 6>

    .. testcleanup:: classify-intersection5

       import make_images
       make_images.classify_intersection5(s, curve1, curve2)

    The :attr:`~.IntersectionClassification.TANGENT_BOTH` classification
    can also occur if the curves are "kissing" but share a zero width
    interior at the point of tangency:

    .. image:: ../../images/classify_intersection9.png
       :align: center

    .. doctest:: classify-intersection9
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.asfortranarray([
       ...     [0.0, 20.0, 40.0],
       ...     [0.0, 40.0,  0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [40.0, 20.0,  0.0],
       ...     [40.0,  0.0, 40.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> intersection = Intersection(0, s, 0, t)
       >>> edge_nodes1 = (nodes1, None, None)
       >>> edge_nodes2 = (nodes2, None, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       <IntersectionClassification.TANGENT_BOTH: 6>

    .. testcleanup:: classify-intersection9

       import make_images
       make_images.classify_intersection9(s, curve1, curve2)

    However, if the `curvature`_ of each curve is identical, we
    don't try to distinguish further:

    .. image:: ../../images/classify_intersection6.png
       :align: center

    .. doctest:: classify-intersection6
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.asfortranarray([
       ...     [-0.125 , -0.125 , 0.375 ],
       ...     [ 0.0625, -0.0625, 0.0625],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [-0.25, -0.25, 0.75],
       ...     [ 0.25, -0.25, 0.25],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> hodograph(curve1, s)
       array([[0.5],
              [0. ]])
       >>> hodograph(curve2, t)
       array([[1.],
              [0.]])
       >>> float(curvature(curve1, s))
       2.0
       >>> float(curvature(curve2, t))
       2.0
       >>> intersection = Intersection(0, s, 0, t)
       >>> edge_nodes1 = (nodes1, None, None)
       >>> edge_nodes2 = (nodes2, None, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       Traceback (most recent call last):
         ...
       NotImplementedError: Tangent curves have same curvature.

    .. testcleanup:: classify-intersection6

       import make_images
       make_images.classify_intersection6(s, curve1, curve2)

    In addition to points of tangency, intersections that happen at
    the end of an edge need special handling:

    .. image:: ../../images/classify_intersection7.png
       :align: center

    .. doctest:: classify-intersection7
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1a = np.asfortranarray([
       ...     [0.0, 4.5, 9.0 ],
       ...     [0.0, 0.0, 2.25],
       ... ])
       >>> curve1a = bezier.Curve(nodes1a, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [11.25, 9.0, 2.75],
       ...     [ 0.0 , 4.5, 1.0 ],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 1.0, 0.375
       >>> curve1a.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> intersection = Intersection(0, s, 0, t)
       >>> edge_nodes1 = (nodes1a, None, None)
       >>> edge_nodes2 = (nodes2, None, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       Traceback (most recent call last):
         ...
       ValueError: ('Intersection occurs at the end of an edge',
                    's', 1.0, 't', 0.375)
       >>>
       >>> nodes1b = np.asfortranarray([
       ...     [9.0, 4.5, 0.0],
       ...     [2.25, 2.375, 2.5],
       ... ])
       >>> curve1b = bezier.Curve(nodes1b, degree=2)
       >>> curve1b.evaluate(0.0) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> intersection = Intersection(1, 0.0, 0, t)
       >>> edge_nodes1 = (nodes1a, nodes1b, None)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       <IntersectionClassification.FIRST: 0>

    .. testcleanup:: classify-intersection7

       import make_images
       make_images.classify_intersection7(s, curve1a, curve1b, curve2)

    As above, some intersections at the end of an edge are part of
    an actual intersection. However, some triangles may just "kiss" at a
    corner intersection:

    .. image:: ../../images/classify_intersection8.png
       :align: center

    .. doctest:: classify-intersection8
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.asfortranarray([
       ...     [0.25, 0.0, 0.0, 0.625, 0.5  , 1.0 ],
       ...     [1.0 , 0.5, 0.0, 0.875, 0.375, 0.75],
       ... ])
       >>> triangle1 = bezier.Triangle(nodes1, degree=2)
       >>> nodes2 = np.asfortranarray([
       ...     [0.0625, -0.25, -1.0, -0.5  , -1.0, -1.0],
       ...     [0.5   ,  1.0 ,  1.0,  0.125,  0.5,  0.0],
       ... ])
       >>> triangle2 = bezier.Triangle(nodes2, degree=2)
       >>> curve1, _, _ = triangle1.edges
       >>> edge_nodes1 = [curve.nodes for curve in triangle1.edges]
       >>> curve2, _, _ = triangle2.edges
       >>> edge_nodes2 = [curve.nodes for curve in triangle2.edges]
       >>> s, t = 0.5, 0.0
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True],
              [ True]])
       >>> intersection = Intersection(0, s, 0, t)
       >>> classify_intersection(intersection, edge_nodes1, edge_nodes2)
       <IntersectionClassification.IGNORED_CORNER: 5>

    .. testcleanup:: classify-intersection8

       import make_images
       make_images.classify_intersection8(
           s, curve1, triangle1, curve2, triangle2)

    .. note::

       This assumes the intersection occurs in :math:`\mathbf{R}^2`
       but doesn't check this.

    .. note::

       This function doesn't allow wiggle room / round-off when checking
       endpoints, nor when checking if the cross product is near zero,
       nor when curvatures are compared. However, the most "correct"
       version of this function likely should allow for some round off.

    Args:
        intersection (Intersection): An intersection object.
        edge_nodes1 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the first triangle being intersected.
        edge_nodes2 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the second triangle being intersected.

    Returns:
        IntersectionClassification: The "inside" curve type, based on
        the classification enum.

    Raises:
        ValueError: If the intersection occurs at the end of either
            curve involved. This is because we want to classify which
            curve to **move forward** on, and we can't move past the
            end of a segment.
    """
    if intersection.s == 1.0 or intersection.t == 1.0:
        raise ValueError(
            "Intersection occurs at the end of an edge",
            "s",
            intersection.s,
            "t",
            intersection.t,
        )

    nodes1 = edge_nodes1[intersection.index_first]
    tangent1 = curve_helpers.evaluate_hodograph(intersection.s, nodes1)
    nodes2 = edge_nodes2[intersection.index_second]
    tangent2 = curve_helpers.evaluate_hodograph(intersection.t, nodes2)
    if ignored_corner(
        intersection, tangent1, tangent2, edge_nodes1, edge_nodes2
    ):
        return CLASSIFICATION_T.IGNORED_CORNER

    # Take the cross product of tangent vectors to determine which one
    # is more "inside" / "to the left".
    cross_prod = _py_helpers.cross_product(
        tangent1.ravel(order="F"), tangent2.ravel(order="F")
    )
    if cross_prod < -ALMOST_TANGENT:
        return CLASSIFICATION_T.FIRST

    elif cross_prod > ALMOST_TANGENT:
        return CLASSIFICATION_T.SECOND

    else:
        # NOTE: A more robust approach would take ||tangent1|| and ||tangent2||
        #       into account when comparing (tangent1 x tangent2) to the
        #       "almost zero" threshold. We (for now) avoid doing this because
        #       normalizing the tangent vectors has a "cost" of ~6 flops each
        #       and that cost would happen for **every** single intersection.
        return classify_tangent_intersection(
            intersection, nodes1, tangent1, nodes2, tangent2
        )


def handle_ends(index1, s, index2, t):
    """Updates intersection parameters if it is on the end of an edge.

    .. note::

       This is a helper used only by :meth:`.Triangle.intersect`.

    Does nothing if the intersection happens in the middle of two
    edges.

    If the intersection occurs at the end of the first curve,
    moves it to the beginning of the next edge. Similar for the
    second curve.

    This function is used as a pre-processing step before passing
    an intersection to :func:`classify_intersection`. There, only
    corners that **begin** an edge are considered, since that
    function is trying to determine which edge to **move forward** on.

    Args:
        index1 (int): The index (among 0, 1, 2) of the first edge in the
            intersection.
        s (float): The parameter along the first curve of the intersection.
        index2 (int): The index (among 0, 1, 2) of the second edge in the
            intersection.
        t (float): The parameter along the second curve of the intersection.

    Returns:
        Tuple[bool, bool, Tuple[int, float, int, float]]: A triple of:

        * flag indicating if the intersection is at the end of an edge
        * flag indicating if the intersection is a "corner"
        * 4-tuple of the "updated" values ``(index1, s, index2, t)``
    """
    edge_end = False
    if s == 1.0:
        s = 0.0
        index1 = (index1 + 1) % 3
        edge_end = True
    # NOTE: This is not a typo, the values can be updated twice if both ``s``
    #       and ``t`` are ``1.0``
    if t == 1.0:
        t = 0.0
        index2 = (index2 + 1) % 3
        edge_end = True
    is_corner = s == 0.0 or t == 0.0
    return edge_end, is_corner, (index1, s, index2, t)


def to_front(intersection, intersections, unused):
    """Rotates a node to the "front".

    .. note::

       This is a helper used only by :func:`basic_interior_combine`, which in
       turn is only used by :func:`combine_intersections`.

    If a node is at the end of a segment, moves it to the beginning
    of the next segment (at the exact same point). We assume that
    callers have pruned ``intersections`` so that there are none
    with ``s == 1.0`` or ``t == 1.0``. Hence, any such intersection
    will be an "artificial" intersection added by :func:`get_next`.

    .. note::

        This method checks for **exact** endpoints, i.e. parameter
        bitwise identical to ``1.0``. But it may make sense to allow
        some wiggle room.

    Args:
        intersection (Intersection): The current intersection.
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of all detected intersections, provided as a reference for
            of potential points to arrive at.
        unused (List[~bezier.hazmat.intersection_helpers.Intersection]): List
            nodes that haven't been used yet in an intersection curved polygon

    Returns:
        Intersection: An intersection to (maybe) move to the beginning
        of the next edge of the triangle.
    """
    if intersection.s == 1.0:
        next_index = (intersection.index_first + 1) % 3
        # Make sure we haven't accidentally ignored an existing intersection.
        for other_int in intersections:
            if other_int.s == 0.0 and other_int.index_first == next_index:
                if other_int in unused:
                    unused.remove(other_int)
                return other_int

        # If we haven't already returned, create **another** artificial
        # intersection.
        return intersection_helpers.Intersection(
            next_index, 0.0, None, None, interior_curve=CLASSIFICATION_T.FIRST
        )

    elif intersection.t == 1.0:
        # NOTE: We assume, but do not check, that ``s == 1.0`` and ``t == 1.0``
        #       are mutually exclusive.
        next_index = (intersection.index_second + 1) % 3
        # Make sure we haven't accidentally ignored an existing intersection.
        for other_int in intersections:
            if other_int.t == 0.0 and other_int.index_second == next_index:
                if other_int in unused:
                    unused.remove(other_int)
                return other_int

        # If we haven't already returned, create **another** artificial
        # intersection.
        return intersection_helpers.Intersection(
            None, None, next_index, 0.0, interior_curve=CLASSIFICATION_T.SECOND
        )

    else:
        return intersection


def get_next_first(intersection, intersections, to_end=True):
    """Gets the next node along the current (first) edge.

    .. note::

       This is a helper used only by :func:`get_next`, which in
       turn is only used by :func:`basic_interior_combine`, which itself
       is only used by :func:`combine_intersections`.

    Along with :func:`get_next_second`, this function does the majority of the
    heavy lifting in :func:`get_next`. **Very** similar to
    :func:`get_next_second`, but this works with the first curve while the
    other function works with the second.

    Args:
        intersection (Intersection): The current intersection.
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of all detected intersections, provided as a reference for
            potential points to arrive at.
        to_end (Optional[bool]): Indicates if the next node should just be
            the end of the first edge or :data:`None`.

    Returns:
        point along a triangle of intersection. This will produce the next
        Optional[~bezier.hazmat.intersection_helpers.Intersection]: The "next"
        intersection along the current (first) edge or the end of the same
        edge. If ``to_end`` is :data:`False` and there are no other
        intersections along the current edge, will return :data:`None` (rather
        than the end of the same edge).
    """
    along_edge = None
    index_first = intersection.index_first
    s = intersection.s
    for other_int in intersections:
        other_s = other_int.s
        if other_int.index_first == index_first and other_s > s:
            if along_edge is None or other_s < along_edge.s:
                along_edge = other_int
    if along_edge is None:
        if to_end:
            # If there is no other intersection on the edge, just return
            # the segment end.
            return intersection_helpers.Intersection(
                index_first,
                1.0,
                None,
                None,
                interior_curve=CLASSIFICATION_T.FIRST,
            )

        else:
            return None

    else:
        return along_edge


def get_next_second(intersection, intersections, to_end=True):
    """Gets the next node along the current (second) edge.

    .. note::

       This is a helper used only by :func:`get_next`, which in
       turn is only used by :func:`basic_interior_combine`, which itself
       is only used by :func:`combine_intersections`.

    Along with :func:`get_next_first`, this function does the majority of the
    heavy lifting in :func:`get_next`. **Very** similar to
    :func:`get_next_first`, but this works with the second curve while the
    other function works with the first.

    Args:
        intersection (Intersection): The current intersection.
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of all detected intersections, provided as a reference for
            potential points to arrive at.
        to_end (Optional[bool]): Indicates if the next node should just be
            the end of the first edge or :data:`None`.

    Returns:
        point along a triangle of intersection. This will produce the next
        Optional[~bezier.hazmat.intersection_helpers.Intersection]: The "next"
        intersection along the current (second) edge or the end of the same
        edge. If ``to_end`` is :data:`False` and there are no other
        intersections along the current edge, will return :data:`None` (rather
        than the end of the same edge).
    """
    along_edge = None
    index_second = intersection.index_second
    t = intersection.t
    for other_int in intersections:
        other_t = other_int.t
        if other_int.index_second == index_second and other_t > t:
            if along_edge is None or other_t < along_edge.t:
                along_edge = other_int
    if along_edge is None:
        if to_end:
            # If there is no other intersection on the edge, just return
            # the segment end.
            return intersection_helpers.Intersection(
                None,
                None,
                index_second,
                1.0,
                interior_curve=CLASSIFICATION_T.SECOND,
            )

        else:
            return None

    else:
        return along_edge


def get_next_coincident(intersection, intersections):
    """Gets the next node along the current (coincident) edge.

    .. note::

       This is a helper used only by :func:`get_next`, which in
       turn is only used by :func:`basic_interior_combine`, which itself
       is only used by :func:`combine_intersections`.

    Along with :func:`get_next_first` and :func:`get_next_second`, this
    function does the majority of the heavy lifting in :func:`get_next`.

    This function moves immediately to the "next" intersection along the
    "current" edge. An intersection labeled as
    :attr:`~.IntersectionClassification.COINCIDENT` can only occur
    at the beginning of a segment. The end will correspond to the ``s`` or
    ``t`` parameter being equal to ``1.0`` and so it will get "rotated"
    to the front of the next edge, where it will be classified according
    to a different edge pair.

    It's also worth noting that some coincident segments will correspond
    to curves moving in opposite directions. In that case, there is
    no "interior" intersection (the curves are facing away) and so that
    segment won't be part of an intersection.

    Args:
        intersection (Intersection): The current intersection.
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of all detected intersections, provided as a reference for
            potential points to arrive at.

    Returns:
        Intersection: The "next" point along a triangle of intersection.
        This will produce the next intersection along the current (second)
        edge or the end of the same edge.
    """
    # Moving along the first or second edge will produce the same result, but
    # since we "rotate" all intersections to the beginning of an edge, the
    # index may not match the first or second (or both).
    along_first = get_next_first(intersection, intersections, to_end=False)
    if along_first is None:
        along_second = get_next_second(
            intersection, intersections, to_end=False
        )
    else:
        return along_first

    if along_second is None:
        return intersection_helpers.Intersection(
            intersection.index_first,
            1.0,
            intersection.index_second,
            1.0,
            interior_curve=CLASSIFICATION_T.COINCIDENT,
        )

    else:
        return along_second


def is_first(classification):
    """Check if a classification is on the "first" curve.

    Args:
        classification (IntersectionClassification): The classification
            being checked.

    Returns:
        bool: Indicating if the classification is on the first curve.
    """
    return classification in (
        CLASSIFICATION_T.FIRST,
        CLASSIFICATION_T.TANGENT_FIRST,
    )


def is_second(classification):
    """Check if a classification is on the "second" curve.

    Args:
        classification (IntersectionClassification): The classification
            being checked.

    Returns:
        bool: Indicating if the classification is on the second curve.
    """
    return classification in (
        CLASSIFICATION_T.SECOND,
        CLASSIFICATION_T.TANGENT_SECOND,
    )


def get_next(intersection, intersections, unused):
    """Gets the next node along a given edge.

    .. note::

       This is a helper used only by :func:`basic_interior_combine`, which in
       turn is only used by :func:`combine_intersections`. This function does
       the majority of the heavy lifting for :func:`basic_interior_combine`.

    .. note::

        This function returns :class:`.Intersection` objects even
        when the point isn't strictly an intersection. This is
        "incorrect" in some sense, but for now, we don't bother
        implementing a class similar to, but different from,
        :class:`.Intersection` to satisfy this need.

    Args:
        intersection (Intersection): The current intersection.
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of all detected intersections, provided as a reference for
            of potential points to arrive at.
        unused (List[~bezier.hazmat.intersection_helpers.Intersection]): List
            nodes that haven't been used yet in an intersection curved polygon

    Returns:
        Intersection: The "next" point along a triangle of intersection.
        This will produce the next intersection along the current edge or
        the end of the current edge.

    Raises:
        ValueError: If the intersection is not classified as
            :attr:`~.IntersectionClassification.FIRST`,
            :attr:`~.IntersectionClassification.TANGENT_FIRST`,
            :attr:`~.IntersectionClassification.SECOND`,
            :attr:`~.IntersectionClassification.TANGENT_SECOND` or
            :attr:`~.IntersectionClassification.COINCIDENT`.
    """
    result = None
    if is_first(intersection.interior_curve):
        result = get_next_first(intersection, intersections)
    elif is_second(intersection.interior_curve):
        result = get_next_second(intersection, intersections)
    elif intersection.interior_curve == CLASSIFICATION_T.COINCIDENT:
        result = get_next_coincident(intersection, intersections)
    else:
        raise ValueError(
            'Cannot get next node if not starting from "FIRST", '
            '"TANGENT_FIRST", "SECOND", "TANGENT_SECOND" or "COINCIDENT".'
        )

    if result in unused:
        unused.remove(result)
    return result


def ends_to_curve(start_node, end_node):
    """Convert a "pair" of intersection nodes to a curve segment.

    .. note::

       This is a helper used only by :func:`basic_interior_combine`, which in
       turn is only used by :func:`combine_intersections`.

    .. note::

       This function could specialize to the first or second segment
       attached to ``start_node`` and ``end_node``. We determine
       first / second based on the classification of ``start_node``,
       but the callers of this function could provide that information /
       isolate the base curve and the two parameters for us.

    .. note::

       This only checks the classification of the ``start_node``.

    Args:
        start_node (Intersection): The beginning of a segment.
        end_node (Intersection): The end of (the same) segment.

    Returns:
        Tuple[int, float, float]: The 3-tuple of:

        * The edge index along the first triangle (if in ``{0, 1, 2}``)
          or the edge index along the second triangle shifted to the right by
          3 (if in ``{3, 4, 5}``)
        * The start parameter along the edge
        * The end parameter along the edge

    Raises:
        ValueError: If the ``start_node`` and ``end_node`` disagree on
            the first curve when classified as
            :attr:`~.IntersectionClassification.FIRST`.
        ValueError: If the ``start_node`` and ``end_node`` disagree on
            the second curve when classified as
            :attr:`~.IntersectionClassification.SECOND`.
        ValueError: If the ``start_node`` and ``end_node`` disagree on
            the both curves when classified as
            :attr:`~.IntersectionClassification.COINCIDENT`.
        ValueError: If the ``start_node`` is not classified as
            :attr:`~.IntersectionClassification.FIRST`,
            :attr:`~.IntersectionClassification.TANGENT_FIRST`,
            :attr:`~.IntersectionClassification.SECOND`,
            :attr:`~.IntersectionClassification.TANGENT_SECOND` or
            :attr:`~.IntersectionClassification.COINCIDENT`.
    """
    if is_first(start_node.interior_curve):
        if end_node.index_first != start_node.index_first:
            raise ValueError(_WRONG_CURVE)

        return start_node.index_first, start_node.s, end_node.s

    elif is_second(start_node.interior_curve):
        if end_node.index_second != start_node.index_second:
            raise ValueError(_WRONG_CURVE)

        return start_node.index_second + 3, start_node.t, end_node.t

    elif start_node.interior_curve == CLASSIFICATION_T.COINCIDENT:
        if end_node.index_first == start_node.index_first:
            return start_node.index_first, start_node.s, end_node.s

        elif end_node.index_second == start_node.index_second:
            return start_node.index_second + 3, start_node.t, end_node.t

        else:
            raise ValueError(_WRONG_CURVE)

    else:
        raise ValueError(
            'Segment start must be classified as "FIRST", "TANGENT_FIRST", '
            '"SECOND", "TANGENT_SECOND" or "COINCIDENT".'
        )


def no_intersections(nodes1, degree1, nodes2, degree2):
    r"""Determine if one triangle is in the other.

    Helper for :func:`combine_intersections` that handles the case
    of no points of intersection. In this case, either the triangles
    are disjoint or one is fully contained in the other.

    To check containment, it's enough to check if one of the corners
    is contained in the other triangle.

    Args:
        nodes1 (numpy.ndarray): The nodes defining the first triangle in
            the intersection (assumed in :math:\mathbf{R}^2`).
        degree1 (int): The degree of the triangle given by ``nodes1``.
        nodes2 (numpy.ndarray): The nodes defining the second triangle in
            the intersection (assumed in :math:\mathbf{R}^2`).
        degree2 (int): The degree of the triangle given by ``nodes2``.

    Returns:
        Tuple[Optional[list], Optional[bool]]: Pair (2-tuple) of

        * Edges info list; will be empty or :data:`None`
        * "Contained" boolean. If not :data:`None`, indicates
          that one of the triangles is contained in the other.
    """
    # NOTE: This is a circular import.
    # pylint: disable=import-outside-toplevel
    from bezier.hazmat import triangle_intersection

    # pylint: enable=import-outside-toplevel

    located = triangle_intersection.locate_point(
        nodes2, degree2, nodes1[0, 0], nodes1[1, 0]
    )
    if located is not None:
        return None, True

    located = triangle_intersection.locate_point(
        nodes1, degree1, nodes2[0, 0], nodes2[1, 0]
    )
    if located is not None:
        return None, False

    return [], None


def tangent_only_intersections(all_types):
    """Determine intersection in the case of only-tangent intersections.

    If the only intersections are tangencies, then either the triangles
    are tangent but don't meet ("kissing" edges) or one triangle is
    internally tangent to the other.

    Thus we expect every intersection to be classified as
    :attr:`~.IntersectionClassification.TANGENT_FIRST`,
    :attr:`~.IntersectionClassification.TANGENT_SECOND`,
    :attr:`~.IntersectionClassification.OPPOSED`,
    :attr:`~.IntersectionClassification.IGNORED_CORNER` or
    :attr:`~.IntersectionClassification.COINCIDENT_UNUSED`.

    What's more, we expect all intersections to be classified the same for
    a given pairing.

    Args:
        all_types (Set[ \
            ~bezier.hazmat.intersection_helpers.IntersectionClassification]):
            The set of all intersection classifications encountered among the
            intersections for the given triangle-triangle pair.

    Returns:
        Tuple[Optional[list], Optional[bool]]: Pair (2-tuple) of

        * Edges info list; will be empty or :data:`None`
        * "Contained" boolean. If not :data:`None`, indicates
          that one of the triangles is contained in the other.

    Raises:
        ValueError: If there are intersections of more than one type among
            :attr:`~.IntersectionClassification.TANGENT_FIRST`,
            :attr:`~.IntersectionClassification.TANGENT_SECOND`,
            :attr:`~.IntersectionClassification.OPPOSED`,
            :attr:`~.IntersectionClassification.IGNORED_CORNER` or
            :attr:`~.IntersectionClassification.COINCIDENT_UNUSED`.
        ValueError: If there is a unique classification, but it isn't one
            of the tangent types.
    """
    if len(all_types) != 1:
        raise ValueError("Unexpected value, types should all match", all_types)

    point_type = all_types.pop()
    if point_type == CLASSIFICATION_T.OPPOSED:
        return [], None

    elif point_type == CLASSIFICATION_T.IGNORED_CORNER:
        return [], None

    elif point_type == CLASSIFICATION_T.TANGENT_FIRST:
        return None, True

    elif point_type == CLASSIFICATION_T.TANGENT_SECOND:
        return None, False

    elif point_type == CLASSIFICATION_T.COINCIDENT_UNUSED:
        return [], None

    else:
        raise ValueError("Point type not for tangency", point_type)


def basic_interior_combine(intersections, max_edges=10):
    """Combine intersections that don't involve tangencies.

    .. note::

       This is a helper used only by :func:`combine_intersections`.

    .. note::

       This helper assumes ``intersections`` isn't empty, but doesn't
       enforce it.

    Args:
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            Intersections from each of the 9 edge-edge pairs from a
            triangle-triangle pairing.
        max_edges (Optional[int]): The maximum number of allowed / expected
            edges per intersection. This is to avoid infinite loops.

    Returns:
        Tuple[Optional[list], Optional[bool]]: Pair (2-tuple) of

        * List of "edge info" lists. Each list represents a curved polygon
          and contains 3-tuples of edge index, start and end (see the
          output of :func:`ends_to_curve`).
        * "Contained" boolean. If not :data:`None`, indicates
          that one of the triangles is contained in the other.

    Raises:
        RuntimeError: If the number of edges in a curved polygon
            exceeds ``max_edges``. This is interpreted as a sign
            that the algorithm failed.
    """
    unused = intersections[:]
    result = []
    while unused:
        start = unused.pop()
        curr_node = start
        next_node = get_next(curr_node, intersections, unused)
        edge_ends = [(curr_node, next_node)]
        while next_node is not start:
            curr_node = to_front(next_node, intersections, unused)
            # NOTE: We also check to break when moving a corner node
            #       to the front. This is because ``intersections``
            #       de-duplicates corners by selecting the one
            #       (of 2 or 4 choices) at the front of segment(s).
            if curr_node is start:
                break

            next_node = get_next(curr_node, intersections, unused)
            edge_ends.append((curr_node, next_node))
            if len(edge_ends) > max_edges:
                raise RuntimeError(
                    "Unexpected number of edges", len(edge_ends)
                )

        edge_info = tuple(
            ends_to_curve(start_node, end_node)
            for start_node, end_node in edge_ends
        )
        result.append(edge_info)
    if len(result) == 1:
        if result[0] in FIRST_TRIANGLE_INFO:
            return None, True

        elif result[0] in SECOND_TRIANGLE_INFO:
            return None, False

    return result, None


def combine_intersections(
    intersections, nodes1, degree1, nodes2, degree2, all_types
):
    """Combine curve-curve intersections into curved polygon(s).

    .. note::

       This is a helper used only by :meth:`.Triangle.intersect`.

    Does so assuming each intersection lies on an edge of one of
    two :class:`.Triangle`-s.

    .. note ::

       This assumes that each ``intersection`` has been classified via
       :func:`classify_intersection` and only the intersections classified
       as :attr:`~.IntersectionClassification.FIRST` and
       :attr:`~.IntersectionClassification.SECOND` were kept.

    Args:
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            Intersections from each of the 9 edge-edge pairs from a
            triangle-triangle pairing.
        nodes1 (numpy.ndarray): The nodes defining the first triangle in
            the intersection (assumed in :math:`\\mathbf{R}^2`).
        degree1 (int): The degree of the triangle given by ``nodes1``.
        nodes2 (numpy.ndarray): The nodes defining the second triangle in
            the intersection (assumed in :math:`\\mathbf{R}^2`).
        degree2 (int): The degree of the triangle given by ``nodes2``.
        all_types (Set[ \
            ~bezier.hazmat.intersection_helpers.IntersectionClassification]):
            The set of all intersection classifications encountered among the
            intersections for the given triangle-triangle pair.

    Returns:
        Tuple[Optional[list], Optional[bool]]: Pair (2-tuple) of

        * List of "edge info" lists. Each list represents a curved polygon
          and contains 3-tuples of edge index, start and end (see the
          output of :func:`ends_to_curve`).
        * "Contained" boolean. If not :data:`None`, indicates
          that one of the triangles is contained in the other.
    """
    if intersections:
        return basic_interior_combine(intersections)

    elif all_types:
        return tangent_only_intersections(all_types)

    else:
        return no_intersections(nodes1, degree1, nodes2, degree2)


def evaluate_barycentric(nodes, degree, lambda1, lambda2, lambda3):
    r"""Compute a point on a triangle.

    Evaluates :math:`B\left(\lambda_1, \lambda_2, \lambda_3\right)` for a
    B |eacute| zier triangle / triangle defined by ``nodes``.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Control point nodes that define the triangle.
        degree (int): The degree of the triangle define by ``nodes``.
        lambda1 (float): Parameter along the reference triangle.
        lambda2 (float): Parameter along the reference triangle.
        lambda3 (float): Parameter along the reference triangle.

    Returns:
        numpy.ndarray: The evaluated point as a ``D x 1`` array (where ``D``
        is the ambient dimension where ``nodes`` reside).
    """
    dimension, num_nodes = nodes.shape
    binom_val = 1.0
    result = np.zeros((dimension, 1), order="F")
    index = num_nodes - 1
    result[:, 0] += nodes[:, index]
    # curve evaluate_multi_barycentric() takes arrays.
    lambda1 = np.asfortranarray([lambda1])
    lambda2 = np.asfortranarray([lambda2])
    for k in range(degree - 1, -1, -1):
        # We want to go from (d C (k + 1)) to (d C k).
        binom_val = (binom_val * (k + 1)) / (degree - k)
        index -= 1  # Step to last element in column.
        #     k = d - 1, d - 2, ...
        # d - k =     1,     2, ...
        # We know column k has (d - k + 1) elements.
        new_index = index - degree + k  # First element in column.
        col_nodes = nodes[:, new_index : index + 1]  # noqa: E203
        col_nodes = np.asfortranarray(col_nodes)
        col_result = curve_helpers.evaluate_multi_barycentric(
            col_nodes, lambda1, lambda2
        )
        result *= lambda3
        result += binom_val * col_result
        # Update index for next iteration.
        index = new_index
    return result


def evaluate_barycentric_multi(nodes, degree, param_vals, dimension):
    r"""Compute multiple points on the triangle.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Control point nodes that define the triangle.
        degree (int): The degree of the triangle define by ``nodes``.
        param_vals (numpy.ndarray): Array of parameter values (as a
            ``N x 3`` array).
        dimension (int): The dimension the triangle lives in.

    Returns:
        numpy.ndarray: The evaluated points, where columns correspond to
        rows of ``param_vals`` and the rows to the dimension of the
        underlying triangle.
    """
    num_vals, _ = param_vals.shape
    result = np.empty((dimension, num_vals), order="F")
    for index, (lambda1, lambda2, lambda3) in enumerate(param_vals):
        result[:, index] = evaluate_barycentric(
            nodes, degree, lambda1, lambda2, lambda3
        )[:, 0]
    return result


def evaluate_cartesian_multi(nodes, degree, param_vals, dimension):
    r"""Compute multiple points on the triangle.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Control point nodes that define the triangle.
        degree (int): The degree of the triangle define by ``nodes``.
        param_vals (numpy.ndarray): Array of parameter values (as a
            ``N x 2`` array).
        dimension (int): The dimension the triangle lives in.

    Returns:
        numpy.ndarray: The evaluated points, where columns correspond to
        rows of ``param_vals`` and the rows to the dimension of the
        underlying triangle.
    """
    num_vals, _ = param_vals.shape
    result = np.empty((dimension, num_vals), order="F")
    for index, (s, t) in enumerate(param_vals):
        result[:, index] = evaluate_barycentric(
            nodes, degree, 1.0 - s - t, s, t
        )[:, 0]
    return result


def compute_edge_nodes(nodes, degree):
    """Compute the nodes of each edges of a triangle.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Control point nodes that define the triangle.
        degree (int): The degree of the triangle define by ``nodes``.

    Returns:
        Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: The nodes in
        the edges of the triangle.
    """
    dimension, _ = np.shape(nodes)
    nodes1 = np.empty((dimension, degree + 1), order="F")
    nodes2 = np.empty((dimension, degree + 1), order="F")
    nodes3 = np.empty((dimension, degree + 1), order="F")
    curr2 = degree
    curr3 = -1
    for i in range(degree + 1):
        nodes1[:, i] = nodes[:, i]
        nodes2[:, i] = nodes[:, curr2]
        nodes3[:, i] = nodes[:, curr3]
        # Update the indices.
        curr2 += degree - i
        curr3 -= i + 2
    return nodes1, nodes2, nodes3


def shoelace_for_area(nodes):
    r"""Compute an auxiliary "shoelace" sum used to compute area.

    .. note::

       This is a helper for :func:`compute_area`.

    Defining :math:`\left[i, j\right] = x_i y_j - y_i x_j` as a shoelace
    term illuminates the name of this helper. On a degree one curve, this
    function will return

    .. math::

       \frac{1}{2}\left[0, 1\right].

    on a degree two curve it will return

    .. math::

       \frac{1}{6}\left(2 \left[0, 1\right] + 2 \left[1, 2\right] +
           \left[0, 2\right]\right)

    and so on.

    For a given :math:`\left[i, j\right]`, the coefficient comes from
    integrating :math:`b_{i, d}, b_{j, d}` on :math:`\left[0, 1\right]` (where
    :math:`b_{i, d}, b_{j, d}` are Bernstein basis polynomials).

    Returns:
        float: The computed sum of shoelace terms.

    Raises:
        UnsupportedDegree: If the degree is not 1, 2, 3 or 4.
    """
    _, num_nodes = nodes.shape
    if num_nodes == 2:
        shoelace = SHOELACE_LINEAR
        scale_factor = 2.0
    elif num_nodes == 3:
        shoelace = SHOELACE_QUADRATIC
        scale_factor = 6.0
    elif num_nodes == 4:
        shoelace = SHOELACE_CUBIC
        scale_factor = 20.0
    elif num_nodes == 5:
        shoelace = SHOELACE_QUARTIC
        scale_factor = 70.0
    else:
        raise _py_helpers.UnsupportedDegree(
            num_nodes - 1, supported=(1, 2, 3, 4)
        )

    result = 0.0
    for multiplier, index1, index2 in shoelace:
        result += multiplier * (
            nodes[0, index1] * nodes[1, index2]
            - nodes[1, index1] * nodes[0, index2]
        )

    return result / scale_factor


def compute_area(edges):
    """Compute the area of a curved polygon.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Uses Green's theorem to compute the area exactly. See
    :attr:`.Triangle.area` and :attr:`.CurvedPolygon.area` for more
    information.

    Args:
        edges (Tuple[numpy.ndarray]): A list of ``2 x N`` arrays, each of which
            are control points for an edge of a curved polygon.

    Returns:
        float: The computed area.
    """
    result = 0.0
    for edge_nodes in edges:
        result += shoelace_for_area(edge_nodes)

    return result
