# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Private helper methods for :mod:`bezier.surface`.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


import collections
import enum
import functools
import operator

import numpy as np
import six

from bezier import _curve_helpers
from bezier import _helpers
from bezier import _intersection_helpers
from bezier import curved_polygon
try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


_MAX_POLY_SUBDIVISIONS = 5
_MAX_LOCATE_SUBDIVISIONS = 20
_LOCATE_EPS = 2.0**(-47)
_SIGN = np.sign  # pylint: disable=no-member
_FLOAT64 = np.float64  # pylint: disable=no-member
_SAME_CURVATURE = 'Tangent curves have same curvature.'
_BAD_TANGENT = (
    'Curves moving in opposite direction but define '
    'overlapping arcs.')
_WRONG_CURVE = 'Start and end node not defined on same curve'
# NOTE: The ``SUBDIVIDE`` matrices are public since used in
#       the ``surface`` module.
LINEAR_SUBDIVIDE = np.array([
    [2, 0, 0],
    [1, 1, 0],
    [0, 2, 0],
    [1, 0, 1],
    [0, 1, 1],
    [0, 0, 2],
], dtype=_FLOAT64) / 2.0
QUADRATIC_SUBDIVIDE = np.array([
    [4, 0, 0, 0, 0, 0],
    [2, 2, 0, 0, 0, 0],
    [1, 2, 1, 0, 0, 0],
    [0, 2, 2, 0, 0, 0],
    [0, 0, 4, 0, 0, 0],
    [2, 0, 0, 2, 0, 0],
    [1, 1, 0, 1, 1, 0],
    [0, 1, 1, 1, 1, 0],
    [0, 0, 2, 0, 2, 0],
    [1, 0, 0, 2, 0, 1],
    [0, 1, 0, 1, 1, 1],
    [0, 0, 1, 0, 2, 1],
    [0, 0, 0, 2, 0, 2],
    [0, 0, 0, 0, 2, 2],
    [0, 0, 0, 0, 0, 4],
], dtype=_FLOAT64) / 4.0
CUBIC_SUBDIVIDE = np.array([
    [8, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 4, 0, 0, 0, 0, 0, 0, 0, 0],
    [2, 4, 2, 0, 0, 0, 0, 0, 0, 0],
    [1, 3, 3, 1, 0, 0, 0, 0, 0, 0],
    [0, 2, 4, 2, 0, 0, 0, 0, 0, 0],
    [0, 0, 4, 4, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 8, 0, 0, 0, 0, 0, 0],
    [4, 0, 0, 0, 4, 0, 0, 0, 0, 0],
    [2, 2, 0, 0, 2, 2, 0, 0, 0, 0],
    [1, 2, 1, 0, 1, 2, 1, 0, 0, 0],
    [0, 1, 2, 1, 1, 2, 1, 0, 0, 0],
    [0, 0, 2, 2, 0, 2, 2, 0, 0, 0],
    [0, 0, 0, 4, 0, 0, 4, 0, 0, 0],
    [2, 0, 0, 0, 4, 0, 0, 2, 0, 0],
    [1, 1, 0, 0, 2, 2, 0, 1, 1, 0],
    [0, 1, 1, 0, 1, 2, 1, 1, 1, 0],
    [0, 0, 1, 1, 0, 2, 2, 1, 1, 0],
    [0, 0, 0, 2, 0, 0, 4, 0, 2, 0],
    [1, 0, 0, 0, 3, 0, 0, 3, 0, 1],
    [0, 1, 0, 0, 1, 2, 0, 2, 1, 1],
    [0, 0, 1, 0, 0, 2, 1, 1, 2, 1],
    [0, 0, 0, 1, 0, 0, 3, 0, 3, 1],
    [0, 0, 0, 0, 2, 0, 0, 4, 0, 2],
    [0, 0, 0, 0, 0, 2, 0, 2, 2, 2],
    [0, 0, 0, 0, 0, 0, 2, 0, 4, 2],
    [0, 0, 0, 0, 0, 0, 0, 4, 0, 4],
    [0, 0, 0, 0, 0, 0, 0, 0, 4, 4],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 8],
], dtype=_FLOAT64) / 8.0
QUARTIC_SUBDIVIDE = np.array([
    [16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [2, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [1, 4, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 2, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 4, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [4, 4, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0],
    [2, 4, 2, 0, 0, 2, 4, 2, 0, 0, 0, 0, 0, 0, 0],
    [1, 3, 3, 1, 0, 1, 3, 3, 1, 0, 0, 0, 0, 0, 0],
    [0, 1, 3, 3, 1, 1, 3, 3, 1, 0, 0, 0, 0, 0, 0],
    [0, 0, 2, 4, 2, 0, 2, 4, 2, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 4, 4, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0],
    [4, 0, 0, 0, 0, 8, 0, 0, 0, 4, 0, 0, 0, 0, 0],
    [2, 2, 0, 0, 0, 4, 4, 0, 0, 2, 2, 0, 0, 0, 0],
    [1, 2, 1, 0, 0, 2, 4, 2, 0, 1, 2, 1, 0, 0, 0],
    [0, 1, 2, 1, 0, 1, 3, 3, 1, 1, 2, 1, 0, 0, 0],
    [0, 0, 1, 2, 1, 0, 2, 4, 2, 1, 2, 1, 0, 0, 0],
    [0, 0, 0, 2, 2, 0, 0, 4, 4, 0, 2, 2, 0, 0, 0],
    [0, 0, 0, 0, 4, 0, 0, 0, 8, 0, 0, 4, 0, 0, 0],
    [2, 0, 0, 0, 0, 6, 0, 0, 0, 6, 0, 0, 2, 0, 0],
    [1, 1, 0, 0, 0, 3, 3, 0, 0, 3, 3, 0, 1, 1, 0],
    [0, 1, 1, 0, 0, 1, 3, 2, 0, 2, 3, 1, 1, 1, 0],
    [0, 0, 1, 1, 0, 0, 2, 3, 1, 1, 3, 2, 1, 1, 0],
    [0, 0, 0, 1, 1, 0, 0, 3, 3, 0, 3, 3, 1, 1, 0],
    [0, 0, 0, 0, 2, 0, 0, 0, 6, 0, 0, 6, 0, 2, 0],
    [1, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0, 4, 0, 1],
    [0, 1, 0, 0, 0, 1, 3, 0, 0, 3, 3, 0, 3, 1, 1],
    [0, 0, 1, 0, 0, 0, 2, 2, 0, 1, 4, 1, 2, 2, 1],
    [0, 0, 0, 1, 0, 0, 0, 3, 1, 0, 3, 3, 1, 3, 1],
    [0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 6, 0, 4, 1],
    [0, 0, 0, 0, 0, 2, 0, 0, 0, 6, 0, 0, 6, 0, 2],
    [0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 4, 0, 4, 2, 2],
    [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 4, 2, 2, 4, 2],
    [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 6, 0, 6, 2],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8, 0, 4],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 4, 4],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 8, 4],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 8],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 16],
], dtype=_FLOAT64) / 16.0
# The Jacobian of a quadratric (in any dimension) as given by
# dB/ds = [-2L1, 2(L1 - L2), 2L2, -2L3, 2L3, 0] * nodes
# dB/dt = [-2L1, -2L2, 0, 2(L1 - L3), 2L2, 2L3] * nodes
# We evaluate this at each of the 6 points in the quadratic
# triangle and then stack them (2 rows * 6 = 12 rows)
# pylint: disable=bad-whitespace
_QUADRATIC_JACOBIAN_HELPER = np.array([
    [-2,  2, 0,  0, 0, 0],
    [-2,  0, 0,  2, 0, 0],
    [-1,  0, 1,  0, 0, 0],
    [-1, -1, 0,  1, 1, 0],
    [ 0, -2, 2,  0, 0, 0],  # noqa: E201
    [ 0, -2, 0,  0, 2, 0],  # noqa: E201
    [-1,  1, 0, -1, 1, 0],
    [-1,  0, 0,  0, 0, 1],
    [ 0, -1, 1, -1, 1, 0],  # noqa: E201
    [ 0, -1, 0, -1, 1, 1],  # noqa: E201
    [ 0,  0, 0, -2, 2, 0],  # noqa: E201
    [ 0,  0, 0, -2, 0, 2],  # noqa: E201
], dtype=_FLOAT64)
_QUADRATIC_TO_BERNSTEIN = np.array([
    [ 2, 0,  0, 0, 0,  0],  # noqa: E201
    [-1, 4, -1, 0, 0,  0],
    [ 0, 0,  2, 0, 0,  0],  # noqa: E201
    [-1, 0,  0, 4, 0, -1],
    [ 0, 0, -1, 0, 4, -1],  # noqa: E201
    [ 0, 0,  0, 0, 0,  2],  # noqa: E201
], dtype=_FLOAT64) / 2.0
# pylint: enable=bad-whitespace
# The Jacobian of a cubic (in any dimension) as given by
# dB/ds = [-3 L1^2, 3 L1(L1 - 2 L2), 3 L2(2 L1 - L2), 3 L2^2, -6 L1 L3,
#          6 L3(L1 - L2), 6 L2 L3, -3 L3^2, 3 L3^2, 0] * nodes
# dB/dt = [-3 L1^2, -6 L1 L2, -3 L2^2, 0, 3 L1(L1 - 2 L3), 6 L2 (L1 - L3),
#          3 L2^2, 3 L3(2 L1 - L3), 6 L2 L3, 3 L3^2] * nodes
# We evaluate this at each of the 15 points in the quartic
# triangle and then stack them (2 rows * 15 = 30 rows)
# pylint: disable=bad-whitespace
_CUBIC_JACOBIAN_HELPER = np.array([
    [-48,  48,   0,  0,   0,   0,  0,   0,  0,  0],
    [-48,   0,   0,  0,  48,   0,  0,   0,  0,  0],
    [-27,   9,  15,  3,   0,   0,  0,   0,  0,  0],
    [-27, -18,  -3,  0,  27,  18,  3,   0,  0,  0],
    [-12, -12,  12, 12,   0,   0,  0,   0,  0,  0],
    [-12, -24, -12,  0,  12,  24, 12,   0,  0,  0],
    [ -3, -15,  -9, 27,   0,   0,  0,   0,  0,  0],  # noqa: E201
    [ -3, -18, -27,  0,   3,  18, 27,   0,  0,  0],  # noqa: E201
    [  0,   0, -48, 48,   0,   0,  0,   0,  0,  0],  # noqa: E201
    [  0,   0, -48,  0,   0,   0, 48,   0,  0,  0],  # noqa: E201
    [-27,  27,   0,  0, -18,  18,  0,  -3,  3,  0],
    [-27,   0,   0,  0,   9,   0,  0,  15,  0,  3],
    [-12,   0,   9,  3, -12,   6,  6,  -3,  3,  0],
    [-12, -12,  -3,  0,   0,   6,  3,   9,  6,  3],
    [ -3,  -9,   0, 12,  -6,  -6, 12,  -3,  3,  0],  # noqa: E201
    [ -3, -12, -12,  0,  -3,   0, 12,   3, 12,  3],  # noqa: E201
    [  0,   0, -27, 27,   0, -18, 18,  -3,  3,  0],  # noqa: E201
    [  0,   0, -27,  0,   0, -18, 27,  -3, 18,  3],  # noqa: E201
    [-12,  12,   0,  0, -24,  24,  0, -12, 12,  0],
    [-12,   0,   0,  0, -12,   0,  0,  12,  0, 12],
    [ -3,  -3,   3,  3, -12,   0, 12, -12, 12,  0],  # noqa: E201
    [ -3,  -6,  -3,  0,  -9,  -6,  3,   0, 12, 12],  # noqa: E201
    [  0,   0, -12, 12,   0, -24, 24, -12, 12,  0],  # noqa: E201
    [  0,   0, -12,  0,   0, -24, 12, -12, 24, 12],  # noqa: E201
    [ -3,   3,   0,  0, -18,  18,  0, -27, 27,  0],  # noqa: E201
    [ -3,   0,   0,  0, -15,   0,  0,  -9,  0, 27],  # noqa: E201
    [  0,   0,  -3,  3,   0, -18, 18, -27, 27,  0],  # noqa: E201
    [  0,   0,  -3,  0,   0, -18,  3, -27, 18, 27],  # noqa: E201
    [  0,   0,   0,  0,   0,   0,  0, -48, 48,  0],  # noqa: E201
    [  0,   0,   0,  0,   0,   0,  0, -48,  0, 48],  # noqa: E201
], dtype=_FLOAT64) / 16.0
# pylint: enable=bad-whitespace
_QUARTIC_TO_BERNSTEIN = np.array([
    [36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [-39, 144, -108, 48, -9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [26, -128, 240, -128, 26, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [-9, 48, -108, 144, -39, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [0, 0, 0, 0, 36, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    [-39, 0, 0, 0, 0, 144, 0, 0, 0, -108, 0, 0, 48, 0, -9],
    [26, -64, -24, 32, -9, -64, 288, -96, 16, -24, -96, 12, 32, 16, -9],
    [-9, 32, -24, -64, 26, 16, -96, 288, -64, 12, -96, -24, 16, 32, -9],
    [0, 0, 0, 0, -39, 0, 0, 0, 144, 0, 0, -108, 0, 48, -9],
    [26, 0, 0, 0, 0, -128, 0, 0, 0, 240, 0, 0, -128, 0, 26],
    [-9, 16, 12, 16, -9, 32, -96, -96, 32, -24, 288, -24, -64, -64, 26],
    [0, 0, 0, 0, 26, 0, 0, 0, -128, 0, 0, 240, 0, -128, 26],
    [-9, 0, 0, 0, 0, 48, 0, 0, 0, -108, 0, 0, 144, 0, -39],
    [0, 0, 0, 0, -9, 0, 0, 0, 48, 0, 0, -108, 0, 144, -39],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 36],
], dtype=_FLOAT64)
# NOTE: We avoid round-off until after ``_QUARTIC_TO_BERNSTEIN``
#       has been applied.
_QUARTIC_BERNSTEIN_FACTOR = 36.0


def polynomial_sign(poly_surface):
    r"""Determine the "sign" of a polynomial on the reference triangle.

    Checks if a polynomial :math:`p(s, t)` is positive, negative
    or mixed sign on the reference triangle.

    Does this by utilizing the B |eacute| zier form of :math:`p`: it is a
    convex combination of the Bernstein basis (real numbers) hence
    if the Bernstein basis is all positive, the polynomial must be.

    If the values are mixed, then we can recursively subdivide
    until we are in a region where the coefficients are all one
    sign.

    Args:
        poly_surface (Surface): A polynomial on the reference triangle
            specified as a surface.

    Returns:
        int: The sign of the polynomial. Will be one of ``-1``, ``1``
        or ``0``. A value of ``0`` indicates a mixed sign or the
        zero polynomial.

    Raises:
        ValueError: If no conclusion is reached after the maximum
            number of subdivisions.
    """
    # The indices where the corner nodes in a surface are.
    corner_indices = (0, poly_surface._degree, -1)
    sub_polys = [poly_surface]
    signs = set()
    for _ in six.moves.xrange(_MAX_POLY_SUBDIVISIONS):
        undecided = []
        for poly in sub_polys:
            # Avoid an unnecessarily copying the nodes.
            nodes = poly._nodes
            # First add all the signs of the corner nodes.
            signs.update(_SIGN(nodes[corner_indices, 0]).astype(int))
            # Then check if the ``nodes`` are **uniformly** one sign.
            if np.all(nodes == 0.0):
                signs.add(0)
            elif np.all(nodes > 0.0):
                signs.add(1)
            elif np.all(nodes < 0.0):
                signs.add(-1)
            else:
                undecided.append(poly)

            if len(signs) > 1:
                return 0

        sub_polys = functools.reduce(
            operator.add, [poly.subdivide() for poly in undecided], ())
        if not sub_polys:
            break

    if len(sub_polys) == 0:
        # NOTE: We are guaranteed that ``len(signs) <= 1``.
        return signs.pop()
    else:
        raise ValueError(
            'Did not reach a conclusion after max subdivisions',
            _MAX_POLY_SUBDIVISIONS)


def _2x2_det(mat):
    r"""Compute the determinant of a 2x2 matrix.

    This is "needed" because :func:`numpy.linalg.det` rounds off
    answers when it doesn't need to. For example:

    .. doctest:: 2-by-2

       >>> import numpy as np
       >>> mat = np.array([
       ...     [-24.0, 3.0],
       ...     [-27.0, 0.0],
       ... ]) / 16.0
       >>> actual_det = -mat[0, 1] * mat[1, 0]
       >>> np_det = np.linalg.det(mat)
       >>> np.abs(actual_det - np_det) == np.spacing(actual_det)
       True

    Args:
        mat (numpy.ndarray): A 2x2 matrix.

    Returns:
        float: The determinant of ``mat``.
    """
    return mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]


def quadratic_jacobian_polynomial(nodes):
    r"""Compute the Jacobian determinant of a quadratic surface.

    Converts :math:`\det(J(s, t))` to a polynomial on the reference
    triangle and represents it as a surface object.

    .. note::

       This assumes that ``nodes`` is 6x2 but doesn't verify this.
       (However, the multiplication by ``_QUADRATIC_JACOBIAN_HELPER``
       would fail if ``nodes`` wasn't 6xN and then the ensuing
       determinants would fail if there weren't 2 columns.)

    Args:
        nodes (numpy.ndarray): A 6x2 array of nodes in a surface.

    Returns:
        numpy.ndarray: Coefficients in Bernstein basis.
    """
    # First evaluate the Jacobian at each of the 6 nodes.
    # pylint: disable=no-member
    jac_parts = _QUADRATIC_JACOBIAN_HELPER.dot(nodes)
    # pylint: enable=no-member
    jac_at_nodes = np.empty((6, 1))
    jac_at_nodes[0, 0] = _2x2_det(jac_parts[:2, :])
    jac_at_nodes[1, 0] = _2x2_det(jac_parts[2:4, :])
    jac_at_nodes[2, 0] = _2x2_det(jac_parts[4:6, :])
    jac_at_nodes[3, 0] = _2x2_det(jac_parts[6:8, :])
    jac_at_nodes[4, 0] = _2x2_det(jac_parts[8:10, :])
    jac_at_nodes[5, 0] = _2x2_det(jac_parts[10:, :])

    # Convert the nodal values to the Bernstein basis...
    bernstein = _QUADRATIC_TO_BERNSTEIN.dot(jac_at_nodes)
    return bernstein


def cubic_jacobian_polynomial(nodes):
    r"""Compute the Jacobian determinant of a cubic surface.

    Converts :math:`\det(J(s, t))` to a polynomial on the reference
    triangle and represents it as a surface object.

    .. note::

       This assumes that ``nodes`` is 10x2 but doesn't verify this.
       (However, the multiplication by ``_CUBIC_JACOBIAN_HELPER``
       would fail if ``nodes`` wasn't 10xN and then the ensuing
       determinants would fail if there weren't 2 columns.)

    Args:
        nodes (numpy.ndarray): A 10x2 array of nodes in a surface.

    Returns:
        numpy.ndarray: 15x1 array, coefficients in Bernstein basis.
    """
    # First evaluate the Jacobian at each of the 15 nodes
    # in the quartic triangle.
    jac_parts = _CUBIC_JACOBIAN_HELPER.dot(nodes)
    jac_at_nodes = np.empty((15, 1))
    jac_at_nodes[0, 0] = _2x2_det(jac_parts[:2, :])
    jac_at_nodes[1, 0] = _2x2_det(jac_parts[2:4, :])
    jac_at_nodes[2, 0] = _2x2_det(jac_parts[4:6, :])
    jac_at_nodes[3, 0] = _2x2_det(jac_parts[6:8, :])
    jac_at_nodes[4, 0] = _2x2_det(jac_parts[8:10, :])
    jac_at_nodes[5, 0] = _2x2_det(jac_parts[10:12, :])
    jac_at_nodes[6, 0] = _2x2_det(jac_parts[12:14, :])
    jac_at_nodes[7, 0] = _2x2_det(jac_parts[14:16, :])
    jac_at_nodes[8, 0] = _2x2_det(jac_parts[16:18, :])
    jac_at_nodes[9, 0] = _2x2_det(jac_parts[18:20, :])
    jac_at_nodes[10, 0] = _2x2_det(jac_parts[20:22, :])
    jac_at_nodes[11, 0] = _2x2_det(jac_parts[22:24, :])
    jac_at_nodes[12, 0] = _2x2_det(jac_parts[24:26, :])
    jac_at_nodes[13, 0] = _2x2_det(jac_parts[26:28, :])
    jac_at_nodes[14, 0] = _2x2_det(jac_parts[28:, :])

    # Convert the nodal values to the Bernstein basis...
    # pylint: disable=no-member
    bernstein = _QUARTIC_TO_BERNSTEIN.dot(jac_at_nodes)
    # pylint: enable=no-member
    bernstein /= _QUARTIC_BERNSTEIN_FACTOR
    return bernstein


def _de_casteljau_one_round(nodes, degree, lambda1, lambda2, lambda3):
    r"""Performs one "round" of the de Casteljau algorithm for surfaces.

    Converts the ``nodes`` into a basis for a surface one degree smaller
    by using the barycentric weights:

    .. math::

       q_{i, j, k} = \lambda_1 \cdot p_{i + 1, j, k} +
           \lambda_2 \cdot p_{i, j + 1, k} + \lambda_2 \cdot p_{i, j, k + 1}

    .. note:

       For degree :math:`d`d, the number of nodes should be
       :math:`(d + 1)(d + 2)/2`, but we don't verify this property.

    Args:
        nodes (numpy.ndarray): The nodes to reduce.
        degree (int): The degree of the surface.
        lambda1 (float): Parameter along the reference triangle.
        lambda2 (float): Parameter along the reference triangle.
        lambda3 (float): Parameter along the reference triangle.

    Returns:
        numpy.ndarray: The converted nodes.
    """
    num_nodes, dimension = nodes.shape
    num_new_nodes = num_nodes - degree - 1

    new_nodes = np.empty((num_new_nodes, dimension))

    index = 0
    # parent_i1 = index + k
    # parent_i2 = index + k + 1
    # parent_i3 = index + degree + 1
    parent_i1 = 0
    parent_i2 = 1
    parent_i3 = degree + 1
    for k in six.moves.xrange(degree):
        for unused_j in six.moves.xrange(degree - k):
            # NOTE: i = (degree - 1) - j - k
            new_nodes[index, :] = (
                lambda1 * nodes[parent_i1, :] +
                lambda2 * nodes[parent_i2, :] +
                lambda3 * nodes[parent_i3, :])
            # Update all the indices.
            parent_i1 += 1
            parent_i2 += 1
            parent_i3 += 1
            index += 1

        # Update the indices that depend on k.
        parent_i1 += 1
        parent_i2 += 1

    return new_nodes


def _make_transform(degree, weights_a, weights_b, weights_c):
    """Compute matrices corresponding to the de Casteljau algorithm.

    Applies the de Casteljau to the identity matrix, thus
    effectively caching the algorithm in a transformation matrix.

    .. note::

       This is premature optimization. It's unclear if the time
       saved from "caching" one round of de Casteljau is cancelled
       out by the extra storage required for the 3 matrices.

    Args:
        degree (int): The degree of a candidate surface.
        weights_a (Tuple[float, float, float]): Barycentric weights for a
            point in the reference triangle
        weights_b (Tuple[float, float, float]): Barycentric weights for a
            point in the reference triangle
        weights_c (Tuple[float, float, float]): Barycentric weights for a
            point in the reference triangle

    Returns:
        Mapping[int, numpy.ndarray]: Mapping from keys to the de Casteljau
        transformation mappings. The keys are ``0`` corresponding to
        ``weights_a``, ``1`` to ``weights_b`` and ``2`` to ``weights_c``.
    """
    num_nodes = ((degree + 1) * (degree + 2)) // 2
    id_mat = np.eye(num_nodes)

    # Pre-compute the matrices that do the reduction so we don't
    # have to **actually** perform the de Casteljau algorithm
    # every time.
    transform = {
        0: de_casteljau_one_round(id_mat, degree, *weights_a),
        1: de_casteljau_one_round(id_mat, degree, *weights_b),
        2: de_casteljau_one_round(id_mat, degree, *weights_c),
    }
    return transform


def _reduced_to_matrix(shape, degree, vals_by_weight):
    r"""Converts a reduced values dictionary into a matrix.

    The ``vals_by_weight`` mapping has keys of the form:
    ``(0, ..., 1, ..., 2, ...)`` where the ``0`` corresponds
    to the number of times the first set of barycentric
    weights was used in the reduction process, and similarly
    for ``1`` and ``2``.

    These points correspond to barycentric weights in their
    own right. For example ``(0, 0, 1, 2)`` corresponds to
    the Barycentric weight
    :math:`\left(\frac{2}{4}, \frac{1}{4}, \frac{1}{4}\right)`.

    Once the keys in ``vals_by_weight`` have been converted
    to Barycentric coordinates, we order them according to
    our rule (bottom to top, left to right) and then return
    them in a single matrix.

    Args:
        shape (tuple): The shape of the result matrix.
        degree (int): The degree of the surface.
        vals_by_weight (Mapping[tuple, numpy.ndarray]): Dictionary
            of reduced nodes according to blending of each of the
            three sets of weights in a reduction.

    Returns:
        numpy.ndarray: The newly created reduced control points.
    """
    result = np.empty(shape)
    index = 0
    for k in six.moves.xrange(degree + 1):
        for j in six.moves.xrange(degree + 1 - k):
            i = degree - j - k
            key = (0,) * i + (1,) * j + (2,) * k
            result[index, :] = vals_by_weight[key]
            index += 1

    return result


def specialize_surface(nodes, degree, weights_a, weights_b, weights_c):
    """Specialize a surface to a reparameterization

    Does so by taking three points (in Barycentric form) within the
    reference triangle and then reparameterizing the surface onto
    the triangle formed by those three points.

    .. note::

       This assumes the surface is degree 1 or greater but doesn't check.

    Args:
        nodes (numpy.ndarray): Control points for a surface.
        degree (int): The degree of the surface.
        weights_a (Tuple[float, float, float]): Barycentric weights for a
            point in the reference triangle
        weights_b (Tuple[float, float, float]): Barycentric weights for a
            point in the reference triangle
        weights_c (Tuple[float, float, float]): Barycentric weights for a
            point in the reference triangle

    Returns:
        numpy.ndarray: The control points for the specialized surface.
    """
    # Uses A-->0, B-->1, C-->2 to represent the specialization used.
    partial_vals = {
        (0,): de_casteljau_one_round(nodes, degree, *weights_a),
        (1,): de_casteljau_one_round(nodes, degree, *weights_b),
        (2,): de_casteljau_one_round(nodes, degree, *weights_c),
    }

    for reduced_deg in six.moves.xrange(degree - 1, 0, -1):
        new_partial = {}
        transform = _make_transform(
            reduced_deg, weights_a, weights_b, weights_c)
        for key, sub_nodes in six.iteritems(partial_vals):
            # Our keys are ascending so we increment from the last value.
            for next_id in six.moves.xrange(key[-1], 2 + 1):
                new_key = key + (next_id,)
                new_partial[new_key] = transform[next_id].dot(sub_nodes)

        partial_vals = new_partial

    return _reduced_to_matrix(nodes.shape, degree, partial_vals)


def _mean_centroid(candidates):
    """Take the mean of all centroids in set of reference triangles.

    Args:
        candidates (List[.Surface]): List of surfaces. We'll only use
            the base ``x`` and ``y`` and the width of the reference
            triangle that each surface represents.

    Returns:
        Tuple[float, float]: The mean of all centroids.
    """
    sum_x = 0.0
    sum_y = 0.0
    sum_width = 0.0
    for candidate in candidates:
        # pylint: disable=protected-access
        sum_x += candidate._base_x
        sum_y += candidate._base_y
        sum_width += candidate._width
        # pylint: enable=protected-access

    denom = float(len(candidates))
    mean_x = sum_x / denom
    mean_y = sum_y / denom
    width_delta = sum_width / (3.0 * denom)

    return mean_x + width_delta, mean_y + width_delta


def jacobian_s(nodes, degree, dimension):
    r"""Compute :math:`\frac{\partial B}{\partial s}`.

    Args:
        nodes (numpy.ndarray): Array of nodes in a surface.
        degree (int): The degree of the surface.
        dimension (int): The dimension the surface lives in.

    Returns:
        numpy.ndarray: Nodes of the Jacobian surface in
            B |eacute| zier form.
    """
    num_nodes = (degree * (degree + 1)) // 2
    result = np.empty((num_nodes, dimension))

    index = 0
    i = 0
    for num_vals in six.moves.xrange(degree, 0, -1):
        for _ in six.moves.xrange(num_vals):
            result[index, :] = nodes[i + 1, :] - nodes[i, :]
            # Update the indices
            index += 1
            i += 1

        # In between each row, the index gains an extra value.
        i += 1

    return float(degree) * result


def jacobian_t(nodes, degree, dimension):
    r"""Compute :math:`\frac{\partial B}{\partial t}`.

    Args:
        nodes (numpy.ndarray): Array of nodes in a surface.
        degree (int): The degree of the surface.
        dimension (int): The dimension the surface lives in.

    Returns:
        numpy.ndarray: Nodes of the Jacobian surface in
            B |eacute| zier form.
    """
    num_nodes = (degree * (degree + 1)) // 2
    result = np.empty((num_nodes, dimension))

    index = 0
    i = 0
    j = degree + 1
    for num_vals in six.moves.xrange(degree, 0, -1):
        for _ in six.moves.xrange(num_vals):
            result[index, :] = nodes[j, :] - nodes[i, :]
            # Update the indices
            index += 1
            i += 1
            j += 1

        # In between each row, the index gains an extra value.
        i += 1

    return float(degree) * result


def newton_refine(surface, x_val, y_val, s, t):
    r"""Refine a solution to :math:`B(s, t) = p` using Newton's method.

    Computes updates via

    .. math::

       \left[\begin{array}{c}
           0 \\ 0 \end{array}\right] \approx
           \left(B\left(s_{\ast}, t_{\ast}\right) -
           \left[\begin{array}{c} x \\ y \end{array}\right]\right) +
           \left[\begin{array}{c c}
               B_s\left(s_{\ast}, t_{\ast}\right) &
               B_t\left(s_{\ast}, t_{\ast}\right) \end{array}\right]
           \left[\begin{array}{c}
               \Delta s \\ \Delta t \end{array}\right]

    For example, (with weights
    :math:`\lambda_1 = 1 - s - t, \lambda_2 = s, \lambda_3 = t`)
    consider the surface

    .. math::

       B(s, t) =
           \left[\begin{array}{c} 0 \\ 0 \end{array}\right] \lambda_1^2 +
           \left[\begin{array}{c} 1 \\ 0 \end{array}\right]
               2 \lambda_1 \lambda_2 +
           \left[\begin{array}{c} 2 \\ 0 \end{array}\right] \lambda_2^2 +
           \left[\begin{array}{c} 2 \\ 1 \end{array}\right]
               2 \lambda_1 \lambda_3 +
           \left[\begin{array}{c} 2 \\ 2 \end{array}\right]
               2 \lambda_2 \lambda_1 +
           \left[\begin{array}{c} 0 \\ 2 \end{array}\right] \lambda_3^2

    and the point
    :math:`B\left(\frac{1}{4}, \frac{1}{2}\right) =
    \frac{1}{4} \left[\begin{array}{c} 5 \\ 5 \end{array}\right]`.

    Starting from the **wrong** point
    :math:`s = \frac{1}{2}, t = \frac{1}{4}`, we have

    .. math::

       \begin{align*}
       \left[\begin{array}{c} x \\ y \end{array}\right] -
          B\left(\frac{1}{2}, \frac{1}{4}\right) &= \frac{1}{4}
          \left[\begin{array}{c} -1 \\ 2 \end{array}\right] \\
       DB\left(\frac{1}{2}, \frac{1}{4}\right) &= \frac{1}{2}
           \left[\begin{array}{c c} 3 & 2 \\ 1 & 6 \end{array}\right] \\
       \Longrightarrow \left[\begin{array}{c}
           \Delta s \\ \Delta t \end{array}\right] &= \frac{1}{32}
           \left[\begin{array}{c} -10 \\ 7 \end{array}\right]
       \end{align*}

    .. image:: images/newton_refine_surface.png
       :align: center

    .. testsetup:: newton-refine-surface

       import numpy as np
       import bezier
       from bezier._surface_helpers import newton_refine

    .. doctest:: newton-refine-surface

       >>> nodes = np.array([
       ...     [0.0, 0.0],
       ...     [1.0, 0.0],
       ...     [2.0, 0.0],
       ...     [2.0, 1.0],
       ...     [2.0, 2.0],
       ...     [0.0, 2.0],
       ... ])
       >>> surface = bezier.Surface(nodes, degree=2)
       >>> surface.is_valid
       True
       >>> (x_val, y_val), = surface.evaluate_cartesian(0.25, 0.5)
       >>> x_val, y_val
       (1.25, 1.25)
       >>> s, t = 0.5, 0.25
       >>> new_s, new_t = newton_refine(surface, x_val, y_val, s, t)
       >>> 32 * (new_s - s)
       -10.0
       >>> 32 * (new_t - t)
       7.0

    .. testcleanup:: newton-refine-surface

       import make_images
       make_images.newton_refine_surface(
           surface, x_val, y_val, s, t, new_s, new_t)

    Args:
        surface (.Surface): A B |eacute| zier surface (assumed to
            be two-dimensional).
        x_val (float): The :math:`x`-coordinate of a point
            on the surface.
        y_val (float): The :math:`y`-coordinate of a point
            on the surface.
        s (float): Approximate :math:`s`-value to be refined.
        t (float): Approximate :math:`t`-value to be refined.

    Returns:
        Tuple[float, float]: The refined :math:`s` and :math:`t` values.
    """
    (surf_x, surf_y), = surface.evaluate_cartesian(s, t, _verify=False)
    if surf_x == x_val and surf_y == y_val:
        # No refinement is needed.
        return s, t

    # Compute Jacobian nodes / stack them horizatonally.
    degree = surface._degree
    jac_both = np.hstack([
        jacobian_s(surface._nodes, degree, surface._dimension),
        jacobian_t(surface._nodes, degree, surface._dimension),
    ])

    lambda1 = 1.0 - s - t
    # The degree of the jacobian is one less.
    for reduced_deg in six.moves.xrange(degree - 1, 0, -1):
        jac_both = de_casteljau_one_round(
            jac_both, reduced_deg, lambda1, s, t)

    # The first column of the jacobian matrix is B_s (i.e. the
    # left-most values in ``jac_both``).
    jac_mat = np.array([
        [jac_both[0, 0], jac_both[0, 2]],
        [jac_both[0, 1], jac_both[0, 3]],
    ])
    rhs = np.array([
        [x_val - surf_x],
        [y_val - surf_y],
    ])
    soln = np.linalg.solve(jac_mat, rhs)
    return s + soln[0, 0], t + soln[1, 0]


def locate_point(surface, x_val, y_val):
    r"""Locate a point on a surface.

    Does so by recursively subdividing the surface and rejecting
    sub-surfaces with bounding boxes that don't contain the point.
    After the sub-surfaces are sufficiently small, uses Newton's
    method to zoom in on the intersection.

    Args:
        surface (.Surface): A B |eacute| zier surface (assumed to
            be two-dimensional).
        x_val (float): The :math:`x`-coordinate of a point
            on the surface.
        y_val (float): The :math:`y`-coordinate of a point
            on the surface.

    Returns:
        Optional[Tuple[float, float]]: The :math:`s` and :math:`t`
        values corresponding to ``x_val`` and ``y_val`` or
        :data:`None` if the point is not on the ``surface``.
    """
    candidates = [surface]
    for _ in six.moves.xrange(_MAX_LOCATE_SUBDIVISIONS + 1):
        next_candidates = []
        for candidate in candidates:
            nodes = candidate._nodes
            if _helpers.contains(nodes, x_val, y_val):
                next_candidates.extend(candidate.subdivide())

        candidates = next_candidates

    if not candidates:
        return None

    # We take the average of all centroids from the candidates
    # that may contain the point.
    s_approx, t_approx = _mean_centroid(candidates)
    s, t = newton_refine(surface, x_val, y_val, s_approx, t_approx)

    actual = surface.evaluate_cartesian(s, t, _verify=False)
    expected = np.array([[x_val, y_val]])
    if not _helpers.vector_close(actual, expected, eps=_LOCATE_EPS):
        s, t = newton_refine(surface, x_val, y_val, s, t)
    return s, t


def classify_intersection(intersection):
    r"""Determine which curve is on the "inside of the intersection".

    This is intended to be a helper for forming a :class:`.CurvedPolygon`
    from the edge intersections of two :class:`.Surface`-s. In order
    to move from one intersection to another (or to the end of an edge),
    the interior edge must be determined at the point of intersection.

    The "typical" case is on the interior of both edges:

    .. image:: images/classify_intersection1.png
       :align: center

    .. testsetup:: classify-intersection1, classify-intersection2,
                   classify-intersection3, classify-intersection4,
                   classify-intersection5, classify-intersection6,
                   classify-intersection7, classify-intersection8

       import numpy as np
       import bezier
       from bezier import _curve_helpers
       from bezier._intersection_helpers import Intersection
       from bezier._surface_helpers import classify_intersection

       def hodograph(curve, s):
           return _curve_helpers.evaluate_hodograph(
               curve._nodes, curve._degree, s)

       def curvature(curve, s):
           nodes = curve._nodes
           tangent = _curve_helpers.evaluate_hodograph(
               nodes, curve._degree, s)
           return _curve_helpers.get_curvature(
               nodes, curve._degree, tangent, s)

    .. doctest:: classify-intersection1
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.array([
       ...     [1.0 , 0.0 ],
       ...     [1.75, 0.25],
       ...     [2.0 , 1.0 ],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.array([
       ...     [0.0   , 0.0   ],
       ...     [1.6875, 0.0625],
       ...     [2.0   , 0.5   ],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.25, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> tangent1 = hodograph(curve1, s)
       >>> tangent1
       array([[ 1.25, 0.75]])
       >>> tangent2 = hodograph(curve2, t)
       >>> tangent2
       array([[ 2. , 0.5]])
       >>> intersection = Intersection(curve1, s, curve2, t)
       >>> classify_intersection(intersection)
       <IntersectionClassification.first: 'first'>

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

    .. image:: images/classify_intersection2.png
       :align: center

    .. doctest:: classify-intersection2
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.array([
       ...     [1.0, 0.0],
       ...     [1.5, 1.0],
       ...     [2.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.array([
       ...     [0.0, 0.0],
       ...     [1.5, 1.0],
       ...     [3.0, 0.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> intersection = Intersection(curve1, s, curve2, t)
       >>> classify_intersection(intersection)
       <IntersectionClassification.tangent_second: 'tangent_second'>

    .. testcleanup:: classify-intersection2

       import make_images
       make_images.classify_intersection2(s, curve1, curve2)

    Depending on the direction of the parameterizations, the interior
    curve may change, but we can use the (signed) `curvature`_ of each
    curve at that point to determine which is on the interior:

    .. _curvature: https://en.wikipedia.org/wiki/Curvature

    .. image:: images/classify_intersection3.png
       :align: center

    .. doctest:: classify-intersection3
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.array([
       ...     [2.0, 0.0],
       ...     [1.5, 1.0],
       ...     [1.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.array([
       ...     [3.0, 0.0],
       ...     [1.5, 1.0],
       ...     [0.0, 0.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> intersection = Intersection(curve1, s, curve2, t)
       >>> classify_intersection(intersection)
       <IntersectionClassification.tangent_first: 'tangent_first'>

    .. testcleanup:: classify-intersection3

       import make_images
       make_images.classify_intersection3(s, curve1, curve2)

    When the curves are moving in opposite directions at a point
    of tangency, there is no side to choose. Either the point of tangency
    is not part of any :class:`.CurvedPolygon` intersection

    .. image:: images/classify_intersection4.png
       :align: center

    .. doctest:: classify-intersection4
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.array([
       ...     [2.0, 0.0],
       ...     [1.5, 1.0],
       ...     [1.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.array([
       ...     [0.0, 0.0],
       ...     [1.5, 1.0],
       ...     [3.0, 0.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> intersection = Intersection(curve1, s, curve2, t)
       >>> classify_intersection(intersection)
       <IntersectionClassification.opposed: 'opposed'>

    .. testcleanup:: classify-intersection4

       import make_images
       make_images.classify_intersection4(s, curve1, curve2)

    or the point of tangency is a "degenerate" part of two
    :class:`.CurvedPolygon` intersections. It is "degenerate"
    because from one direction, the point should be classified
    as ``1`` and from another as ``0``:

    .. image:: images/classify_intersection5.png
       :align: center

    .. doctest:: classify-intersection5
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.array([
       ...     [1.0, 0.0],
       ...     [1.5, 1.0],
       ...     [2.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.array([
       ...     [3.0, 0.0],
       ...     [1.5, 1.0],
       ...     [0.0, 0.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> intersection = Intersection(curve1, s, curve2, t)
       >>> classify_intersection(intersection)
       Traceback (most recent call last):
         ...
       NotImplementedError: Curves moving in opposite direction
                            but define overlapping arcs.

    .. testcleanup:: classify-intersection5

       import make_images
       make_images.classify_intersection5(s, curve1, curve2)

    However, if the `curvature`_ of each curve is identical, we
    don't try to distinguish further:

    .. image:: images/classify_intersection6.png
       :align: center

    .. doctest:: classify-intersection6
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.array([
       ...     [ 0.375,  0.0625],
       ...     [-0.125, -0.0625],
       ...     [-0.125,  0.0625],
       ... ])
       >>> curve1 = bezier.Curve(nodes1, degree=2)
       >>> nodes2 = np.array([
       ...     [ 0.75,  0.25],
       ...     [-0.25, -0.25],
       ...     [-0.25,  0.25],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> hodograph(curve1, s)
       array([[-0.5, 0. ]])
       >>> hodograph(curve2, t)
       array([[-1., 0.]])
       >>> curvature(curve1, s)
       -2.0
       >>> curvature(curve2, t)
       -2.0
       >>> intersection = Intersection(curve1, s, curve2, t)
       >>> classify_intersection(intersection)
       Traceback (most recent call last):
         ...
       NotImplementedError: Tangent curves have same curvature.

    .. testcleanup:: classify-intersection6

       import make_images
       make_images.classify_intersection6(s, curve1, curve2)

    In addition to points of tangency, intersections that happen at
    the end of an edge need special handling:

    .. image:: images/classify_intersection7.png
       :align: center

    .. doctest:: classify-intersection7
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1a = np.array([
       ...     [0.0, 0.0 ],
       ...     [4.5, 0.0 ],
       ...     [9.0, 2.25],
       ... ])
       >>> curve1a = bezier.Curve(nodes1a, degree=2)
       >>> nodes2 = np.array([
       ...     [11.25, 0.0],
       ...     [ 9.0 , 4.5],
       ...     [ 2.75, 1.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2, degree=2)
       >>> s, t = 1.0, 0.375
       >>> curve1a.evaluate(s) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> intersection = Intersection(curve1a, s, curve2, t)
       >>> classify_intersection(intersection)
       Traceback (most recent call last):
         ...
       ValueError: ('Intersection occurs at the end of an edge',
                    's', 1.0, 't', 0.375)
       >>>
       >>> nodes1b = np.array([
       ...     [9.0, 2.25 ],
       ...     [4.5, 2.375],
       ...     [0.0, 2.5  ],
       ... ])
       >>> curve1b = bezier.Curve(nodes1b, degree=2)
       >>> curve1b.evaluate(0.0) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> intersection = Intersection(curve1b, 0.0, curve2, t)
       >>> classify_intersection(intersection)
       <IntersectionClassification.first: 'first'>

    .. testcleanup:: classify-intersection7

       import make_images
       make_images.classify_intersection7(s, curve1a, curve1b, curve2)

    As above, some intersections at the end of an edge are part of
    an actual intersection. However, some surfaces may just "kiss" at a
    corner intersection:

    .. image:: images/classify_intersection8.png
       :align: center

    .. doctest:: classify-intersection8
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1 = np.array([
       ...     [0.25 , 1.0  ],
       ...     [0.0  , 0.5  ],
       ...     [0.0  , 0.0  ],
       ...     [0.625, 0.875],
       ...     [0.5  , 0.375],
       ...     [1.0  , 0.75 ],
       ... ])
       >>> surface1 = bezier.Surface(nodes1, degree=2)
       >>> nodes2 = np.array([
       ...     [ 0.0625, 0.5  ],
       ...     [-0.25  , 1.0  ],
       ...     [-1.0   , 1.0  ],
       ...     [-0.5   , 0.125],
       ...     [-1.0   , 0.5  ],
       ...     [-1.0   , 0.0  ],
       ... ])
       >>> surface2 = bezier.Surface(nodes2, degree=2)
       >>> curve1, _, _ = surface1.edges
       >>> curve2, _, _ = surface2.edges
       >>> s, t = 0.5, 0.0
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([[ True, True]], dtype=bool)
       >>> intersection = Intersection(curve1, s, curve2, t)
       >>> classify_intersection(intersection)
       <IntersectionClassification.ignored_corner: 'ignored_corner'>

    .. testcleanup:: classify-intersection8

       import make_images
       make_images.classify_intersection8(
           s, curve1, surface1, curve2, surface2)

    .. note::

       This assumes the intersection occurs in :math:`\mathbf{R}^2`
       but doesn't check this.

    .. note::

       This function doesn't allow wiggle room / round-off when checking
       endpoints, nor when checking if the cross product is near zero,
       nor when curvatures are compared. However, the most "correct"
       version of this function likely should allow for some round off.

    Args:
        intersection (.Intersection): An intersection object.

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
        raise ValueError('Intersection occurs at the end of an edge',
                         's', intersection.s, 't', intersection.t)

    tangent1 = _curve_helpers.evaluate_hodograph(
        intersection.first._nodes, intersection.first._degree, intersection.s)
    tangent2 = _curve_helpers.evaluate_hodograph(
        intersection.second._nodes, intersection.second._degree,
        intersection.t)

    if _ignored_corner(intersection, tangent1, tangent2):
        return IntersectionClassification.ignored_corner

    # Take the cross product of tangent vectors to determine which one
    # is more "inside" / "to the left".
    cross_prod = _helpers.cross_product(tangent1, tangent2)
    if cross_prod < 0:
        return IntersectionClassification.first
    elif cross_prod > 0:
        return IntersectionClassification.second
    else:
        return _classify_tangent_intersection(
            intersection, tangent1, tangent2)


def _classify_tangent_intersection(intersection, tangent1, tangent2):
    """Helper for func:`classify_intersection` at tangencies.

    Args:
        intersection (.Intersection): An intersection object.
        tangent1 (numpy.ndarray): The tangent vector to the first curve
            at the intersection.
        tangent2 (numpy.ndarray): The tangent vector to the second curve
            at the intersection.

    Returns:
        IntersectionClassification: The "inside" curve type, based on
        the classification enum. Will either be ``opposed`` or one
        of the ``tangent`` values.

    Raises:
        NotImplementedError: If the curves are tangent, moving in opposite
            directions, but enclose overlapping arcs.
        NotImplementedError: If the curves are tangent at the intersection
            and have the same curvature.
    """
    # Each array is 1x2 (i.e. a row vector).
    dot_prod = tangent1.dot(tangent2.T)
    # Unpack 1x1 array into a scalar (and assert size).
    (dot_prod,), = dot_prod
    # NOTE: When computing curvatures we assume that we don't have lines
    #       here, because lines that are tangent at an intersection are
    #       parallel and we don't handle that case.
    curvature1 = _curve_helpers.get_curvature(
        intersection.first._nodes, intersection.first._degree,
        tangent1, intersection.s)
    curvature2 = _curve_helpers.get_curvature(
        intersection.second._nodes, intersection.second._degree,
        tangent2, intersection.t)
    if dot_prod < 0:
        # If the tangent vectors are pointing in the opposite direction,
        # then the curves are facing opposite directions.
        sign1, sign2 = _SIGN([curvature1, curvature2])
        if sign1 == sign2:
            # If both curvatures are positive, since the curves are
            # moving in opposite directions, the tangency isn't part of
            # the surface intersection.
            if sign1 == 1.0:
                return IntersectionClassification.opposed
            else:
                raise NotImplementedError(_BAD_TANGENT)
        else:
            delta_c = abs(curvature1) - abs(curvature2)
            if delta_c == 0.0:
                raise NotImplementedError(_SAME_CURVATURE)
            elif sign1 == _SIGN(delta_c):
                return IntersectionClassification.opposed
            else:
                raise NotImplementedError(_BAD_TANGENT)
    else:
        if curvature1 > curvature2:
            return IntersectionClassification.tangent_first
        elif curvature1 < curvature2:
            return IntersectionClassification.tangent_second
        else:
            raise NotImplementedError(_SAME_CURVATURE)


def edge_cycle(edge1, edge2, edge3):
    """Make edges follow the cycle ``1->2->3->1``.

    Does this by setting the "previous" and "next" edge
    values on each curved edge.

    Args:
        edge1 (.Curve): First curve in cycle.
        edge2 (.Curve): Second curve in cycle.
        edge3 (.Curve): Third curve in cycle.
    """
    # pylint: disable=protected-access
    edge1._edge_index = 0
    edge1._next_edge = edge2
    edge1._previous_edge = edge3

    edge2._edge_index = 1
    edge2._next_edge = edge3
    edge2._previous_edge = edge1

    edge3._edge_index = 2
    edge3._next_edge = edge1
    edge3._previous_edge = edge2
    # pylint: enable=protected-access


def _ignored_edge_corner(edge_tangent, corner_tangent, corner_previous_edge):
    """Check ignored when a corner lies **inside** another edge.

    Helper for :func:`_ignored_corner` where one of ``s`` and
    ``t`` are ``0``, but **not both**.

    Args:
        edge_tangent (numpy.ndarray): Tangent vector along the edge
            that the intersection occurs in the middle of.
        corner_tangent (numpy.ndarray): Tangent vector at the corner
            where intersection occurs (at the beginning of edge).
        corner_previous_edge (.Curve): Edge that ends at the corner
            intersection (whereas ``corner_tangent`` comes from the edge
            that **begins** at the corner intersection).

    Returns:
        bool: Indicates if the corner intersection should be ignored.
    """
    cross_prod = _helpers.cross_product(edge_tangent, corner_tangent)
    # A negative cross product indicates that ``edge_tangent`` is
    # "inside" / "to the left" of ``corner_tangent`` (due to right-hand rule).
    if cross_prod > 0.0:
        return False

    # Do the same for the **other** tangent at the corner.
    alt_corner_tangent = _curve_helpers.evaluate_hodograph(
        corner_previous_edge._nodes, corner_previous_edge._degree, 1.0)
    # Change the direction of the "in" tangent so that it points "out".
    alt_corner_tangent *= -1.0
    cross_prod = _helpers.cross_product(edge_tangent, alt_corner_tangent)
    return cross_prod <= 0.0


def _ignored_double_corner(intersection, tangent_s, tangent_t):
    """Check if an intersection is an "ignored" double corner.

    Helper for :func:`_ignored_corner` where both ``s`` and
    ``t`` are ``0``.

    Does so by checking if either edge through the ``t`` corner goes
    through the interior of the other surface. An interior check
    is done by checking that a few cross products are positive.

    Args:
        intersection (.Intersection): An intersection to "diagnose".
        tangent_s (numpy.ndarray): The tangent vector to the first curve
            at the intersection.
        tangent_t (numpy.ndarray): The tangent vector to the second curve
            at the intersection.

    Returns:
        bool: Indicates if the corner is to be ignored.
    """
    # Compute the other edge for the ``s`` surface.
    # pylint: disable=protected-access
    prev_edge = intersection.first._previous_edge
    # pylint: enable=protected-access
    alt_tangent_s = _curve_helpers.evaluate_hodograph(
        prev_edge._nodes, prev_edge._degree, 1.0)

    # First check if ``tangent_t`` is interior to the ``s`` surface.
    cross_prod1 = _helpers.cross_product(tangent_s, tangent_t)
    # A positive cross product indicates that ``tangent_t`` is
    # interior to ``tangent_s``. Similar for ``alt_tangent_s``.
    # If ``tangent_t`` is interior to both, then the surfaces
    # do more than just "kiss" at the corner, so the corner should
    # not be ignored.
    if cross_prod1 >= 0.0:
        # Only compute ``cross_prod2`` if we need to.
        cross_prod2 = _helpers.cross_product(alt_tangent_s, tangent_t)
        if cross_prod2 >= 0.0:
            return False

    # If ``tangent_t`` is not interior, we check the other ``t``
    # edge that ends at the corner.
    # pylint: disable=protected-access
    prev_edge = intersection.second._previous_edge
    # pylint: enable=protected-access
    alt_tangent_t = _curve_helpers.evaluate_hodograph(
        prev_edge._nodes, prev_edge._degree, 1.0)
    # Change the direction of the "in" tangent so that it points "out".
    alt_tangent_t *= -1.0

    cross_prod3 = _helpers.cross_product(tangent_s, alt_tangent_t)
    if cross_prod3 >= 0.0:
        # Only compute ``cross_prod4`` if we need to.
        cross_prod4 = _helpers.cross_product(alt_tangent_s, alt_tangent_t)
        if cross_prod4 >= 0.0:
            return False

    # If neither of ``tangent_t`` or ``alt_tangent_t`` are interior
    # to the ``s`` surface, one of two things is true. Either
    # the two surfaces have no interior intersection (1) or the
    # ``s`` surface is bounded by both edges of the ``t`` surface
    # at the corner intersection (2). To detect (2), we only need
    # check if ``tangent_s`` is interior to both ``tangent_t``
    # and ``alt_tangent_t``. ``cross_prod1`` contains
    # (tangent_s) x (tangent_t), so it's negative will tell if
    # ``tangent_s`` is interior. Similarly, ``cross_prod3``
    # contains (tangent_s) x (alt_tangent_t), but we also reversed
    # the sign on ``alt_tangent_t`` so switching the sign back
    # and reversing the arguments in the cross product cancel out.
    return not (cross_prod1 <= 0.0 and cross_prod3 >= 0.0)


def _ignored_corner(intersection, tangent_s, tangent_t):
    """Check if an intersection is an "ignored" corner.

    An "ignored" corner is one where the surfaces just "kiss" at
    the point of intersection but their interiors do not meet.

    We can determine this by comparing the tangent lines from
    the point of intersection.

    .. note::

       This assumes the ``intersection`` has been shifted to the
       beginning of a curve so only checks if ``s == 0.0`` or ``t == 0.0``
       (rather than also checking for ``1.0``).

    .. note::

       This assumes the first and second curves in ``intersection`` are edges
       in a surface, so the code relies on ``previous_edge`` being valid.

    Args:
        intersection (.Intersection): An intersection to "diagnose".
        tangent_s (numpy.ndarray): The tangent vector to the first curve
            at the intersection.
        tangent_t (numpy.ndarray): The tangent vector to the second curve
            at the intersection.

    Returns:
        bool: Indicates if the corner is to be ignored.
    """
    if intersection.s == 0.0:
        # pylint: disable=protected-access
        if intersection.t == 0.0:
            # Double corner.
            return _ignored_double_corner(
                intersection, tangent_s, tangent_t)
        else:
            # s-only corner.
            prev_edge = intersection.first._previous_edge
            return _ignored_edge_corner(tangent_t, tangent_s, prev_edge)
        # pylint: enable=protected-access
    elif intersection.t == 0.0:
        # t-only corner.
        # pylint: disable=protected-access
        prev_edge = intersection.second._previous_edge
        # pylint: enable=protected-access
        return _ignored_edge_corner(tangent_s, tangent_t, prev_edge)
    else:
        # Not a corner.
        return False


def handle_corners(intersection):
    """Mutates an intersection if it is on a corner.

    Does nothing if the intersection happens in the middle of two
    edges.

    If the intersection occurs at the end of the first curve,
    moves it to the beginning of the next edge. Similar for the
    second curve.

    This function is used as a pre-processing step before passing
    an intersection to :func:`classify_intersection`. There, only
    corners that **begin** an edge are considered, since that
    function is trying to determine which edge to **move forward** on.

    .. note::

       This assumes the first and second curves in ``intersection`` are edges
       in a surface, so the code (may) rely on ``next_edge`` and / or
       ``previous_edge`` being valid.

    Args:
        intersection (.Intersection): An intersection to mutate.

    Returns:
        bool: Indicating if the object was changed.
    """
    changed = False
    if intersection.s == 1.0:
        intersection.s = 0.0
        # pylint: disable=protected-access
        intersection.first = intersection.first._next_edge
        # pylint: enable=protected-access
        changed = True
    if intersection.t == 1.0:
        intersection.t = 0.0
        # pylint: disable=protected-access
        intersection.second = intersection.second._next_edge
        # pylint: enable=protected-access
        changed = True

    return changed


def _identifier(intersection):
    """Provide a simple value to identify an intersection.

    Args:
        intersection (.Intersection): The current intersection.

    Returns:
        Tuple[int, float, int, float]: The edge indices for the intersection
        and the parameter values:

        * first edge index
        * first parameter
        * second edge index
        * second parameter
    """
    return (
        intersection.first._edge_index,  # pylint: disable=protected-access
        intersection.s,
        intersection.second._edge_index,  # pylint: disable=protected-access
        intersection.t,
    )


def verify_duplicates(duplicates, uniques):
    """Verify that a set of intersections had expected duplicates.

    .. note::

       This method only considers two intersections as equal if the
       ``s`` and ``t`` values are bit-wise identical.

    Args:
        duplicates (List[.Intersection]): List of intersections
            corresponding to duplicates that were filtered out.
        uniques (List[.Intersection]): List of "final" intersections
            with duplicates filtered out.

    Raises:
        ValueError: If the ``uniques`` are not actually all unique.
        ValueError: If one of the ``duplicates`` does not correspond to
            an intersection in ``uniques``.
        ValueError: If a duplicate occurs only once but does not have
            exactly one of ``s`` and ``t`` equal to ``0.0``.
        ValueError: If a duplicate occurs three times but does not have
            exactly both ``s == t == 0.0``.
        ValueError: If a duplicate occurs a number other than one or three
            times.
    """
    counter = collections.Counter(_identifier(dupe) for dupe in duplicates)
    uniques_keys = set(_identifier(uniq) for uniq in uniques)
    if len(uniques_keys) < len(uniques):
        raise ValueError('Non-unique intersection')

    for key, count in six.iteritems(counter):
        if key not in uniques_keys:
            raise ValueError('Duplicate not among uniques', key)

        if count == 1:
            if (key[1], key[3]).count(0.0) != 1:
                raise ValueError('Count == 1 should be a single corner', key)
        elif count == 3:
            if (key[1], key[3]) != (0.0, 0.0):
                raise ValueError('Count == 3 should be a double corner', key)
        else:
            raise ValueError('Unexpected duplicate count', count)


def _to_front(intersection, intersections, unused):
    """Rotates a node to the "front".

    Helper for :func:`combine_intersections`.

    If a node is at the end of a segment, moves it to the beginning
    of the next segment (at the exact same point).

    .. note::

        This method checks for **exact** endpoints, i.e. parameter
        bitwise identical to ``1.0``. But we should probably allow
        some wiggle room.

    Args:
        intersection (.Intersection): The current intersection.
        intersections (List[.Intersection]): List of all detected
            intersections, provided as a reference for potential
            points to arrive at.
        unused (List[.Intersection]): List of nodes that haven't been
            used yet in an intersection curved polygon

    Returns:
        .Intersection: An intersection to (maybe) move to the beginning
        of the next segment(s).
    """
    changed = False
    if intersection.s == 1.0:
        changed = True
        # pylint: disable=protected-access
        new_intersection = _intersection_helpers.Intersection(
            intersection.first._next_edge, 0.0,
            intersection.second, intersection.t,
            interior_curve=intersection.interior_curve)
        # pylint: enable=protected-access
        intersection = new_intersection

    if intersection.t == 1.0:
        changed = True
        # pylint: disable=protected-access
        new_intersection = _intersection_helpers.Intersection(
            intersection.first, intersection.s,
            intersection.second._next_edge, 0.0,
            interior_curve=intersection.interior_curve)
        # pylint: enable=protected-access
        intersection = new_intersection

    if changed:
        # Make sure we haven't accidentally ignored an existing intersection.
        for other_int in intersections:
            if (other_int.s == intersection.s and
                    other_int.first is intersection.first):
                intersection = other_int
                break

            if (other_int.t == intersection.t and
                    other_int.second is intersection.second):
                intersection = other_int
                break

    if intersection in unused:
        unused.remove(intersection)
    return intersection


def _get_next_first(intersection, intersections):
    """Gets the next node along the current (first) edge.

    Helper for :func:`_get_next` and along with :func:`_get_next_second`, this
    function does the majority of the heavy lifting. **Very** similar to
    :func:`_get_next_second`, but this works with the first curve while the
    other function works with the second.

    Args:
        intersection (.Intersection): The current intersection.
        intersections (List[.Intersection]): List of all detected
            intersections, provided as a reference for potential
            points to arrive at.

    Returns:
        .Intersection: The "next" point along a surface of intersection.
        This will produce the next intersection along the current (first)
        edge or the end of the same edge.
    """
    along_edge = None
    first = intersection.first
    s = intersection.s
    for other_int in intersections:
        other_s = other_int.s
        if other_int.first is first and other_s > s:
            # NOTE: We skip tangent intersections that don't occur
            #       at a corner.
            if (other_s < 1.0 and
                    other_int.interior_curve not in _ACCEPTABLE):
                continue
            if along_edge is None or other_s < along_edge.s:
                along_edge = other_int

    if along_edge is None:
        # If there is no other intersection on the edge, just return
        # the segment end.
        new_intersection = _intersection_helpers.Intersection(
            first, 1.0, None, None,
            interior_curve=IntersectionClassification.first)
        return new_intersection
    else:
        return along_edge


def _get_next_second(intersection, intersections):
    """Gets the next node along the current (second) edge.

    Helper for :func:`_get_next` and along with :func:`_get_next_first`, this
    function does the majority of the heavy lifting. **Very** similar to
    :func:`_get_next_first`, but this works with the second curve while the
    other function works with the first.

    Args:
        intersection (.Intersection): The current intersection.
        intersections (List[.Intersection]): List of all detected
            intersections, provided as a reference for potential
            points to arrive at.

    Returns:
        .Intersection: The "next" point along a surface of intersection.
        This will produce the next intersection along the current (second)
        edge or the end of the same edge.
    """
    along_edge = None
    second = intersection.second
    t = intersection.t
    for other_int in intersections:
        other_t = other_int.t
        if other_int.second is second and other_t > t:
            # NOTE: We skip tangent intersections that don't occur
            #       at a corner.
            if (other_t < 1.0 and
                    other_int.interior_curve not in _ACCEPTABLE):
                continue
            if along_edge is None or other_t < along_edge.t:
                along_edge = other_int

    if along_edge is None:
        # If there is no other intersection on the edge, just return
        # the segment end.
        new_intersection = _intersection_helpers.Intersection(
            None, None, second, 1.0,
            interior_curve=IntersectionClassification.second)
        return new_intersection
    else:
        return along_edge


def _get_next(intersection, intersections, unused):
    """Gets the next node along a given edge.

    Helper for :func:`combine_intersections`, but does the majority
    of the heavy lifting.

    .. note::

        This function returns :class:`.Intersection` objects even
        when the point isn't strictly an intersection. This is
        "incorrect" in some sense, but for now, we don't bother
        implementing a class similar to, but different from,
        :class:`.Intersection` to satisfy this need.

    Args:
        intersection (.Intersection): The current intersection.
        intersections (List[.Intersection]): List of all detected
            intersections, provided as a reference for potential
            points to arrive at.
        unused (List[.Intersection]): List of nodes that haven't been
            used yet in an intersection curved polygon

    Returns:
        .Intersection: The "next" point along a surface of intersection.
        This will produce the next intersection along the current edge or
        the end of the current edge.

    Raises:
        ValueError: If the intersection is not classified as
            :attr:`~.IntersectionClassification.first` or
            :attr:`~.IntersectionClassification.second`.
    """
    result = None
    if intersection.interior_curve is IntersectionClassification.first:
        result = _get_next_first(intersection, intersections)
    elif intersection.interior_curve is IntersectionClassification.second:
        result = _get_next_second(intersection, intersections)
    else:
        raise ValueError('Cannot get next node if not starting from '
                         '"first" or "second".')

    if result in unused:
        unused.remove(result)
    return result


def _ends_to_curve(start_node, end_node):
    """Convert a "pair" of intersection nodes to a curve segment.

    .. note::

       This function could specialize to the first or second segment
       attached to ``start_node`` and ``end_node``. We determine
       first / second based on the classification of ``start_node``,
       but the callers of this function could provide that information /
       isolate the base curve and the two parameters for us.

    .. note::

       This only checks the classification of the ``start_node``.

    Args:
        start_node (.Intersection): The beginning of a segment.
        end_node (.Intersection): The end of (the same) segment.

    Returns:
        .Curve: The segment between the nodes.

    Raises:
        ValueError: If the ``start_node`` and ``end_node`` disagree on
            the first curve when classified as "first" or disagree on
            the second curve when classified as "second".
        ValueError: If the ``start_node`` is not classified as
            :attr:`~.IntersectionClassification.first` or
            :attr:`~.IntersectionClassification.second`.
    """
    if start_node.interior_curve is IntersectionClassification.first:
        first = start_node.first
        if end_node.first is not first:
            raise ValueError(_WRONG_CURVE)
        return first.specialize(start_node.s, end_node.s)
    elif start_node.interior_curve is IntersectionClassification.second:
        second = start_node.second
        if end_node.second is not second:
            raise ValueError(_WRONG_CURVE)
        return second.specialize(start_node.t, end_node.t)
    else:
        raise ValueError('Segment start must be classified as '
                         '"first" or "second".')


def _to_curved_polygon(surface):
    """Convert a surface to a curved polygon.

    This is a helper for :func:`combine_intersections` to allow
    returning a curved polygon when the intersection is the
    entirety of one surface.

    Args:
        surface (.Surface): The surface to convert.

    Returns:
        .CurvedPolygon: The converted object.
    """
    edges = surface._get_edges()  # pylint: disable=protected-access
    return curved_polygon.CurvedPolygon(*edges, _verify=False)


def _no_intersections(surface1, surface2):
    r"""Determine if one surface is in the other.

    Helper for :func:`combine_intersections` that handles the case
    of no points of intersection. In this case, either the surfaces
    are disjoint or one is fully contained in the other.

    To check containment, it's enough to check if one of the corners
    is contained in the other surface.

    Args:
        surface1 (.Surface): First surface in intersection (assumed in
            :math:\mathbf{R}^2`).
        surface2 (.Surface): Second surface in intersection (assumed in
            :math:\mathbf{R}^2`).

    Returns:
        list: Either an empty list if one surface isn't contained
        in the other. Otherwise, the list will have a single
        :class:`.CurvedPolygon` corresponding to the internal surface.
    """
    # NOTE: We want the nodes to be 1x2 but accessing ``nodes1[[0], :]``
    #       and ``nodes2[[0], :]`` makes a copy while the accesses
    #       below **do not** copy. See
    #       (https://docs.scipy.org/doc/numpy-1.6.0/reference/
    #        arrays.indexing.html#advanced-indexing)
    corner1 = surface1._nodes[0, :].reshape((1, 2))
    if surface2.locate(corner1) is not None:
        return [_to_curved_polygon(surface1)]

    corner2 = surface2._nodes[0, :].reshape((1, 2))
    if surface1.locate(corner2) is not None:
        return [_to_curved_polygon(surface2)]

    return []


def _tangent_only_intersections(intersections, surface1, surface2):
    """Determine intersection in the case of only-tangent intersections.

    If the only intersections are tangencies, then either the surfaces
    are tangent but don't meet ("kissing" edges) or one surface is
    internally tangent to the other.

    Thus we expect every intersection in ``intersections`` to be
    classified as :attr:`~.IntersectionClassification.tangent_first`,
    :attr:`~.IntersectionClassification.tangent_second` or
    :attr:`~.IntersectionClassification.opposed`.

    What's more, we expect all intersections to be classified the same for
    a given pairing.

    Args:
        intersections (list): A list of :class:`.Intersection` objects
            produced by :func:`.all_intersections` applied to each of
            the 9 edge-edge pairs from a surface-surface pairing.
        surface1 (.Surface): First surface in intersection.
        surface2 (.Surface): Second surface in intersection.

    Returns:
        list: Either an empty list if one surface isn't contained
        in the other. Otherwise, the list will have a single
        :class:`.CurvedPolygon` corresponding to the internal surface.

    Raises:
        ValueError: If there are intersections of more than one type among
            :attr:`~.IntersectionClassification.tangent_first`,
            :attr:`~.IntersectionClassification.tangent_second` or
            :attr:`~.IntersectionClassification.opposed`.
        ValueError: If there is a unique classification, but it isn't one
            of the tangent types.
    """
    all_types = set([intersection.interior_curve
                     for intersection in intersections])
    if len(all_types) != 1:
        raise ValueError('Unexpected value, types should all match',
                         all_types)
    point_type = all_types.pop()
    ignored_types = (
        IntersectionClassification.opposed,
        IntersectionClassification.ignored_corner,
    )
    if point_type in ignored_types:
        return []
    elif point_type is IntersectionClassification.tangent_first:
        return [_to_curved_polygon(surface1)]
    elif point_type is IntersectionClassification.tangent_second:
        return [_to_curved_polygon(surface2)]
    else:
        raise ValueError('Point type not for tangency', point_type)


def _basic_interior_combine(intersections):
    """Combine intersections that don't involve tangencies.

    Helper for :func:`combine_intersections`.

    .. note::

       This helper assumes ``intersections`` isn't empty, but doesn't
       enforce it.

    Args:
        intersections (list): A list of :class:`.Intersection` objects
            produced by :func:`.all_intersections` applied to each of
            the 9 edge-edge pairs from a surface-surface pairing.

    Returns:
        List[.CurvedPolygon]: All of the intersections encountered.
    """
    unused = [intersection for intersection in intersections
              if intersection.interior_curve in _ACCEPTABLE]
    result = []
    while unused:
        start = unused.pop()
        curr_node = start
        next_node = _get_next(start, intersections, unused)
        edge_ends = [(curr_node, next_node)]
        while next_node is not start:
            curr_node = _to_front(next_node, intersections, unused)
            # NOTE: We also check to break when moving a corner node
            #       to the front. This is because ``intersections``
            #       de-duplicates corners by selecting the one
            #       (of 2 or 4 choices) at the front of segment(s).
            if curr_node is start:
                break
            next_node = _get_next(curr_node, intersections, unused)
            edge_ends.append((curr_node, next_node))

        edge_gen = (_ends_to_curve(*pair) for pair in edge_ends)
        result.append(curved_polygon.CurvedPolygon(*edge_gen, _verify=False))

    return result


def combine_intersections(intersections, surface1, surface2):
    """Combine curve-curve intersections into curved polygon(s).

    Does so assuming each intersection lies on an edge of one of
    two :class:`.Surface`-s.

    .. note ::

       This assumes that each ``intersection`` has been classified via
       :func:`classify_intersection`.

    Args:
        intersections (list): A list of :class:`.Intersection` objects
            produced by :func:`.all_intersections` applied to each of
            the 9 edge-edge pairs from a surface-surface pairing.
        surface1 (.Surface): First surface in intersection.
        surface2 (.Surface): Second surface in intersection.

    Returns:
        List[~bezier.curved_polygon.CurvedPolygon]: A list of curved polygons
        that compose the intersected objects.
    """
    if len(intersections) == 0:
        return _no_intersections(surface1, surface2)

    result = _basic_interior_combine(intersections)
    if len(result) > 0:
        return result

    return _tangent_only_intersections(intersections, surface1, surface2)


def evaluate_barycentric(nodes, degree, lambda1, lambda2, lambda3):
    r"""Compute a point on a surface.

    Evaluates :math:`B\left(\lambda_1, \lambda_2, \lambda_3\right)` for a
    B |eacute| zier surface / triangle defined by ``nodes``.

    Args:
        nodes (numpy.ndarray): Control point nodes that define the surface.
        degree (int): The degree of the surface define by ``nodes``.
        lambda1 (float): Parameter along the reference triangle.
        lambda2 (float): Parameter along the reference triangle.
        lambda3 (float): Parameter along the reference triangle.

    Returns:
        numpy.ndarray: The evaluate point as a ``1xD`` array (where ``D``
        is the ambient dimension where ``nodes`` reside).
    """
    if degree == 1:
        weights = np.array([
            [lambda1, lambda2, lambda3],
        ])
    elif degree == 2:
        weights = np.array([
            [
                lambda1 * lambda1,
                2.0 * lambda1 * lambda2,
                lambda2 * lambda2,
                2.0 * lambda1 * lambda3,
                2.0 * lambda2 * lambda3,
                lambda3 * lambda3,
            ]
        ])
    elif degree == 3:
        weights = np.array([
            [
                lambda1 * lambda1 * lambda1,
                3.0 * lambda1 * lambda1 * lambda2,
                3.0 * lambda1 * lambda2 * lambda2,
                lambda2 * lambda2 * lambda2,
                3.0 * lambda1 * lambda1 * lambda3,
                6.0 * lambda1 * lambda2 * lambda3,
                3.0 * lambda2 * lambda2 * lambda3,
                3.0 * lambda1 * lambda3 * lambda3,
                3.0 * lambda2 * lambda3 * lambda3,
                lambda3 * lambda3 * lambda3,
            ]
        ])
    else:
        result = nodes
        for reduced_deg in six.moves.xrange(degree, 0, -1):
            result = de_casteljau_one_round(
                result, reduced_deg, lambda1, lambda2, lambda3)
        return result

    return weights.dot(nodes)  # pylint: disable=no-member


class IntersectionClassification(enum.Enum):
    """Enum classifying the "interior" curve in an intersection.

    Provided as the output values for :func:`.classify_intersection`.
    """

    first = 'first'
    """The first curve is on the interior."""

    second = 'second'
    """The second curve is on the interior."""

    opposed = 'opposed'
    """Tangent intersection with opposed interiors."""

    tangent_first = 'tangent_first'
    """Tangent intersection, first curve is on the interior."""

    tangent_second = 'tangent_second'
    """Tangent intersection, second curve is on the interior."""

    ignored_corner = 'ignored_corner'
    """Intersection at a corner, interiors don't intersect."""


# NOTE: This constant must be define **after** ``IntersectionClassification``
#       is, hence we can't define it at the top with the other constants.
_ACCEPTABLE = (
    IntersectionClassification.first,
    IntersectionClassification.second,
)

# pylint: disable=invalid-name
if _speedup is None:  # pragma: NO COVER
    de_casteljau_one_round = _de_casteljau_one_round
else:
    de_casteljau_one_round = _speedup.speedup.de_casteljau_one_round
# pylint: enable=invalid-name
