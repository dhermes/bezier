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


_MAX_POLY_SUBDIVISIONS = 5
_MAX_LOCATE_SUBDIVISIONS = 20
_LOCATE_EPS = 2.0**(-47)
_SIGN = np.sign  # pylint: disable=no-member
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
], dtype=float) / 2.0
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
], dtype=float) / 4.0
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
], dtype=float) / 8.0
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
], dtype=float) / 16.0
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
], dtype=float)
_QUADRATIC_TO_BERNSTEIN = np.array([
    [ 2, 0,  0, 0, 0,  0],  # noqa: E201
    [-1, 4, -1, 0, 0,  0],
    [ 0, 0,  2, 0, 0,  0],  # noqa: E201
    [-1, 0,  0, 4, 0, -1],
    [ 0, 0, -1, 0, 4, -1],  # noqa: E201
    [ 0, 0,  0, 0, 0,  2],  # noqa: E201
], dtype=float) / 2.0
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
], dtype=float) / 16.0
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
], dtype=float)
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
    sub_polys = [poly_surface]
    signs = set()
    for _ in six.moves.xrange(_MAX_POLY_SUBDIVISIONS):
        undecided = []
        for poly in sub_polys:
            # Avoid an unnecessarily copying the nodes.
            # pylint: disable=protected-access
            nodes = poly._nodes
            # pylint: enable=protected-access
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
       >>> actual_det == np_det
       False
       >>> abs(actual_det - np_det) < 1e-16
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
    # pylint: disable=no-member
    bernstein = _QUADRATIC_TO_BERNSTEIN.dot(jac_at_nodes)
    # pylint: enable=no-member
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
        numpy.ndarray: Coefficients in Bernstein basis.
    """
    # First evaluate the Jacobian at each of the 15 nodes
    # in the quartic triangle.
    # pylint: disable=no-member
    jac_parts = _CUBIC_JACOBIAN_HELPER.dot(nodes)
    # pylint: enable=no-member
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


def de_casteljau_one_round(nodes, degree, lambda1, lambda2, lambda3):
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
        sum_x += candidate.base_x
        sum_y += candidate.base_y
        sum_width += candidate.width

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

       >>> surface = bezier.Surface(np.array([
       ...     [0.0, 0.0],
       ...     [1.0, 0.0],
       ...     [2.0, 0.0],
       ...     [2.0, 1.0],
       ...     [2.0, 2.0],
       ...     [0.0, 2.0],
       ... ]))
       >>> surface.is_valid
       True
       >>> x_val, y_val = surface.evaluate_cartesian(0.25, 0.5)
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
    surf_x, surf_y = surface.evaluate_cartesian(s, t)
    if surf_x == x_val and surf_y == y_val:
        # No refinement is needed.
        return s, t

    nodes = surface._nodes  # pylint: disable=protected-access
    # Compute Jacobian nodes / stack them horizatonally.
    jac_both = np.hstack([
        jacobian_s(nodes, surface.degree, surface.dimension),
        jacobian_t(nodes, surface.degree, surface.dimension),
    ])

    lambda1 = 1.0 - s - t
    # The degree of the jacobian is one less.
    for reduced_deg in six.moves.xrange(surface.degree - 1, 0, -1):
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
            nodes = candidate._nodes  # pylint: disable=protected-access
            if _helpers.contains(nodes, x_val, y_val):
                next_candidates.extend(candidate.subdivide())

        candidates = next_candidates

    if not candidates:
        return None

    # We take the average of all centroids from the candidates
    # that may contain the point.
    s_approx, t_approx = _mean_centroid(candidates)
    s, t = newton_refine(surface, x_val, y_val, s_approx, t_approx)

    actual = surface.evaluate_cartesian(s, t)
    expected = np.array([x_val, y_val])
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
                   classify-intersection7

       import numpy as np
       import bezier
       from bezier import _curve_helpers
       from bezier._intersection_helpers import Intersection
       from bezier._surface_helpers import classify_intersection

       def hodograph(curve, s):
           return _curve_helpers.evaluate_hodograph(
               curve._nodes, curve.degree, s)

       def curvature(curve, s):
           nodes = curve._nodes
           tangent = _curve_helpers.evaluate_hodograph(
               nodes, curve.degree, s)
           return _curve_helpers.get_curvature(
               nodes, curve.degree, tangent, s)

    .. doctest:: classify-intersection1
       :options: +NORMALIZE_WHITESPACE

       >>> curve1 = bezier.Curve(np.array([
       ...     [1.0 , 0.0 ],
       ...     [1.75, 0.25],
       ...     [2.0 , 1.0 ],
       ... ]))
       >>> curve2 = bezier.Curve(np.array([
       ...     [0.0   , 0.0   ],
       ...     [1.6875, 0.0625],
       ...     [2.0   , 0.5   ],
       ... ]))
       >>> s, t = 0.25, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([ True, True], dtype=bool)
       >>> tangent1 = hodograph(curve1, s)
       >>> tangent1
       array([ 1.25, 0.75])
       >>> tangent2 = hodograph(curve2, t)
       >>> tangent2
       array([ 2. , 0.5])
       >>> intersection = Intersection(curve1, s, curve2, t)
       >>> classify_intersection(intersection)
       <IntersectionClassification.first: 'first'>

    .. testcleanup:: classify-intersection1

       import make_images
       make_images.classify_intersection1(
           s, curve1, tangent1, curve2, tangent2)

    We determine the interior / left one by using the `right-hand rule`_:
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

    If the cross-product quantity
    :math:`B_1'(s) \times B_2'(t) = x_1'(s) y_2'(t) - x_2'(t) y_1'(s)`
    is positive, then the first curve is "to the right", i.e. the second
    curve is interior. If the cross-product is negative, the first
    curve is interior.

    When :math:`B_1'(s) \times B_2'(t) = 0`, the tangent
    vectors are parallel, i.e. the intersection is a point of tangency:

    .. image:: images/classify_intersection2.png
       :align: center

    .. doctest:: classify-intersection2
       :options: +NORMALIZE_WHITESPACE

       >>> curve1 = bezier.Curve(np.array([
       ...     [1.0, 0.0],
       ...     [1.5, 1.0],
       ...     [2.0, 0.0],
       ... ]))
       >>> curve2 = bezier.Curve(np.array([
       ...     [0.0, 0.0],
       ...     [1.5, 1.0],
       ...     [3.0, 0.0],
       ... ]))
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([ True, True], dtype=bool)
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

       >>> curve1 = bezier.Curve(np.array([
       ...     [2.0, 0.0],
       ...     [1.5, 1.0],
       ...     [1.0, 0.0],
       ... ]))
       >>> curve2 = bezier.Curve(np.array([
       ...     [3.0, 0.0],
       ...     [1.5, 1.0],
       ...     [0.0, 0.0],
       ... ]))
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([ True, True], dtype=bool)
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

       >>> curve1 = bezier.Curve(np.array([
       ...     [2.0, 0.0],
       ...     [1.5, 1.0],
       ...     [1.0, 0.0],
       ... ]))
       >>> curve2 = bezier.Curve(np.array([
       ...     [0.0, 0.0],
       ...     [1.5, 1.0],
       ...     [3.0, 0.0],
       ... ]))
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([ True, True], dtype=bool)
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

       >>> curve1 = bezier.Curve(np.array([
       ...     [1.0, 0.0],
       ...     [1.5, 1.0],
       ...     [2.0, 0.0],
       ... ]))
       >>> curve2 = bezier.Curve(np.array([
       ...     [3.0, 0.0],
       ...     [1.5, 1.0],
       ...     [0.0, 0.0],
       ... ]))
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([ True, True], dtype=bool)
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

       >>> curve1 = bezier.Curve(np.array([
       ...     [ 0.375,  0.0625],
       ...     [-0.125, -0.0625],
       ...     [-0.125,  0.0625],
       ... ]))
       >>> curve2 = bezier.Curve(np.array([
       ...     [ 0.75,  0.25],
       ...     [-0.25, -0.25],
       ...     [-0.25,  0.25],
       ... ]))
       >>> s, t = 0.5, 0.5
       >>> curve1.evaluate(s) == curve2.evaluate(t)
       array([ True, True], dtype=bool)
       >>> hodograph(curve1, s)
       array([-0.5, 0. ])
       >>> hodograph(curve2, t)
       array([-1., 0.])
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

       >>> curve1a = bezier.Curve(np.array([
       ...     [0.0, 0.0 ],
       ...     [4.5, 0.0 ],
       ...     [9.0, 2.25],
       ... ]))
       >>> curve2 = bezier.Curve(np.array([
       ...     [11.25, 0.0],
       ...     [ 9.0 , 4.5],
       ...     [ 2.75, 1.0],
       ... ]))
       >>> s, t = 1.0, 0.375
       >>> curve1a.evaluate(s) == curve2.evaluate(t)
       array([ True, True], dtype=bool)
       >>> intersection = Intersection(curve1a, s, curve2, t)
       >>> classify_intersection(intersection)
       Traceback (most recent call last):
         ...
       ValueError: ('Intersection occurs at the end of an edge',
                    's', 1.0, 't', 0.375)
       >>>
       >>> curve1b = bezier.Curve(np.array([
       ...     [9.0, 2.25 ],
       ...     [4.5, 2.375],
       ...     [0.0, 2.5  ],
       ... ]))
       >>> curve1b.evaluate(0.0) == curve2.evaluate(t)
       array([ True, True], dtype=bool)
       >>> intersection = Intersection(curve1b, 0.0, curve2, t)
       >>> classify_intersection(intersection)
       <IntersectionClassification.first: 'first'>

    .. testcleanup:: classify-intersection7

       import make_images
       make_images.classify_intersection7(s, curve1a, curve1b, curve2)

    .. note::

       This assumes the intersection occurs in :math:`\mathbf{R}^2`
       but doesn't check this.

    .. note::

       This function doesn't allow wiggle room / round-off when checking
       endpoints, nor when checking if the cross-product is near zero,
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
        intersection.left._nodes,  # pylint: disable=protected-access
        intersection.left.degree, intersection.s)
    tangent2 = _curve_helpers.evaluate_hodograph(
        intersection.right._nodes,  # pylint: disable=protected-access
        intersection.right.degree, intersection.t)

    # Take the cross-product of tangent vectors to determine which one
    # is more "to the left".
    cross_prod = _helpers.cross_product(
        np.array([tangent1]), np.array([tangent2]))
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
        tangent1 (numpy.ndarray): The tangent vector on the ``left`` curve
            at the intersection.
        tangent2 (numpy.ndarray): The tangent vector on the ``right`` curve
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
    dot_prod = tangent1.dot(tangent2)
    # NOTE: When computing curvatures we assume that we don't have lines
    #       here, because lines that are tangent at an intersection are
    #       parallel and we don't handle that case.
    curvature1 = _curve_helpers.get_curvature(
        intersection.left._nodes,  # pylint: disable=protected-access
        intersection.left.degree, tangent1, intersection.s)
    curvature2 = _curve_helpers.get_curvature(
        intersection.right._nodes,  # pylint: disable=protected-access
        intersection.right.degree, tangent2, intersection.t)
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


def handle_corners(intersection):
    """Mutates an intersection if it is on a corner.

    Does nothing if the intersection happens in the middle of two
    edges.

    If the intersection occurs at the end of the ``left`` curve,
    moves it to the beginning of the next edge. Similar for the
    ``right`` curve.

    This function is used as a pre-processing step before passing
    an intersection to :func:`classify_intersection`. There, only
    corners that **begin** an edge are considered, since that
    function is trying to determine which edge to **move forward** on.

    .. note::

       This assumes the ``left`` and ``right`` curves are edges
       in a surface, so the code (may) rely on ``next_edge``
       and / or ``previous_edge`` being valid.

    Args:
        intersection (.Intersection): An intersection to mutate.

    Returns:
        bool: Indicating if the object was changed.
    """
    changed = False
    if intersection.s == 1.0:
        # pylint: disable=protected-access
        intersection._s_val = 0.0
        intersection._left = intersection.left.next_edge
        # pylint: enable=protected-access
        changed = True
    if intersection.t == 1.0:
        # pylint: disable=protected-access
        intersection._t_val = 0.0
        intersection._right = intersection.right.next_edge
        # pylint: enable=protected-access
        changed = True

    return changed


def _identifier(intersection):
    """Provide a simple value to identify an intersection.

    Args:
        intersection (.Intersection): The current intersection.

    Returns:
        Tuple[int, float, int, float]: The edge indices (left first,
        right third) for the intersection and the parameter values
        (left second, right fourth).
    """
    return (intersection.left.edge_index, intersection.s,
            intersection.right.edge_index, intersection.t)


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
        new_intersection = _intersection_helpers.Intersection(
            intersection.left.next_edge, 0.0,
            intersection.right, intersection.t)
        new_intersection.interior_curve = intersection.interior_curve
        intersection = new_intersection

    if intersection.t == 1.0:
        changed = True
        new_intersection = _intersection_helpers.Intersection(
            intersection.left, intersection.s,
            intersection.right.next_edge, 0.0)
        new_intersection.interior_curve = intersection.interior_curve
        intersection = new_intersection

    if changed:
        # Make sure we haven't accidentally ignored an existing intersection.
        for other_int in intersections:
            if (other_int.s == intersection.s and
                    other_int.left is intersection.left):
                intersection = other_int
                break

            if (other_int.t == intersection.t and
                    other_int.right is intersection.right):
                intersection = other_int
                break

    if intersection in unused:
        unused.remove(intersection)
    return intersection


# pylint: disable=too-many-branches
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
    acceptable = (IntersectionClassification.first,
                  IntersectionClassification.second)

    result = None
    if intersection.interior_curve is IntersectionClassification.first:
        along_edge = None
        left = intersection.left
        s = intersection.s
        for other_int in intersections:
            other_s = other_int.s
            if other_int.left is left and other_s > s:
                # NOTE: We skip tangent intersections that don't occur
                #       at a corner.
                if (other_s < 1.0 and
                        other_int.interior_curve not in acceptable):
                    continue
                if along_edge is None or other_s < along_edge.s:
                    along_edge = other_int

        if along_edge is None:
            # Just return the segment end.
            new_intersection = _intersection_helpers.Intersection(
                intersection.left, 1.0, None, None)
            new_intersection.interior_curve = IntersectionClassification.first
            result = new_intersection
        else:
            result = along_edge
    elif intersection.interior_curve is IntersectionClassification.second:
        along_edge = None
        right = intersection.right
        t = intersection.t
        for other_int in intersections:
            other_t = other_int.t
            if other_int.right is right and other_t > t:
                # NOTE: We skip tangent intersections that don't occur
                #       at a corner.
                if (other_t < 1.0 and
                        other_int.interior_curve not in acceptable):
                    continue
                if along_edge is None or other_t < along_edge.t:
                    along_edge = other_int

        if along_edge is None:
            # Just return the segment end.
            new_intersection = _intersection_helpers.Intersection(
                None, None, intersection.right, 1.0)
            # pylint: disable=redefined-variable-type
            new_intersection.interior_curve = IntersectionClassification.second
            # pylint: enable=redefined-variable-type
            result = new_intersection
        else:
            result = along_edge
    else:
        raise ValueError('Cannot get next node if not starting from '
                         '"first" or "second".')

    if result in unused:
        unused.remove(result)
    return result
# pylint: enable=too-many-branches


def _ends_to_curve(start_node, end_node):
    """Convert a "pair" of intersection nodes to a curve segment.

    .. note::

       This function determines "left" or "right" curve based on the
       classification of ``start_node``, but the callers of this
       function could provide that information / isolate the
       base curve and the two parameters for us.

    .. note::

       This only checks the classification of the ``start_node``.

    Args:
        start_node (.Intersection): The beginning of a segment.
        end_node (.Intersection): The end of (the same) segment.

    Returns:
        .Curve: The segment between the nodes.

    Raises:
        ValueError: If the ``start_node`` and ``end_node`` disagree on
            the "left" curve when classified as "first" or disagree on
            the "right" curve when classified as "second".
        ValueError: If the ``start_node`` is not classified as
            :attr:`~.IntersectionClassification.first` or
            :attr:`~.IntersectionClassification.second`.
    """
    if start_node.interior_curve is IntersectionClassification.first:
        left = start_node.left
        if end_node.left is not left:
            raise ValueError(_WRONG_CURVE)
        return left.specialize(start_node.s, end_node.s)
    elif start_node.interior_curve is IntersectionClassification.second:
        right = start_node.right
        if end_node.right is not right:
            raise ValueError(_WRONG_CURVE)
        return right.specialize(start_node.t, end_node.t)
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
    return curved_polygon.CurvedPolygon(*edges)


def _no_intersections(surface1, surface2):
    """Determine if one surface is in the other.

    Helper for :func:`combine_intersections` that handles the case
    of no points of intersection. In this case, either the surfaces
    are disjoint or one is fully contained in the other.

    To check containment, it's enough to check if one of the corners
    is contained in the other surface.

    Args:
        surface1 (.Surface): First surface in intersection.
        surface2 (.Surface): Second surface in intersection.

    Returns:
        list: Either an empty list if one surface isn't contained
        in the other. Otherwise, the list will have a single
        :class:`.CurvedPolygon` corresponding to the internal surface.
    """
    nodes1 = surface1._nodes  # pylint: disable=protected-access
    corner1 = nodes1[[0], :]
    if surface2.locate(corner1) is not None:
        return [_to_curved_polygon(surface1)]

    nodes2 = surface2._nodes  # pylint: disable=protected-access
    corner2 = nodes2[[0], :]
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
    if point_type is IntersectionClassification.opposed:
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
    acceptable = (IntersectionClassification.first,
                  IntersectionClassification.second)
    unused = [intersection for intersection in intersections
              if intersection.interior_curve in acceptable]
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

        result.append(curved_polygon.CurvedPolygon(*(
            _ends_to_curve(*pair)
            for pair in edge_ends
        )))

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


class IntersectionClassification(enum.Enum):
    """Enum classifying the "interior" curve in an intersection."""
    first = 'first'
    second = 'second'
    opposed = 'opposed'
    tangent_first = 'tangent_first'
    tangent_second = 'tangent_second'
