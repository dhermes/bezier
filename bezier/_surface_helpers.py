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

"""Private helper methods for :mod:`bezier.surface`."""


import functools
import operator

import numpy as np
import six


MAX_SUBDIVISIONS = 5
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
LINEAR_JACOBIAN_HELPER = np.array([
    [-1, 1, 0],
    [-1, 0, 1],
], dtype=float)
# The Jacobian of a quadratric (in any dimension) as given by
# dB/ds = [-2L1, 2(L1 - L2), 2L2, -2L3, 2L3, 0] * nodes
# dB/dt = [-2L1, -2L2, 0, 2(L1 - L3), 2L2, 2L3] * nodes
# We evaluate this at each of the 6 points in the quadratic
# triangle and then stack them (2 rows * 6 = 12 rows)
# pylint: disable=bad-whitespace
QUADRATIC_JACOBIAN_HELPER = np.array([
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
QUADRATIC_TO_BERNSTEIN = np.array([
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
CUBIC_JACOBIAN_HELPER = np.array([
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
QUARTIC_TO_BERNSTEIN = np.array([
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
], dtype=float) / 36.0
# pylint: enable=bad-whitespace


def polynomial_sign(poly_surface):
    r"""Determine the "sign" of a polynomial on the reference triangle.

    Checks if a polynomial :math:`p(s, t)` is positive, negative
    or mixed sign on the reference triangle.

    Does this by utilizing the Bezier form of :math:`p`: it is a
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
    for _ in six.moves.xrange(MAX_SUBDIVISIONS):
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
            MAX_SUBDIVISIONS)


def quadratic_jacobian_polynomial(nodes):
    r"""Compute the Jacobian determinant of a quadratic surface.

    Converts :math:`\det(J(s, t))` to a polynomial on the reference
    triangle and represents it as a surface object.

    .. note::

        This assumes that ``nodes`` is 6x2 but doesn't verify this.
        (However, the multiplication by ``QUADRATIC_JACOBIAN_HELPER``
        would fail if ``nodes`` wasn't 6xN and then the ensuing
        determinants would fail if there weren't 2 columns.)

    Args:
        nodes (numpy.ndarray): A 6x2 array of nodes in a surface.

    Returns:
        numpy.ndarray: Coefficients in Bernstein basis.
    """
    # First evaluate the Jacobian at each of the 6 nodes.
    # pylint: disable=no-member
    jac_parts = QUADRATIC_JACOBIAN_HELPER.dot(nodes)
    # pylint: enable=no-member
    jac_at_nodes = np.empty((6, 1))
    jac_at_nodes[0, 0] = np.linalg.det(jac_parts[:2, :])
    jac_at_nodes[1, 0] = np.linalg.det(jac_parts[2:4, :])
    jac_at_nodes[2, 0] = np.linalg.det(jac_parts[4:6, :])
    jac_at_nodes[3, 0] = np.linalg.det(jac_parts[6:8, :])
    jac_at_nodes[4, 0] = np.linalg.det(jac_parts[8:10, :])
    jac_at_nodes[5, 0] = np.linalg.det(jac_parts[10:, :])

    # Convert the nodal values to the Bernstein basis...
    # pylint: disable=no-member
    bernstein = QUADRATIC_TO_BERNSTEIN.dot(jac_at_nodes)
    # pylint: enable=no-member
    return bernstein


def cubic_jacobian_polynomial(nodes):
    r"""Compute the Jacobian determinant of a cubic surface.

    Converts :math:`\det(J(s, t))` to a polynomial on the reference
    triangle and represents it as a surface object.

    .. note::

        This assumes that ``nodes`` is 10x2 but doesn't verify this.
        (However, the multiplication by ``CUBIC_JACOBIAN_HELPER``
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
    jac_parts = CUBIC_JACOBIAN_HELPER.dot(nodes)
    # pylint: enable=no-member
    jac_at_nodes = np.empty((15, 1))
    jac_at_nodes[0, 0] = np.linalg.det(jac_parts[:2, :])
    jac_at_nodes[1, 0] = np.linalg.det(jac_parts[2:4, :])
    jac_at_nodes[2, 0] = np.linalg.det(jac_parts[4:6, :])
    jac_at_nodes[3, 0] = np.linalg.det(jac_parts[6:8, :])
    jac_at_nodes[4, 0] = np.linalg.det(jac_parts[8:10, :])
    jac_at_nodes[5, 0] = np.linalg.det(jac_parts[10:12, :])
    jac_at_nodes[6, 0] = np.linalg.det(jac_parts[12:14, :])
    jac_at_nodes[7, 0] = np.linalg.det(jac_parts[14:16, :])
    jac_at_nodes[8, 0] = np.linalg.det(jac_parts[16:18, :])
    jac_at_nodes[9, 0] = np.linalg.det(jac_parts[18:20, :])
    jac_at_nodes[10, 0] = np.linalg.det(jac_parts[20:22, :])
    jac_at_nodes[11, 0] = np.linalg.det(jac_parts[22:24, :])
    jac_at_nodes[12, 0] = np.linalg.det(jac_parts[24:26, :])
    jac_at_nodes[13, 0] = np.linalg.det(jac_parts[26:28, :])
    jac_at_nodes[14, 0] = np.linalg.det(jac_parts[28:, :])

    # Convert the nodal values to the Bernstein basis...
    # pylint: disable=no-member
    bernstein = QUARTIC_TO_BERNSTEIN.dot(jac_at_nodes)
    # pylint: enable=no-member
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
