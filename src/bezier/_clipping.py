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

r"""Proof-of-concept for B |eacute| zier clipping.

.. _algorithm: https://dx.doi.org/10.1016/0010-4485(90)90039-F

The B |eacute| zier clipping `algorithm`_ is used to intersect
two planar B |eacute| zier curves. It proceeds by using "fat lines"
to recursively prune the region of accepted parameter ranges until
the ranges converge to points. (A "fat line" is a rectangular region of a
bounded distance from the line connecting the start and end points of a
B |eacute| zier curve.)

It has quadratic convergence. It can be used to find tangent intersections,
which is the primary usage within ``bezier``.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


import numpy as np
import six


def compute_implicit_line(nodes):
    """Compute the implicit form of the line connecting curve endpoints.

    .. note::

       This assumes, but does not check, that the first and last nodes
       in ``nodes`` are different.

    Computes :math:`a, b` and :math:`c` in the normalized implicit equation
    for the line

    .. math::

       ax + by + c = 0

    where :math:`a^2 + b^2 = 1` (only unique up to sign).

    Args:
        nodes (numpy.ndarray): ``2 x N`` array of nodes in a curve.
            The line will be (directed) from the first to last
            node in ``nodes``.

    Returns:
        Tuple[float, float, float]: The triple of

        * The :math:`x` coefficient :math:`a`
        * The :math:`y` coefficient :math:`b`
        * The constant :math:`c`
    """
    delta = nodes[:, -1] - nodes[:, 0]
    length = np.linalg.norm(delta, ord=2)
    # Normalize and rotate 90 degrees to the "left".
    coeff_a = -delta[1] / length
    coeff_b = delta[0] / length
    # c = - ax - by = (delta[1] x - delta[0] y) / L
    # NOTE: We divide by ``length`` at the end to "put off" rounding.
    coeff_c = (delta[1] * nodes[0, 0] - delta[0] * nodes[1, 0]) / length

    return coeff_a, coeff_b, coeff_c


def compute_fat_line(nodes):
    """Compute the "fat line" around a B |eacute| zier curve.

    Both computes the implicit (normalized) form

    .. math::

       ax + by + c = 0

    for the line connecting the first and last node in ``nodes``.
    Also computes the maximum and minimum distances to that line
    from each control point.

    Args:
        nodes (numpy.ndarray): ``2 x N`` array of nodes in a curve.

    Returns:
        Tuple[float, float, float, float, float]: The 5-tuple of

        * The :math:`x` coefficient :math:`a`
        * The :math:`y` coefficient :math:`b`
        * The constant :math:`c`
        * The "minimum" distance to the fat line among the control points.
        * The "maximum" distance to the fat line among the control points.
    """
    coeff_a, coeff_b, coeff_c = compute_implicit_line(nodes)

    # NOTE: This assumes, but does not check, that there are two rows.
    _, num_nodes = nodes.shape
    d_min = 0.0
    d_max = 0.0
    for index in six.moves.xrange(1, num_nodes - 1):  # Only interior nodes.
        curr_dist = (
            coeff_a * nodes[0, index] + coeff_b * nodes[1, index] + coeff_c)
        if curr_dist < d_min:
            d_min = curr_dist
        elif curr_dist > d_max:
            d_max = curr_dist

    return coeff_a, coeff_b, coeff_c, d_min, d_max
