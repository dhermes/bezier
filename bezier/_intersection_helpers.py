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

"""Private helper methods for intersecting B |eacute| zier shapes.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


import numpy as np


_WEIGHTS_QUADRATIC = np.array([
    [1.0, -2.0, 1.0],
])
_WEIGHTS_CUBIC = np.array([
    [1.0, -2.0, 1.0, 0.0],
    [0.0, 1.0, -2.0, 1.0],
])
_WEIGHTS_QUARTIC = np.array([
    [1.0, -2.0, 1.0, 0.0, 0.0],
    [0.0, 1.0, -2.0, 1.0, 0.0],
    [0.0, 0.0, 1.0, -2.0, 1.0],
])


def bbox_intersect(nodes1, nodes2):
    r"""Bounding box intersection predicate.

    Determines if the bounding box of two sets of control points
    intersects in :math:`\mathbf{R}^2` with non-trivial
    intersection (i.e. tangent bounding boxes are insufficient).

    .. note::

        Though we assume (and the code relies on this fact) that
        the nodes are two-dimensional, we don't check it.

    Args:
        nodes1 (numpy.ndarray): Set of control points for a
            B |eacute| zier shape.
        nodes2 (numpy.ndarray): Set of control points for a
            B |eacute| zier shape.

    Returns:
        bool: Predicate indicating if the bounding boxes intersect.
    """
    left1, bottom1 = np.min(nodes1, axis=0)
    right1, top1 = np.max(nodes1, axis=0)
    left2, bottom2 = np.min(nodes2, axis=0)
    right2, top2 = np.max(nodes2, axis=0)

    return (right2 > left1 and right1 > left2 and
            top2 > bottom1 and top1 > bottom2)


def linearization_error(curve):
    r"""Compute the maximum error of a linear approximation.

    We use the line

    .. math::

       L(s) = v_0 (1 - s) + v_n s

    and compute a bound on the maximum error

    .. math::

       \max_{s \in \left[0, 1\right]} \|B(s) - L(s)\|.

    Rather than computing the actual maximum (a tight bound), we
    use an upper bound via the remainder from Lagrange interpolation
    in each component. This leaves us with :math:`\frac{s(s - 1)}{2!}`
    times the second derivative in each component.

    The second derivative curve is degree :math:`d = n - 2` and
    is given by

    .. math::

       B''(s) = n(n - 1) \sum_{j = 0}^{d} \binom{d}{j} s^j
       (1 - s)^{d - j} \cdot \Delta^2 v_j

    Due to this form (and the convex combination property of
    B |eacute| zier Curves) we know each component of the second derivative
    will be bounded by the maximum of that component among the
    :math:`\Delta^2 v_j`.

    Args:
        curve (~bezier.curve.Curve): A curve to be approximated by a line.

    Returns:
        float: The maximum error between the curve and the
        linear approximation.
    """
    if curve.degree == 1:
        return 0.0

    if curve.degree == 2:
        weights = _WEIGHTS_QUADRATIC
    elif curve.degree == 3:
        weights = _WEIGHTS_CUBIC
    elif curve.degree == 4:
        weights = _WEIGHTS_QUARTIC
    else:
        weights = np.zeros((curve.degree - 1, curve.degree + 1))
        eye = np.eye(curve.degree - 1)
        weights[:, :-2] = eye
        weights[:, 1:-1] += -2 * eye
        weights[:, 2:] += eye

    nodes = curve._nodes  # pylint: disable=protected-access
    second_deriv = weights.dot(nodes)  # pylint: disable=no-member
    worst_case = np.max(np.abs(second_deriv), axis=0)

    # max_{0 <= s <= 1} s(1 - s)/2 = 1/8 = 0.125
    multiplier = 0.125 * curve.degree * (curve.degree - 1)
    return multiplier * np.linalg.norm(worst_case, ord=2)
