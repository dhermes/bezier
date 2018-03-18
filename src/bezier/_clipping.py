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

from bezier import _geometric_intersection
from bezier import _helpers

NO_PARALLEL = 'Parallel lines not supported during clipping.'
DEFAULT_S_MIN = 1.0
DEFAULT_S_MAX = 0.0


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
            coeff_a * nodes[0, index] + coeff_b * nodes[1, index] + coeff_c
        )
        if curr_dist < d_min:
            d_min = curr_dist
        elif curr_dist > d_max:
            d_max = curr_dist
    return coeff_a, coeff_b, coeff_c, d_min, d_max


def _update_parameters(s_min, s_max, start0, end0, start1, end1):
    """Update clipped parameter range.

    .. note::

       This is a helper for :func:`clip_range`.

    Does so by intersecting one of the two fat lines with an edge
    of the convex hull of the distance polynomial of the curve being
    clipped.

    If both of ``s_min`` and ``s_max`` are "unset", then any :math:`s`
    value that is valid for ``s_min`` would also be valid for ``s_max``.
    Rather than adding a special case to handle this scenario, **only**
    ``s_min`` will be updated.

    In cases where a given parameter :math:`s` would be a valid update
    for both ``s_min`` and ``s_max``
    This function **only** updates ``s_min``

    Args:
        s_min (float): Current start of clipped interval. If "unset", this
            value will be ``DEFAULT_S_MIN``.
        s_max (float): Current end of clipped interval. If "unset", this
            value will be ``DEFAULT_S_MAX``.
        start0 (numpy.ndarray): A 1D NumPy ``2``-array that is the start
            vector of one of the two fat lines.
        end0 (numpy.ndarray): A 1D NumPy ``2``-array that is the end
            vector of one of the two fat lines.
        start1 (numpy.ndarray): A 1D NumPy ``2``-array that is the start
            vector of an edge of the convex hull of the distance
            polynomial :math:`d(t)` as an explicit B |eacute| zier curve.
        end1 (numpy.ndarray): A 1D NumPy ``2``-array that is the end
            vector of an edge of the convex hull of the distance
            polynomial :math:`d(t)` as an explicit B |eacute| zier curve.

    Returns:
        Tuple[float, float]: The (possibly updated) start and end
        of the clipped parameter range.

    Raises:
        NotImplementedError: If the two line segments are parallel. (This
            case will be supported at some point, just not now.)
    """
    s, t, success = _geometric_intersection.segment_intersection(
        start0, end0, start1, end1
    )
    if not success:
        raise NotImplementedError(NO_PARALLEL)

    if _helpers.in_interval(t, 0.0, 1.0):
        if _helpers.in_interval(s, 0.0, s_min):
            return s, s_max

        elif _helpers.in_interval(s, s_max, 1.0):
            return s_min, s

    return s_min, s_max


def _check_parameter_range(s_min, s_max):
    r"""Performs a final check on a clipped parameter range.

    .. note::

       This is a helper for :func:`clip_range`.

    If both values are unchanged from the "unset" default, this returns
    the whole interval :math:`\left[0.0, 1.0\right]`.

    If only one of the values is set to some parameter :math:`s`, this
    returns the "degenerate" interval :math:`\left[s, s\right]`. (We rely
    on the fact that ``s_min`` must be the only set value, based on how
    :func:`_update_parameters` works.)

    Otherwise, this simply returns ``[s_min, s_max]``.

    Args:
        s_min (float): Current start of clipped interval. If "unset", this
            value will be ``DEFAULT_S_MIN``.
        s_max (float): Current end of clipped interval. If "unset", this
            value will be ``DEFAULT_S_MAX``.

    Returns:
        Tuple[float, float]: The (possibly updated) start and end
        of the clipped parameter range.
    """
    if s_min == DEFAULT_S_MIN:
        # Based on the way ``_update_parameters`` works, we know
        # both parameters must be unset if ``s_min``.
        return 0.0, 1.0

    if s_max == DEFAULT_S_MAX:
        return s_min, s_min

    return s_min, s_max


def clip_range(nodes1, nodes2):
    r"""Reduce the parameter range where two curves can intersect.

    Does so by using the "fat line" for ``nodes1`` and computing the
    distance polynomial against ``nodes2``.

    .. note::

       This assumes, but does not check that the curves being considered
       will only have one intersection in the parameter ranges
       :math:`s \in \left[0, 1\right]`, :math:`t \in \left[0, 1\right]`.
       This assumption is based on the fact that B |eacute| zier clipping
       is meant to be used to find tangent intersections for already
       subdivided (i.e. sufficiently zoomed in) curve segments.

    Args:
        nodes1 (numpy.ndarray): ``2 x N1`` array of nodes in a curve which
            will define the clipping region.
        nodes2 (numpy.ndarray): ``2 x N2`` array of nodes in a curve which
            will be clipped.

    Returns:
        Tuple[float, float]: The pair of

        * The start parameter of the clipped range.
        * The end parameter of the clipped range.
    """
    # pylint: disable=too-many-locals
    coeff_a, coeff_b, coeff_c, d_min, d_max = compute_fat_line(nodes1)
    # NOTE: This assumes, but does not check, that there are two rows.
    _, num_nodes2 = nodes2.shape
    polynomial = np.empty((2, num_nodes2), order='F')
    denominator = float(num_nodes2 - 1)
    for index in six.moves.xrange(num_nodes2):
        polynomial[0, index] = index / denominator
        polynomial[1, index] = (
            coeff_a * nodes2[0, index] + coeff_b * nodes2[1, index] + coeff_c
        )
    # Define segments for the top and the bottom of the region
    # bounded by the fat line.
    start_bottom = np.asfortranarray([0.0, d_min])
    end_bottom = np.asfortranarray([1.0, d_min])
    start_top = np.asfortranarray([0.0, d_max])
    end_top = np.asfortranarray([1.0, d_max])
    s_min = DEFAULT_S_MIN
    s_max = DEFAULT_S_MAX
    # NOTE: We avoid computing the convex hull and just compute where
    #       all segments connecting two control points intersect the
    #       fat lines.
    for start_index in six.moves.xrange(num_nodes2 - 1):
        for end_index in six.moves.xrange(start_index + 1, num_nodes2):
            s_min, s_max = _update_parameters(
                s_min,
                s_max,
                start_bottom,
                end_bottom,
                polynomial[:, start_index],
                polynomial[:, end_index],
            )
            s_min, s_max = _update_parameters(
                s_min,
                s_max,
                start_top,
                end_top,
                polynomial[:, start_index],
                polynomial[:, end_index],
            )
    return _check_parameter_range(s_min, s_max)
    # pylint: enable=too-many-locals
