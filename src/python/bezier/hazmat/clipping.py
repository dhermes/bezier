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
.. _quadratic convergence: https://doi.org/10.1016/j.cagd.2007.12.006

The B |eacute| zier clipping `algorithm`_ is used to intersect
two planar B |eacute| zier curves. It proceeds by using "fat lines"
to recursively prune the region of accepted parameter ranges until
the ranges converge to points. (A "fat line" is a rectangular region of a
bounded distance from the line connecting the start and end points of a
B |eacute| zier curve.)

It has `quadratic convergence`_. It can be used to find tangent intersections,
which is the primary usage within ``bezier``.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import numpy as np

from bezier import _helpers
from bezier.hazmat import geometric_intersection


NO_PARALLEL = "Parallel lines not supported during clipping."
DEFAULT_S_MIN = 1.0
DEFAULT_S_MAX = 0.0


def compute_implicit_line(nodes):
    r"""Compute the implicit form of the line connecting curve endpoints.

    .. note::

       This assumes, but does not check, that the first and last node
       in ``nodes`` are different.

    Computes :math:`a, b` and :math:`c` in the implicit equation
    for the line

    .. math::

       ax + by + c = 0

    In cases where the line is normalized (:math:`a^2 + b^2 = 1`, only unique
    up to sign), the function :math:`d(x, y) = ax + by + c` can be used to
    determine the point-line distance from a point
    :math:`\left[\begin{array}{c} x \\ y \end{array}\right]` to the line.

    .. image:: ../../images/compute_implicit_line.png
       :align: center

    .. testsetup:: compute-implicit-line

       import numpy as np
       from bezier.hazmat.clipping import compute_implicit_line

    .. doctest:: compute-implicit-line
       :options: +NORMALIZE_WHITESPACE

       >>> import numpy as np
       >>> nodes = np.asfortranarray([
       ...     [0.0, 1.0, 3.0, 4.0],
       ...     [0.0, 2.5, 0.5, 3.0],
       ... ])
       >>> coeff_a, coeff_b, coeff_c = compute_implicit_line(nodes)
       >>> float(coeff_a), float(coeff_b), float(coeff_c)
       (-3.0, 4.0, 0.0)

    .. testcleanup:: compute-implicit-line

       import make_images
       make_images.compute_implicit_line(nodes)

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
    # Rotate 90 degrees to the "left".
    coeff_a = -delta[1]
    coeff_b = delta[0]
    # c = - ax - by = delta[1] x - delta[0] y
    coeff_c = delta[1] * nodes[0, 0] - delta[0] * nodes[1, 0]
    return coeff_a, coeff_b, coeff_c


def compute_fat_line(nodes):
    """Compute the "fat line" around a B |eacute| zier curve.

    Both computes the implicit form

    .. math::

       ax + by + c = 0

    for the line connecting the first and last node in ``nodes``.
    Also computes the maximum and minimum (non-normalized) distances to that
    line from each control point where (non-normalized) distance :math:`d` is
    computed as :math:`d_i = a x_i + b y_i + c`. (To make :math:`d` a true
    measure of Euclidean distance the line would need to be normalized so that
    :math:`a^2 + b^2 = 1`.)

    .. image:: ../../images/compute_fat_line.png
       :align: center

    .. testsetup:: compute-fat-line

       import numpy as np
       from bezier.hazmat.clipping import compute_fat_line

    .. doctest:: compute-fat-line
       :options: +NORMALIZE_WHITESPACE

       >>> nodes = np.asfortranarray([
       ...     [0.0, 1.0, 3.0, 4.0],
       ...     [2.0, 4.5, 2.5, 5.0],
       ... ])
       >>> coeff_a, coeff_b, coeff_c, d_min, d_max = compute_fat_line(nodes)
       >>> float(coeff_a), float(coeff_b), float(coeff_c)
       (-3.0, 4.0, -8.0)
       >>> float(d_min), float(d_max)
       (-7.0, 7.0)

    .. testcleanup:: compute-fat-line

       import make_images
       info = coeff_a, coeff_b, coeff_c, d_min, d_max
       make_images.compute_fat_line(nodes, info)

    Args:
        nodes (numpy.ndarray): ``2 x N`` array of nodes in a curve.

    Returns:
        Tuple[float, float, float, float, float]: The 5-tuple of

        * The :math:`x` coefficient :math:`a`
        * The :math:`y` coefficient :math:`b`
        * The constant :math:`c`
        * The "minimum" (non-normalized) distance to the fat line among the
          control points.
        * The "maximum" (non-normalized) distance to the fat line among the
          control points.
    """
    coeff_a, coeff_b, coeff_c = compute_implicit_line(nodes)
    # NOTE: This assumes, but does not check, that there are two rows.
    _, num_nodes = nodes.shape
    d_min = 0.0
    d_max = 0.0
    for index in range(1, num_nodes - 1):  # Only interior nodes.
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
    s, t, success = geometric_intersection.segment_intersection(
        start0, end0, start1, end1
    )

    if not success:
        raise NotImplementedError(NO_PARALLEL)

    # NOTE: We can only **widen** the interval with a real intersection.
    #       I.e. we can push the minimum to the left or the maximum to the
    #       right.

    if _helpers.in_interval(t, 0.0, 1.0):
        if _helpers.in_interval(s, 0.0, s_min):
            s_min = s

        if _helpers.in_interval(s, s_max, 1.0):
            s_max = s

    return s_min, s_max


def _clip_range_polynomial(nodes, coeff_a, coeff_b, coeff_c):
    r"""Compute control points for a polynomial used to clip range.

    Args:
        nodes (numpy.ndarray): ``2 x N`` array of nodes in a curve.
            The line will be (directed) from the first to last
            node in ``nodes``.
        coeff_a (float): The :math:`a` coefficient in a line
            :math:`ax + by + c = 0`.
        coeff_b (float): The :math:`b` coefficient in a line
            :math:`ax + by + c = 0`.
        coeff_c (float): The :math:`c` coefficient in a line
            :math:`ax + by + c = 0`.

    Returns:
        Tuple[int, numpy.ndarray]: Pair of

        * The degree :math:`d` of ``nodes`` (also the width of the
          ``x``-interval :math:`\left[0, d\right]`)
        * ``2 x N`` array of polynomial curve with distances
          :math:`d_i = a x_i + b y_i + c` as the control points (and the
          ``x``-coordinates evenly spaced).
    """
    _, num_nodes = nodes.shape
    polynomial = np.empty((2, num_nodes), order="F")
    for index in range(num_nodes):
        polynomial[0, index] = index
        polynomial[1, index] = (
            coeff_a * nodes[0, index] + coeff_b * nodes[1, index] + coeff_c
        )

    return num_nodes - 1, polynomial


def clip_range(nodes1, nodes2):
    r"""Reduce the parameter range where two curves can intersect.

    .. note::

       This assumes, but does not check that the curves being considered
       will only have one intersection in the parameter ranges
       :math:`s \in \left[0, 1\right]`, :math:`t \in \left[0, 1\right]`.
       This assumption is based on the fact that B |eacute| zier clipping
       is meant to be used (within this library) to find tangent intersections
       for already subdivided (i.e. sufficiently zoomed in) curve segments.

    Two B |eacute| zier curves :math:`B_1(s)` and :math:`B_2(t)` are defined by
    ``nodes1`` and ``nodes2``. The "fat line" (see :func:`compute_fat_line`)
    for :math:`B_1(s)` is used to narrow the range of possible :math:`t`-values
    in an intersection by considering the (non-normalized) distance polynomial
    for :math:`B_2(t)`:

    .. math::

       d(t) = \sum_{j = 0}^m \binom{m}{j} t^j (1 - t)^{m - j} \cdot d_j

    Here :math:`d_j = a x_j + b y_j + c` are the (non-normalized) distances of
    each control point
    :math:`\left[\begin{array}{c} x_j \\ y_j \end{array}\right]`
    of :math:`B_2(t)` to the implicit line for :math:`B_1(s)`.

    Consider the following pair of B |eacute| zier curves and the distances
    from **all** of the control points to the implicit line for :math:`B_1(s)`:

    .. testsetup:: clip-range-start, clip-range

       import numpy as np
       from bezier.hazmat.clipping import clip_range

       nodes1 = np.asfortranarray([
           [2.0, 4.5, 2.5, 5.0],
           [0.0, 1.0, 3.0, 4.0],
       ])
       nodes2 = np.asfortranarray([
           [-0.25 , 3.75 , 7.0  ],
           [ 3.125, 0.875, 3.125],
       ])

    .. doctest:: clip-range-start
       :options: +NORMALIZE_WHITESPACE

       >>> nodes1
       array([[2. , 4.5, 2.5, 5. ],
              [0. , 1. , 3. , 4. ]])
       >>> nodes2
       array([[-0.25 ,  3.75 ,  7.   ],
              [ 3.125,  0.875,  3.125]])

    .. image:: ../../images/clip_range.png
       :align: center

    The (non-normalized) distances from the control points of :math:`B_2(t)`
    define the distance polynomial :math:`d(t)`. By writing this polynomial as
    a B |eacute| zier curve, a convex hull can be formed. The intersection of
    this convex hull with the "fat line" of :math:`B_1(s)` determines the
    extreme :math:`t` values possible and allows clipping the range of
    :math:`B_2(t)`:

    .. image:: ../../images/clip_range_distances.png
       :align: center

    .. doctest:: clip-range
       :options: +NORMALIZE_WHITESPACE

       >>> s_min, s_max = clip_range(nodes1, nodes2)
       >>> float(s_min), float(s_max)
       (0.25, 0.875)

    .. testcleanup:: clip-range

       import make_images
       make_images.clip_range(nodes1, nodes2)
       make_images.clip_range_distances(nodes1, nodes2)

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
    # NOTE: There is no corresponding "pylint: enable", but the disable only
    #       applies in this lexical scope.
    # pylint: disable=too-many-locals
    coeff_a, coeff_b, coeff_c, d_min, d_max = compute_fat_line(nodes1)
    degree2, polynomial = _clip_range_polynomial(
        nodes2, coeff_a, coeff_b, coeff_c
    )
    # Define segments for the top and the bottom of the region
    # bounded by the fat line.
    start_bottom = np.asfortranarray([0.0, d_min])
    end_bottom = np.asfortranarray([degree2, d_min])
    start_top = np.asfortranarray([0.0, d_max])
    end_top = np.asfortranarray([degree2, d_max])

    s_min = DEFAULT_S_MIN
    s_max = DEFAULT_S_MAX
    # NOTE: We avoid computing the convex hull of ``d(t)`` / ``polynomial`` and
    #       just compute where all segments connecting two control points
    #       intersect the fat lines. In order to account for intersections
    #       at ``t = 0`` or ``t = 1``, we just check the height of ``d(t)``.
    if d_min <= polynomial[1, 0] <= d_max:
        s_min = 0.0  # NOT YET
    if d_min <= polynomial[1, degree2] <= d_max:
        s_max = 1.0  # NOT YET
    for start_index in range(degree2):
        for end_index in range(start_index + 1, degree2 + 1):
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
    return s_min, s_max
