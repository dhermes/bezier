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

"""Helpers for intersecting B |eacute| zier curves via geometric methods.

The functions here are pure Python and many have equivalent implementations
written in Fortran and exposed via a Cython wrapper.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import itertools

import numpy as np

from bezier.hazmat import curve_helpers
from bezier.hazmat import helpers as _py_helpers
from bezier.hazmat import intersection_helpers


# Set the threshold for exponent at half the bits available, this way one round
# of Newton's method can (usually) finish the job by squaring the error.
_ERROR_VAL = 0.5**26
_MAX_INTERSECT_SUBDIVISIONS = 20
_MAX_CANDIDATES = 64
_UNHANDLED_LINES = (
    "If both curves are lines, the intersection should have "
    "been computed already."
)
_TOO_MANY_TEMPLATE = (
    "The number of candidate intersections is too high.\n"
    "{:d} candidate pairs."
)
_NO_CONVERGE_TEMPLATE = (
    "Curve intersection failed to converge to approximately linear "
    "subdivisions after {:d} iterations."
)
_MIN_INTERVAL_WIDTH = 0.5**40
_SPLIT_POINT = np.asfortranarray([[0.5], [0.5]])


def bbox_intersect(nodes1, nodes2):
    r"""Bounding box intersection predicate.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

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
        int: Enum from :class:`.BoxIntersectionType` indicating the type of
        bounding box intersection.
    """
    left1, right1, bottom1, top1 = _py_helpers.bbox(nodes1)
    left2, right2, bottom2, top2 = _py_helpers.bbox(nodes2)
    if right2 < left1 or right1 < left2 or top2 < bottom1 or top1 < bottom2:
        return BoxIntersectionType.DISJOINT

    if (
        right2 == left1
        or right1 == left2
        or top2 == bottom1
        or top1 == bottom2
    ):
        return BoxIntersectionType.TANGENT

    else:
        return BoxIntersectionType.INTERSECTION


def linearization_error(nodes):
    r"""Compute the maximum error of a linear approximation.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    .. note::

       This is a helper for :class:`.Linearization`, which is used during the
       curve-curve intersection process.

    We use the line

    .. math::

       L(s) = v_0 (1 - s) + v_n s

    and compute a bound on the maximum error

    .. math::

       \max_{s \in \left[0, 1\right]} \|B(s) - L(s)\|_2.

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

    For example, the curve

    .. math::

       B(s) = \left[\begin{array}{c} 0 \\ 0 \end{array}\right] (1 - s)^2
           + \left[\begin{array}{c} 3 \\ 1 \end{array}\right] 2s(1 - s)
           + \left[\begin{array}{c} 9 \\ -2 \end{array}\right] s^2

    has
    :math:`B''(s) \equiv \left[\begin{array}{c} 6 \\ -8 \end{array}\right]`
    which has norm :math:`10` everywhere, hence the maximum error is

    .. math::

       \left.\frac{s(1 - s)}{2!} \cdot 10\right|_{s = \frac{1}{2}}
          = \frac{5}{4}.

    .. image:: ../../images/linearization_error.png
       :align: center

    .. testsetup:: linearization-error, linearization-error-fail

       import numpy as np
       from bezier.hazmat.geometric_intersection import linearization_error

    .. doctest:: linearization-error

       >>> import numpy as np
       >>> nodes = np.asfortranarray([
       ...     [0.0, 3.0, 9.0],
       ...     [0.0, 1.0, -2.0],
       ... ])
       >>> float(linearization_error(nodes))
       1.25

    .. testcleanup:: linearization-error

       import make_images
       make_images.linearization_error(nodes)

    As a **non-example**, consider a "pathological" set of control points:

    .. math::

       B(s) = \left[\begin{array}{c} 0 \\ 0 \end{array}\right] (1 - s)^3
           + \left[\begin{array}{c} 5 \\ 12 \end{array}\right] 3s(1 - s)^2
           + \left[\begin{array}{c} 10 \\ 24 \end{array}\right] 3s^2(1 - s)
           + \left[\begin{array}{c} 30 \\ 72 \end{array}\right] s^3

    By construction, this lies on the line :math:`y = \frac{12x}{5}`, but
    the parametrization is cubic:
    :math:`12 \cdot x(s) = 5 \cdot y(s) = 180s(s^2 + 1)`. Hence, the fact
    that the curve is a line is not accounted for and we take the worse
    case among the nodes in:

    .. math::

       B''(s) = 3 \cdot 2 \cdot \left(
           \left[\begin{array}{c} 0 \\ 0 \end{array}\right] (1 - s)
           + \left[\begin{array}{c} 15 \\ 36 \end{array}\right] s\right)

    which gives a nonzero maximum error:

    .. doctest:: linearization-error-fail

       >>> nodes = np.asfortranarray([
       ...     [0.0, 5.0, 10.0, 30.0],
       ...     [0.0, 12.0, 24.0, 72.0],
       ... ])
       >>> float(linearization_error(nodes))
       29.25

    Though it may seem that ``0`` is a more appropriate answer, consider
    the **goal** of this function. We seek to linearize curves and then
    intersect the linear approximations. Then the :math:`s`-values from
    the line-line intersection is lifted back to the curves. Thus
    the error :math:`\|B(s) - L(s)\|_2` is more relevant than the
    underlying algebraic curve containing :math:`B(s)`.

    .. note::

       It may be more appropriate to use a **relative** linearization error
       rather than the **absolute** error provided here. It's unclear if
       the domain :math:`\left[0, 1\right]` means the error is **already**
       adequately scaled or if the error should be scaled by the arc
       length of the curve or the (easier-to-compute) length of the line.

    Args:
        nodes (numpy.ndarray): Nodes of a curve.

    Returns:
        float: The maximum error between the curve and the
        linear approximation.
    """
    _, num_nodes = nodes.shape
    degree = num_nodes - 1
    if degree == 1:
        return 0.0

    second_deriv = nodes[:, :-2] - 2.0 * nodes[:, 1:-1] + nodes[:, 2:]
    worst_case = np.max(np.abs(second_deriv), axis=1)
    # max_{0 <= s <= 1} s(1 - s)/2 = 1/8 = 0.125
    multiplier = 0.125 * degree * (degree - 1)
    # NOTE: worst_case is 1D due to np.max(), so this is the vector norm.
    return multiplier * np.linalg.norm(worst_case, ord=2)


def segment_intersection(start0, end0, start1, end1):
    r"""Determine the intersection of two line segments.

    Assumes each line is parametric

    .. math::

       \begin{alignat*}{2}
        L_0(s) &= S_0 (1 - s) + E_0 s &&= S_0 + s \Delta_0 \\
        L_1(t) &= S_1 (1 - t) + E_1 t &&= S_1 + t \Delta_1.
       \end{alignat*}

    To solve :math:`S_0 + s \Delta_0 = S_1 + t \Delta_1`, we use the
    cross product:

    .. math::

       \left(S_0 + s \Delta_0\right) \times \Delta_1 =
           \left(S_1 + t \Delta_1\right) \times \Delta_1 \Longrightarrow
       s \left(\Delta_0 \times \Delta_1\right) =
           \left(S_1 - S_0\right) \times \Delta_1.

    Similarly

    .. math::

       \Delta_0 \times \left(S_0 + s \Delta_0\right) =
           \Delta_0 \times \left(S_1 + t \Delta_1\right) \Longrightarrow
       \left(S_1 - S_0\right) \times \Delta_0 =
           \Delta_0 \times \left(S_0 - S_1\right) =
           t \left(\Delta_0 \times \Delta_1\right).

    .. note::

       Since our points are in :math:`\mathbf{R}^2`, the "traditional"
       cross product in :math:`\mathbf{R}^3` will always point in the
       :math:`z` direction, so in the above we mean the :math:`z`
       component of the cross product, rather than the entire vector.

    For example, the diagonal lines

    .. math::

       \begin{align*}
        L_0(s) &= \left[\begin{array}{c} 0 \\ 0 \end{array}\right] (1 - s) +
                  \left[\begin{array}{c} 2 \\ 2 \end{array}\right] s \\
        L_1(t) &= \left[\begin{array}{c} -1 \\ 2 \end{array}\right] (1 - t) +
                  \left[\begin{array}{c} 1 \\ 0 \end{array}\right] t
       \end{align*}

    intersect at :math:`L_0\left(\frac{1}{4}\right) =
    L_1\left(\frac{3}{4}\right) =
    \frac{1}{2} \left[\begin{array}{c} 1 \\ 1 \end{array}\right]`.

    .. image:: ../../images/segment_intersection1.png
       :align: center

    .. testsetup:: segment-intersection1, segment-intersection2

       import numpy as np
       from bezier.hazmat.geometric_intersection import segment_intersection

    .. doctest:: segment-intersection1
       :options: +NORMALIZE_WHITESPACE

       >>> start0 = np.asfortranarray([0.0, 0.0])
       >>> end0 = np.asfortranarray([2.0, 2.0])
       >>> start1 = np.asfortranarray([-1.0, 2.0])
       >>> end1 = np.asfortranarray([1.0, 0.0])
       >>> s, t, _ = segment_intersection(start0, end0, start1, end1)
       >>> float(s)
       0.25
       >>> float(t)
       0.75

    .. testcleanup:: segment-intersection1

       import make_images
       make_images.segment_intersection1(start0, end0, start1, end1, s)

    Taking the parallel (but different) lines

    .. math::

       \begin{align*}
        L_0(s) &= \left[\begin{array}{c} 1 \\ 0 \end{array}\right] (1 - s) +
                  \left[\begin{array}{c} 0 \\ 1 \end{array}\right] s \\
        L_1(t) &= \left[\begin{array}{c} -1 \\ 3 \end{array}\right] (1 - t) +
                  \left[\begin{array}{c} 3 \\ -1 \end{array}\right] t
       \end{align*}

    we should be able to determine that the lines don't intersect, but
    this function is not meant for that check:

    .. image:: ../../images/segment_intersection2.png
       :align: center

    .. doctest:: segment-intersection2
       :options: +NORMALIZE_WHITESPACE

       >>> start0 = np.asfortranarray([1.0, 0.0])
       >>> end0 = np.asfortranarray([0.0, 1.0])
       >>> start1 = np.asfortranarray([-1.0, 3.0])
       >>> end1 = np.asfortranarray([3.0, -1.0])
       >>> _, _, success = segment_intersection(start0, end0, start1, end1)
       >>> success
       False

    .. testcleanup:: segment-intersection2

       import make_images
       make_images.segment_intersection2(start0, end0, start1, end1)

    Instead, we use :func:`parallel_lines_parameters`:

    .. testsetup:: segment-intersection2-continued

       import numpy as np
       from bezier.hazmat.geometric_intersection import (
           parallel_lines_parameters
       )

       start0 = np.asfortranarray([1.0, 0.0])
       end0 = np.asfortranarray([0.0, 1.0])
       start1 = np.asfortranarray([-1.0, 3.0])
       end1 = np.asfortranarray([3.0, -1.0])

    .. doctest:: segment-intersection2-continued

       >>> disjoint, _ = parallel_lines_parameters(start0, end0, start1, end1)
       >>> disjoint
       True

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        start0 (numpy.ndarray): A 1D NumPy ``2``-array that is the start
            vector :math:`S_0` of the parametric line :math:`L_0(s)`.
        end0 (numpy.ndarray): A 1D NumPy ``2``-array that is the end
            vector :math:`E_0` of the parametric line :math:`L_0(s)`.
        start1 (numpy.ndarray): A 1D NumPy ``2``-array that is the start
            vector :math:`S_1` of the parametric line :math:`L_1(s)`.
        end1 (numpy.ndarray): A 1D NumPy ``2``-array that is the end
            vector :math:`E_1` of the parametric line :math:`L_1(s)`.

    Returns:
        Tuple[float, float, bool]: Pair of :math:`s_{\ast}` and
        :math:`t_{\ast}` such that the lines intersect:
        :math:`L_0\left(s_{\ast}\right) = L_1\left(t_{\ast}\right)` and then
        a boolean indicating if an intersection was found (i.e. if the lines
        aren't parallel).
    """
    delta0 = end0 - start0
    delta1 = end1 - start1
    cross_d0_d1 = _py_helpers.cross_product(delta0, delta1)
    if cross_d0_d1 == 0.0:
        return None, None, False

    else:
        start_delta = start1 - start0
        s = _py_helpers.cross_product(start_delta, delta1) / cross_d0_d1
        t = _py_helpers.cross_product(start_delta, delta0) / cross_d0_d1
        return s, t, True


def parallel_lines_parameters(start0, end0, start1, end1):
    r"""Checks if two parallel lines ever meet.

    Meant as a back-up when :func:`segment_intersection` fails.

    .. note::

       This function assumes but never verifies that the lines
       are parallel.

    In the case that the segments are parallel and lie on **different**
    lines, then there is a **guarantee** of no intersection. However, if
    they are on the exact same line, they may define a shared segment
    coincident to both lines.

    In :func:`segment_intersection`, we utilized the normal form of the
    lines (via the cross product):

    .. math::

       \begin{align*}
       L_0(s) \times \Delta_0 &\equiv S_0 \times \Delta_0 \\
       L_1(t) \times \Delta_1 &\equiv S_1 \times \Delta_1
       \end{align*}

    So, we can detect if :math:`S_1` is on the first line by
    checking if

    .. math::

       S_0 \times \Delta_0 \stackrel{?}{=} S_1 \times \Delta_0.

    If it is not on the first line, then we are done, the
    segments don't meet:

    .. image:: ../../images/parallel_lines_parameters1.png
       :align: center

    .. testsetup:: parallel-different1, parallel-different2

       import numpy as np
       from bezier.hazmat.geometric_intersection import (
           parallel_lines_parameters
       )

    .. doctest:: parallel-different1

       >>> # Line: y = 1
       >>> start0 = np.asfortranarray([0.0, 1.0])
       >>> end0 = np.asfortranarray([1.0, 1.0])
       >>> # Vertical shift up: y = 2
       >>> start1 = np.asfortranarray([-1.0, 2.0])
       >>> end1 = np.asfortranarray([3.0, 2.0])
       >>> disjoint, _ = parallel_lines_parameters(start0, end0, start1, end1)
       >>> disjoint
       True

    .. testcleanup:: parallel-different1

       import make_images
       make_images.helper_parallel_lines(
           start0, end0, start1, end1, "parallel_lines_parameters1.png")

    If :math:`S_1` **is** on the first line, we want to check that
    :math:`S_1` and :math:`E_1` define parameters outside of
    :math:`\left[0, 1\right]`. To compute these parameters:

    .. math::

       L_1(t) = S_0 + s_{\ast} \Delta_0 \Longrightarrow
           s_{\ast} = \frac{\Delta_0^T \left(
               L_1(t) - S_0\right)}{\Delta_0^T \Delta_0}.

    For example, the intervals :math:`\left[0, 1\right]` and
    :math:`\left[\frac{3}{2}, 2\right]` (via
    :math:`S_1 = S_0 + \frac{3}{2} \Delta_0` and
    :math:`E_1 = S_0 + 2 \Delta_0`) correspond to segments that
    don't meet:

    .. image:: ../../images/parallel_lines_parameters2.png
       :align: center

    .. doctest:: parallel-different2

       >>> start0 = np.asfortranarray([1.0, 0.0])
       >>> delta0 = np.asfortranarray([2.0, -1.0])
       >>> end0 = start0 + 1.0 * delta0
       >>> start1 = start0 + 1.5 * delta0
       >>> end1 = start0 + 2.0 * delta0
       >>> disjoint, _ = parallel_lines_parameters(start0, end0, start1, end1)
       >>> disjoint
       True

    .. testcleanup:: parallel-different2

       import make_images
       make_images.helper_parallel_lines(
           start0, end0, start1, end1, "parallel_lines_parameters2.png")

    but if the intervals overlap, like :math:`\left[0, 1\right]` and
    :math:`\left[-1, \frac{1}{2}\right]`, the segments meet:

    .. image:: ../../images/parallel_lines_parameters3.png
       :align: center

    .. testsetup:: parallel-different3, parallel-different4

       import numpy as np
       from bezier.hazmat.geometric_intersection import (
           parallel_lines_parameters
       )

       start0 = np.asfortranarray([1.0, 0.0])
       delta0 = np.asfortranarray([2.0, -1.0])
       end0 = start0 + 1.0 * delta0

    .. doctest:: parallel-different3

       >>> start1 = start0 - 1.5 * delta0
       >>> end1 = start0 + 0.5 * delta0
       >>> disjoint, parameters = parallel_lines_parameters(
       ...     start0, end0, start1, end1)
       >>> disjoint
       False
       >>> parameters
       array([[0.  , 0.5 ],
              [0.75, 1.  ]])

    .. testcleanup:: parallel-different3

       import make_images
       make_images.helper_parallel_lines(
           start0, end0, start1, end1, "parallel_lines_parameters3.png")

    Similarly, if the second interval completely contains the first,
    the segments meet:

    .. image:: ../../images/parallel_lines_parameters4.png
       :align: center

    .. doctest:: parallel-different4

       >>> start1 = start0 + 4.5 * delta0
       >>> end1 = start0 - 3.5 * delta0
       >>> disjoint, parameters = parallel_lines_parameters(
       ...     start0, end0, start1, end1)
       >>> disjoint
       False
       >>> parameters
       array([[1.    , 0.    ],
              [0.4375, 0.5625]])

    .. testcleanup:: parallel-different4

       import make_images
       make_images.helper_parallel_lines(
           start0, end0, start1, end1, "parallel_lines_parameters4.png")

    .. note::

       This function doesn't currently allow wiggle room around the
       desired value, i.e. the two values must be bitwise identical.
       However, the most "correct" version of this function likely
       should allow for some round off.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        start0 (numpy.ndarray): A 1D NumPy ``2``-array that is the start
            vector :math:`S_0` of the parametric line :math:`L_0(s)`.
        end0 (numpy.ndarray): A 1D NumPy ``2``-array that is the end
            vector :math:`E_0` of the parametric line :math:`L_0(s)`.
        start1 (numpy.ndarray): A 1D NumPy ``2``-array that is the start
            vector :math:`S_1` of the parametric line :math:`L_1(s)`.
        end1 (numpy.ndarray): A 1D NumPy ``2``-array that is the end
            vector :math:`E_1` of the parametric line :math:`L_1(s)`.

    Returns:
        Tuple[bool, Optional[numpy.ndarray]]: A pair of

        * Flag indicating if the lines are disjoint.
        * An optional ``2 x 2`` matrix of ``s-t`` parameters only present if
          the lines aren't disjoint. The first column will contain the
          parameters at the beginning of the shared segment and the second
          column will correspond to the end of the shared segment.
    """
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-branches
    delta0 = end0 - start0
    line0_const = _py_helpers.cross_product(start0, delta0)
    start1_against = _py_helpers.cross_product(start1, delta0)
    if line0_const != start1_against:
        return True, None

    # Each array is a 1D vector, so we can use the vector dot product.
    norm0_sq = np.vdot(delta0, delta0)
    #                S1 = L1(0) = S0 + sA D0
    # <==>        sA D0 = S1 - S0
    #  ==> sA (D0^T D0) = D0^T (S1 - S0)
    s_val0 = np.vdot(start1 - start0, delta0) / norm0_sq
    #                E1 = L1(1) = S0 + sB D0
    # <==>        sB D0 = E1 - S0
    #  ==> sB (D0^T D0) = D0^T (E1 - S0)
    s_val1 = np.vdot(end1 - start0, delta0) / norm0_sq
    # s = s_val0 + t (s_val1 - s_val0)
    # t = 0                                <==> s = s_val0
    # t = 1                                <==> s = s_val1
    # t = -s_val0 / (s_val1 - s_val0)      <==> s = 0
    # t = (1 - s_val0) / (s_val1 - s_val0) <==> s = 1
    if s_val0 <= s_val1:
        # In this branch the segments are moving in the same direction, i.e.
        # (t=0<-->s=s_val0) are both less than (t=1<-->s_val1).
        if 1.0 < s_val0:
            return True, None

        elif s_val0 < 0.0:
            start_s = 0.0
            start_t = -s_val0 / (s_val1 - s_val0)
        else:
            start_s = s_val0
            start_t = 0.0
        if s_val1 < 0.0:
            return True, None

        elif 1.0 < s_val1:
            end_s = 1.0
            end_t = (1.0 - s_val0) / (s_val1 - s_val0)
        else:
            end_s = s_val1
            end_t = 1.0
    else:
        # In this branch the segments are moving in opposite directions, i.e.
        # in (t=0<-->s=s_val0) and (t=1<-->s_val1) we have 0 < 1
        # but ``s_val0 > s_val1``.
        if s_val0 < 0.0:
            return True, None

        elif 1.0 < s_val0:
            start_s = 1.0
            start_t = (s_val0 - 1.0) / (s_val0 - s_val1)
        else:
            start_s = s_val0
            start_t = 0.0
        if 1.0 < s_val1:
            return True, None

        elif s_val1 < 0.0:
            end_s = 0.0
            end_t = s_val0 / (s_val0 - s_val1)
        else:
            end_s = s_val1
            end_t = 1.0
    parameters = np.asfortranarray([[start_s, end_s], [start_t, end_t]])
    return False, parameters


def line_line_collide(line1, line2):
    """Determine if two line segments meet.

    This is a helper for :func:`convex_hull_collide` in the
    special case that the two convex hulls are actually
    just line segments. (Even in this case, this is only
    problematic if both segments are on a single line.)

    Args:
        line1 (numpy.ndarray): ``2 x 2`` array of start and end nodes.
        line2 (numpy.ndarray): ``2 x 2`` array of start and end nodes.

    Returns:
        bool: Indicating if the line segments collide.
    """
    s, t, success = segment_intersection(
        line1[:, 0], line1[:, 1], line2[:, 0], line2[:, 1]
    )
    if success:
        return _py_helpers.in_interval(
            s, 0.0, 1.0
        ) and _py_helpers.in_interval(t, 0.0, 1.0)

    else:
        disjoint, _ = parallel_lines_parameters(
            line1[:, 0], line1[:, 1], line2[:, 0], line2[:, 1]
        )
        return not disjoint


def convex_hull_collide(nodes1, nodes2):
    """Determine if the convex hulls of two curves collide.

    .. note::

       This is a helper for :func:`from_linearized`.

    Args:
        nodes1 (numpy.ndarray): Control points of a first curve.
        nodes2 (numpy.ndarray): Control points of a second curve.

    Returns:
        bool: Indicating if the convex hulls collide.
    """
    polygon1 = _py_helpers.simple_convex_hull(nodes1)
    _, polygon_size1 = polygon1.shape
    polygon2 = _py_helpers.simple_convex_hull(nodes2)
    _, polygon_size2 = polygon2.shape
    if polygon_size1 == 2 and polygon_size2 == 2:
        return line_line_collide(polygon1, polygon2)

    else:
        return _py_helpers.polygon_collide(polygon1, polygon2)


def from_linearized(first, second, intersections):
    """Determine curve-curve intersection from pair of linearizations.

    .. note::

       This assumes that at least one of ``first`` and ``second`` is
       not a line. The line-line case should be handled "early"
       by :func:`check_lines`.

    .. note::

       This assumes the caller has verified that the bounding boxes
       for ``first`` and ``second`` actually intersect.

    If there is an intersection along the segments, adds that intersection
    to ``intersections``. Otherwise, returns without doing anything.

    Args:
        first (Linearization): First curve being intersected.
        second (Linearization): Second curve being intersected.
        intersections (list): A list of existing intersections.

    Raises:
        ValueError: If ``first`` and ``second`` both have linearization error
            of ``0.0`` (i.e. they are both lines). This is because this
            function expects the caller to have used :func:`check_lines`
            already.
    """
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-return-statements
    s, t, success = segment_intersection(
        first.start_node, first.end_node, second.start_node, second.end_node
    )
    bad_parameters = False
    if success:
        if not (
            _py_helpers.in_interval(s, 0.0, 1.0)
            and _py_helpers.in_interval(t, 0.0, 1.0)
        ):
            bad_parameters = True
    else:
        if first.error == 0.0 and second.error == 0.0:
            raise ValueError(_UNHANDLED_LINES)

        # Just fall back to a Newton iteration starting in the middle of
        # the given intervals.
        bad_parameters = True
        s = 0.5
        t = 0.5
    if bad_parameters:
        # In the unlikely case that we have parallel segments or segments
        # that intersect outside of [0, 1] x [0, 1], we can still exit
        # if the convex hulls don't intersect.
        if not convex_hull_collide(first.curve.nodes, second.curve.nodes):
            return

    # Now, promote ``s`` and ``t`` onto the original curves.
    orig_s = (1 - s) * first.curve.start + s * first.curve.end
    orig_t = (1 - t) * second.curve.start + t * second.curve.end
    refined_s, refined_t = intersection_helpers.full_newton(
        orig_s, first.curve.original_nodes, orig_t, second.curve.original_nodes
    )
    refined_s, success = _py_helpers.wiggle_interval(refined_s)
    if not success:
        return

    refined_t, success = _py_helpers.wiggle_interval(refined_t)
    if not success:
        return

    add_intersection(refined_s, refined_t, intersections)


def add_intersection(s, t, intersections):
    r"""Adds an intersection to list of ``intersections``.

    .. note::

       This is a helper for :func:`from_linearized` and :func:`endpoint_check`.
       These functions are used (directly or indirectly) by
       :func:`all_intersections` exclusively, and that function has a
       Fortran equivalent.

    Accounts for repeated intersection points. If the intersection has already
    been found, does nothing.

    If ``s`` is below :math:`2^{-10}`, it will be replaced with ``1 - s``
    and compared against ``1 - s'`` for all ``s'`` already in
    ``intersections``. (Similar if ``t`` is below the
    :attr:`~bezier.hazmat.intersection_helpers.ZERO_THRESHOLD`.)
    This is perfectly "appropriate" since evaluating a B |eacute| zier curve
    requires using both ``s`` and ``1 - s``, so both values are equally
    relevant.

    Compares :math:`\|p - q\|` to :math:`\|p\|` where :math:`p = (s, t)` is
    current candidate intersection (or the "normalized" version, such as
    :math:`p = (1 - s, t)`) and :math:`q` is one of the already added
    intersections. If the difference is below :math:`2^{-36}` (i.e.
    :attr:`~bezier.hazmat.intersection_helpers.NEWTON_ERROR_RATIO`)
    then the intersection is considered to be duplicate.

    Args:
        s (float): The first parameter in an intersection.
        t (float): The second parameter in an intersection.
        intersections (list): List of existing intersections.
    """
    if not intersections:
        intersections.append((s, t))
        return

    if s < intersection_helpers.ZERO_THRESHOLD:
        candidate_s = 1.0 - s
    else:
        candidate_s = s
    if t < intersection_helpers.ZERO_THRESHOLD:
        candidate_t = 1.0 - t
    else:
        candidate_t = t
    norm_candidate = np.linalg.norm([candidate_s, candidate_t], ord=2)
    for existing_s, existing_t in intersections:
        # NOTE: |(1 - s1) - (1 - s2)| = |s1 - s2| in exact arithmetic, so
        #       we just compute ``s1 - s2`` rather than using
        #       ``candidate_s`` / ``candidate_t``. Due to round-off, these
        #       differences may be slightly different, but only up to machine
        #       precision.
        delta_s = s - existing_s
        delta_t = t - existing_t
        norm_update = np.linalg.norm([delta_s, delta_t], ord=2)
        if (
            norm_update
            < intersection_helpers.NEWTON_ERROR_RATIO * norm_candidate
        ):
            return

    intersections.append((s, t))


def endpoint_check(
    first, node_first, s, second, node_second, t, intersections
):
    r"""Check if curve endpoints are identical.

    .. note::

       This is a helper for :func:`tangent_bbox_intersection`. These
       functions are used (directly or indirectly) by
       :func:`all_intersections` exclusively, and that function has a
       Fortran equivalent.

    Args:
        first (SubdividedCurve): First curve being intersected (assumed in
            :math:`\mathbf{R}^2`).
        node_first (numpy.ndarray): 1D ``2``-array, one of the endpoints
            of ``first``.
        s (float): The parameter corresponding to ``node_first``, so
             expected to be one of ``0.0`` or ``1.0``.
        second (SubdividedCurve): Second curve being intersected (assumed in
            :math:`\mathbf{R}^2`).
        node_second (numpy.ndarray): 1D ``2``-array, one of the endpoints
            of ``second``.
        t (float): The parameter corresponding to ``node_second``, so
             expected to be one of ``0.0`` or ``1.0``.
        intersections (list): A list of already encountered
            intersections. If these curves intersect at their tangency,
            then those intersections will be added to this list.
    """
    if _py_helpers.vector_close(node_first, node_second):
        orig_s = (1 - s) * first.start + s * first.end
        orig_t = (1 - t) * second.start + t * second.end
        add_intersection(orig_s, orig_t, intersections)


def tangent_bbox_intersection(first, second, intersections):
    r"""Check if two curves with tangent bounding boxes intersect.

    .. note::

       This is a helper for :func:`intersect_one_round`. These
       functions are used (directly or indirectly) by
       :func:`all_intersections` exclusively, and that function has a
       Fortran equivalent.

    If the bounding boxes are tangent, intersection can
    only occur along that tangency.

    If the curve is **not** a line, the **only** way the curve can touch
    the bounding box is at the endpoints. To see this, consider the
    component

    .. math::

       x(s) = \sum_j W_j x_j.

    Since :math:`W_j > 0` for :math:`s \in \left(0, 1\right)`, if there
    is some :math:`k` with :math:`x_k < M = \max x_j`, then for any
    interior :math:`s`

    .. math::

       x(s) < \sum_j W_j M = M.

    If all :math:`x_j = M`, then :math:`B(s)` falls on the line
    :math:`x = M`. (A similar argument holds for the other three
    component-extrema types.)

    .. note::

       This function assumes callers will not pass curves that can be
       linearized / are linear. In :func:`all_intersections`, curves
       are pre-processed to do any linearization before the
       subdivision / intersection process begins.

    Args:
        first (SubdividedCurve): First curve being intersected (assumed in
            :math:`\mathbf{R}^2`).
        second (SubdividedCurve): Second curve being intersected (assumed in
            :math:`\mathbf{R}^2`).
        intersections (list): A list of already encountered
            intersections. If these curves intersect at their tangency,
            then those intersections will be added to this list.
    """
    node_first1 = first.nodes[:, 0]
    node_first2 = first.nodes[:, -1]
    node_second1 = second.nodes[:, 0]
    node_second2 = second.nodes[:, -1]
    endpoint_check(
        first, node_first1, 0.0, second, node_second1, 0.0, intersections
    )
    endpoint_check(
        first, node_first1, 0.0, second, node_second2, 1.0, intersections
    )
    endpoint_check(
        first, node_first2, 1.0, second, node_second1, 0.0, intersections
    )
    endpoint_check(
        first, node_first2, 1.0, second, node_second2, 1.0, intersections
    )


def bbox_line_intersect(nodes, line_start, line_end):
    r"""Determine intersection of a bounding box and a line.

    We do this by first checking if either the start or end node of the
    segment are contained in the bounding box. If they aren't, then
    checks if the line segment intersects any of the four sides of the
    bounding box.

    .. note::

       This function is "half-finished". It makes no distinction between
       "tangent" intersections of the box and segment and other types
       of intersection. However, the distinction is worthwhile, so this
       function should be "upgraded" at some point.

    Args:
        nodes (numpy.ndarray): Points (``2 x N``) that determine a
            bounding box.
        line_start (numpy.ndarray): Beginning of a line segment (1D
            ``2``-array).
        line_end (numpy.ndarray): End of a line segment (1D ``2``-array).

    Returns:
        int: Enum from :class:`.BoxIntersectionType` indicating the type of
        bounding box intersection.
    """
    left, right, bottom, top = _py_helpers.bbox(nodes)
    if _py_helpers.in_interval(
        line_start[0], left, right
    ) and _py_helpers.in_interval(line_start[1], bottom, top):
        return BoxIntersectionType.INTERSECTION

    if _py_helpers.in_interval(
        line_end[0], left, right
    ) and _py_helpers.in_interval(line_end[1], bottom, top):
        return BoxIntersectionType.INTERSECTION

    # NOTE: We allow ``segment_intersection`` to fail below (i.e.
    #       ``success=False``). At first, this may appear to "ignore"
    #       some potential intersections of parallel lines. However,
    #       no intersections will be missed. If parallel lines don't
    #       overlap, then there is nothing to miss. If they do overlap,
    #       then either the segment will have endpoints on the box (already
    #       covered by the checks above) or the segment will contain an
    #       entire side of the box, which will force it to intersect the 3
    #       edges that meet at the two ends of those sides. The parallel
    #       edge will be skipped, but the other two will be covered.
    # Bottom Edge
    s_bottom, t_bottom, success = segment_intersection(
        np.asfortranarray([left, bottom]),
        np.asfortranarray([right, bottom]),
        line_start,
        line_end,
    )
    if (
        success
        and _py_helpers.in_interval(s_bottom, 0.0, 1.0)
        and _py_helpers.in_interval(t_bottom, 0.0, 1.0)
    ):
        return BoxIntersectionType.INTERSECTION

    # Right Edge
    s_right, t_right, success = segment_intersection(
        np.asfortranarray([right, bottom]),
        np.asfortranarray([right, top]),
        line_start,
        line_end,
    )
    if (
        success
        and _py_helpers.in_interval(s_right, 0.0, 1.0)
        and _py_helpers.in_interval(t_right, 0.0, 1.0)
    ):
        return BoxIntersectionType.INTERSECTION

    # Top Edge
    s_top, t_top, success = segment_intersection(
        np.asfortranarray([right, top]),
        np.asfortranarray([left, top]),
        line_start,
        line_end,
    )
    if (
        success
        and _py_helpers.in_interval(s_top, 0.0, 1.0)
        and _py_helpers.in_interval(t_top, 0.0, 1.0)
    ):
        return BoxIntersectionType.INTERSECTION

    # NOTE: We skip the "last" edge. This is because any curve
    #       that doesn't have an endpoint on a curve must cross
    #       at least two, so we will have already covered such curves
    #       in one of the branches above.
    return BoxIntersectionType.DISJOINT


def intersect_one_round(candidates, intersections):
    """Perform one step of the intersection process.

    .. note::

       This is a helper for :func:`all_intersections` and that function
       has a Fortran equivalent.

    Checks if the bounding boxes of each pair in ``candidates``
    intersect. If the bounding boxes do not intersect, the pair
    is discarded. Otherwise, the pair is "accepted". Then we
    attempt to linearize each curve in an "accepted" pair and
    track the overall linearization error for every curve
    encountered.

    Args:
        candidates (Union[list, itertools.chain]): An iterable of
            pairs of curves (or linearized curves).
        intersections (list): A list of already encountered
            intersections. If any intersections can be readily determined
            during this round of subdivision, then they will be added
            to this list.

    Returns:
        list: Returns a list of the next round of ``candidates``.
    """
    next_candidates = []
    # NOTE: In the below we replace ``isinstance(a, B)`` with
    #       ``a.__class__ is B``, which is a 3-3.5x speedup.
    for first, second in candidates:
        both_linearized = False
        if first.__class__ is Linearization:
            if second.__class__ is Linearization:
                both_linearized = True
                bbox_int = bbox_intersect(
                    first.curve.nodes, second.curve.nodes
                )
            else:
                bbox_int = bbox_line_intersect(
                    second.nodes, first.start_node, first.end_node
                )
        else:
            if second.__class__ is Linearization:
                bbox_int = bbox_line_intersect(
                    first.nodes, second.start_node, second.end_node
                )
            else:
                bbox_int = bbox_intersect(first.nodes, second.nodes)
        if bbox_int == BoxIntersectionType.DISJOINT:
            continue

        if bbox_int == BoxIntersectionType.TANGENT and not both_linearized:
            # NOTE: Ignore tangent bounding boxes in the linearized case
            #       because ``tangent_bbox_intersection()`` assumes that both
            #       curves are not linear.
            tangent_bbox_intersection(first, second, intersections)
            continue

        if both_linearized:
            # If both ``first`` and ``second`` are linearizations, then
            # we can intersect them immediately.
            from_linearized(first, second, intersections)
            continue

        # If we haven't ``continue``-d, add the accepted pair.
        # NOTE: This may be a wasted computation, e.g. if ``first``
        #       or ``second`` occur in multiple accepted pairs (the caller
        #       only passes one pair at a time). However, in practice
        #       the number of such pairs will be small so this cost
        #       will be low.
        lin1 = map(Linearization.from_shape, first.subdivide())
        lin2 = map(Linearization.from_shape, second.subdivide())
        next_candidates.extend(itertools.product(lin1, lin2))
    return next_candidates


def prune_candidates(candidates):
    """Reduce number of candidate intersection pairs.

    .. note::

       This is a helper for :func:`all_intersections`.

    Uses more strict bounding box intersection predicate by forming the
    actual convex hull of each candidate curve segment and then checking
    if those convex hulls collide.

    Args:
        candidates (List[Union[SubdividedCurve, Linearization]]): An iterable
            of pairs of curves (or linearized curves).

    Returns:
        List[Union[SubdividedCurve, Linearization]]: A pruned list of curve
        pairs.
    """
    pruned = []
    # NOTE: In the below we replace ``isinstance(a, B)`` with
    #       ``a.__class__ is B``, which is a 3-3.5x speedup.
    for first, second in candidates:
        if first.__class__ is Linearization:
            nodes1 = first.curve.nodes
        else:
            nodes1 = first.nodes
        if second.__class__ is Linearization:
            nodes2 = second.curve.nodes
        else:
            nodes2 = second.nodes
        if convex_hull_collide(nodes1, nodes2):
            pruned.append((first, second))
    return pruned


def make_same_degree(nodes1, nodes2):
    """Degree-elevate a curve so two curves have matching degree.

    Args:
        nodes1 (numpy.ndarray): Set of control points for a
            B |eacute| zier curve.
        nodes2 (numpy.ndarray): Set of control points for a
            B |eacute| zier curve.

    Returns:
        Tuple[numpy.ndarray, numpy.ndarray]: The potentially degree-elevated
        nodes passed in.
    """
    _, num_nodes1 = nodes1.shape
    _, num_nodes2 = nodes2.shape
    for _ in range(num_nodes2 - num_nodes1):
        nodes1 = curve_helpers.elevate_nodes(nodes1)
    for _ in range(num_nodes1 - num_nodes2):
        nodes2 = curve_helpers.elevate_nodes(nodes2)
    return nodes1, nodes2


def coincident_parameters(nodes1, nodes2):
    r"""Check if two B |eacute| zier curves are coincident.

    Does so by projecting each segment endpoint onto the other curve

    .. math::

       B_1(s_0) = B_2(0) \\
       B_1(s_m) = B_2(1) \\
       B_1(0) = B_2(t_0) \\
       B_1(1) = B_2(t_n)

    and then finding the "shared interval" where both curves are defined.
    If such an interval can't be found (e.g. if one of the endpoints can't be
    located on the other curve), returns :data:`None`.

    If such a "shared interval" does exist, then this will specialize
    each curve onto that shared interval and check if the new control points
    agree.

    Args:
        nodes1 (numpy.ndarray): Set of control points for a
            B |eacute| zier curve.
        nodes2 (numpy.ndarray): Set of control points for a
            B |eacute| zier curve.

    Returns:
        Optional[Tuple[Tuple[float, float], ...]]: A ``2 x 2`` array of
        parameters where the two coincident curves meet. If they are not
        coincident, returns :data:`None`.
    """
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-return-statements,too-many-branches
    nodes1, nodes2 = make_same_degree(nodes1, nodes2)
    s_initial = curve_helpers.locate_point(
        nodes1, nodes2[:, 0].reshape((2, 1), order="F")
    )
    s_final = curve_helpers.locate_point(
        nodes1, nodes2[:, -1].reshape((2, 1), order="F")
    )
    if s_initial is not None and s_final is not None:
        # In this case, if the curves were coincident, then ``curve2``
        # would be "fully" contained in ``curve1``, so we specialize
        # ``curve1`` down to that interval to check.
        specialized1 = curve_helpers.specialize_curve(
            nodes1, s_initial, s_final
        )
        if _py_helpers.vector_close(
            specialized1.ravel(order="F"), nodes2.ravel(order="F")
        ):
            return ((s_initial, 0.0), (s_final, 1.0))

        else:
            return None

    t_initial = curve_helpers.locate_point(
        nodes2, nodes1[:, 0].reshape((2, 1), order="F")
    )
    t_final = curve_helpers.locate_point(
        nodes2, nodes1[:, -1].reshape((2, 1), order="F")
    )
    if t_initial is None and t_final is None:
        # An overlap must have two endpoints and since at most one of the
        # endpoints of ``curve2`` lies on ``curve1`` (as indicated by at
        # least one of the ``s``-parameters being ``None``), we need (at least)
        # one endpoint of ``curve1`` on ``curve2``.
        return None

    if t_initial is not None and t_final is not None:
        # In this case, if the curves were coincident, then ``curve1``
        # would be "fully" contained in ``curve2``, so we specialize
        # ``curve2`` down to that interval to check.
        specialized2 = curve_helpers.specialize_curve(
            nodes2, t_initial, t_final
        )
        if _py_helpers.vector_close(
            nodes1.ravel(order="F"), specialized2.ravel(order="F")
        ):
            return ((0.0, t_initial), (1.0, t_final))

        else:
            return None

    if s_initial is None and s_final is None:
        # An overlap must have two endpoints and since exactly one of the
        # endpoints of ``curve1`` lies on ``curve2`` (as indicated by exactly
        # one of the ``t``-parameters being ``None``), we need (at least)
        # one endpoint of ``curve1`` on ``curve2``.
        return None

    # At this point, we know exactly one of the ``s``-parameters and exactly
    # one of the ``t``-parameters is not ``None``.
    if s_initial is None:
        if t_initial is None:
            # B1(s_final) = B2(1) AND B1(1) = B2(t_final)
            start_s = s_final
            end_s = 1.0
            start_t = 1.0
            end_t = t_final
        else:
            # B1(0) = B2(t_initial) AND B1(s_final) = B2(1)
            start_s = 0.0
            end_s = s_final
            start_t = t_initial
            end_t = 1.0
    else:
        if t_initial is None:
            # B1(s_initial) = B2(0) AND B1(1 ) = B2(t_final)
            start_s = s_initial
            end_s = 1.0
            start_t = 0.0
            end_t = t_final
        else:
            # B1(0) = B2(t_initial) AND B1(s_initial) = B2(0)
            start_s = 0.0
            end_s = s_initial
            start_t = t_initial
            end_t = 0.0
    width_s = abs(start_s - end_s)
    width_t = abs(start_t - end_t)
    if width_s < _MIN_INTERVAL_WIDTH and width_t < _MIN_INTERVAL_WIDTH:
        return None

    specialized1 = curve_helpers.specialize_curve(nodes1, start_s, end_s)
    specialized2 = curve_helpers.specialize_curve(nodes2, start_t, end_t)
    if _py_helpers.vector_close(
        specialized1.ravel(order="F"), specialized2.ravel(order="F")
    ):
        return ((start_s, start_t), (end_s, end_t))

    else:
        return None


def check_lines(first, second):
    """Checks if two curves are lines and tries to intersect them.

    .. note::

       This is a helper for :func:`.all_intersections`.

    If they are not lines / not linearized, immediately returns :data:`False`
    with no "return value".

    If they are lines, attempts to intersect them (even if they are parallel
    and share a coincident segment).

    Args:
        first (Union[SubdividedCurve, Linearization]): First curve being
            intersected.
        second (Union[SubdividedCurve, Linearization]): Second curve being
            intersected.

    Returns:
        Tuple[bool, Optional[Tuple[numpy.ndarray, bool]]]: A pair of

        * Flag indicating if both candidates in the pair are lines.
        * Optional "result" populated only if both candidates are lines.
          When this result is populated, it will be a pair of

          * array of parameters of intersection
          * flag indicating if the two candidates share a coincident segment
    """
    # NOTE: In the below we replace ``isinstance(a, B)`` with
    #       ``a.__class__ is B``, which is a 3-3.5x speedup.
    if not (
        first.__class__ is Linearization
        and second.__class__ is Linearization
        and first.error == 0.0
        and second.error == 0.0
    ):
        return False, None

    s, t, success = segment_intersection(
        first.start_node, first.end_node, second.start_node, second.end_node
    )
    if success:
        if _py_helpers.in_interval(s, 0.0, 1.0) and _py_helpers.in_interval(
            t, 0.0, 1.0
        ):
            intersections = np.asfortranarray([[s], [t]])
            result = intersections, False
        else:
            result = np.empty((2, 0), order="F"), False
    else:
        disjoint, params = parallel_lines_parameters(
            first.start_node,
            first.end_node,
            second.start_node,
            second.end_node,
        )
        if disjoint:
            result = np.empty((2, 0), order="F"), False
        else:
            result = params, True
    return True, result


def all_intersections(nodes_first, nodes_second):
    r"""Find the points of intersection among a pair of curves.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    .. note::

       This assumes both curves are in :math:`\mathbf{R}^2`, but does not
       **explicitly** check this. However, functions used here will fail if
       that assumption fails, e.g. :func:`bbox_intersect` and
       :func:`newton_refine() <.hazmat.intersection_helpers.newton_refine>`.

    Args:
        nodes_first (numpy.ndarray): Control points of a curve to be
            intersected with ``nodes_second``.
        nodes_second (numpy.ndarray): Control points of a curve to be
            intersected with ``nodes_first``.

    Returns:
        Tuple[numpy.ndarray, bool]: An array and a flag:

        * A ``2 x N`` array of intersection parameters.
          Each row contains a pair of values :math:`s` and :math:`t`
          (each in :math:`\left[0, 1\right]`) such that the curves
          intersect: :math:`B_1(s) = B_2(t)`.
        * Flag indicating if the curves are coincident.

    Raises:
        ValueError: If the subdivision iteration does not terminate
            before exhausting the maximum number of subdivisions.
        NotImplementedError: If the subdivision process picks up too
            many candidate pairs. This typically indicates tangent
            curves or coincident curves (though there are mitigations for
            those cases in place).
    """
    curve_first = SubdividedCurve(nodes_first, nodes_first)
    curve_second = SubdividedCurve(nodes_second, nodes_second)
    candidate1 = Linearization.from_shape(curve_first)
    candidate2 = Linearization.from_shape(curve_second)
    # Handle the line-line intersection case as a one-off.
    both_linear, result = check_lines(candidate1, candidate2)
    if both_linear:
        return result

    candidates = [(candidate1, candidate2)]
    intersections = []
    coincident = False
    for _ in range(_MAX_INTERSECT_SUBDIVISIONS):
        candidates = intersect_one_round(candidates, intersections)
        if len(candidates) > _MAX_CANDIDATES:
            candidates = prune_candidates(candidates)
            # If pruning didn't fix anything, we check if the curves are
            # coincident and "fail" if they aren't.
            if len(candidates) > _MAX_CANDIDATES:
                params = coincident_parameters(nodes_first, nodes_second)
                if params is None:
                    raise NotImplementedError(
                        _TOO_MANY_TEMPLATE.format(len(candidates))
                    )

                intersections = params
                coincident = True
                # Artificially empty out candidates so that this
                # function exits.
                candidates = []
        # If none of the candidate pairs have been accepted, then there are
        # no more intersections to find.
        if not candidates:
            if intersections:
                # NOTE: The transpose of a C-ordered array is Fortran-ordered,
                #       i.e. this is on purpose.
                return np.array(intersections, order="C").T, coincident

            return np.empty((2, 0), order="F"), coincident

    msg = _NO_CONVERGE_TEMPLATE.format(_MAX_INTERSECT_SUBDIVISIONS)
    raise ValueError(msg)


class BoxIntersectionType:  # pylint: disable=too-few-public-methods
    """Enum representing all possible bounding box intersections.

    .. note::

       This class would be more "correct" as an :class:`enum.Enum`, but it we
       keep the values integers to make interfacing with Fortran easier.
    """

    INTERSECTION = 0
    """Bounding boxes overlap with positive area."""
    TANGENT = 1
    """Bounding boxes are tangent but do not overlap."""
    DISJOINT = 2
    """Bounding boxes do not intersect or touch."""


class SubdividedCurve:  # pylint: disable=too-few-public-methods
    """A data wrapper for a B |eacute| zier curve

    To be used for intersection algorithm via repeated subdivision,
    where the ``start`` and ``end`` parameters must be tracked.

    Args:
        nodes (numpy.ndarray): The control points of the current
            subdivided curve
        original_nodes (numpy.ndarray): The control points of the original
            curve used to define the current one (before subdivision began).
        start (Optional[float]): The start parameter after subdivision.
        end (Optional[float]): The end parameter after subdivision.
    """

    __slots__ = ("nodes", "original_nodes", "start", "end")

    def __init__(self, nodes, original_nodes, start=0.0, end=1.0):
        self.nodes = nodes
        self.original_nodes = original_nodes
        self.start = start
        self.end = end

    @property
    def __dict__(self):
        """dict: Dictionary of current subdivided curve's property namespace.

        This is just a stand-in property for the usual ``__dict__``. This
        class defines ``__slots__`` so by default would not provide a
        ``__dict__``.

        This also means that the current object can't be modified by the
        returned dictionary.
        """
        return {
            "nodes": self.nodes,
            "original_nodes": self.original_nodes,
            "start": self.start,
            "end": self.end,
        }

    def subdivide(self):
        """Split the curve into a left and right half.

        See :meth:`.Curve.subdivide` for more information.

        Returns:
            Tuple[SubdividedCurve, SubdividedCurve]: The left and right
            sub-curves.
        """
        left_nodes, right_nodes = curve_helpers.subdivide_nodes(self.nodes)
        midpoint = 0.5 * (self.start + self.end)
        left = SubdividedCurve(
            left_nodes, self.original_nodes, start=self.start, end=midpoint
        )
        right = SubdividedCurve(
            right_nodes, self.original_nodes, start=midpoint, end=self.end
        )
        return left, right


class Linearization:
    """A linearization of a curve.

    This class is provided as a stand-in for a curve, so it
    provides a similar interface.

    Args:
        curve (SubdividedCurve): A curve that is linearized.
        error (float): The linearization error. Expected to have been
            computed via :func:`linearization_error`.
    """

    __slots__ = ("curve", "error", "start_node", "end_node")

    def __init__(self, curve, error):
        self.curve = curve
        """SubdividedCurve: The curve that this linearization approximates."""
        self.error = error
        """float: The linearization error for the linearized curve."""
        self.start_node = curve.nodes[:, 0]
        """numpy.ndarray: The 1D start vector of this linearization."""
        self.end_node = curve.nodes[:, -1]
        """numpy.ndarray: The 1D end vector of this linearization."""

    @property
    def __dict__(self):
        """dict: Dictionary of current linearization's property namespace.

        This is just a stand-in property for the usual ``__dict__``. This
        class defines ``__slots__`` so by default would not provide a
        ``__dict__``.

        This also means that the current object can't be modified by the
        returned dictionary.
        """
        return {
            "curve": self.curve,
            "error": self.error,
            "start_node": self.start_node,
            "end_node": self.end_node,
        }

    def subdivide(self):
        """Do-nothing method to match the :class:`.Curve` interface.

        Returns:
            Tuple[~bezier.hazmat.geometric_intersection.Linearization]: List of
            all subdivided parts, which is just the current object.
        """
        return (self,)

    @classmethod
    def from_shape(cls, shape):
        """Try to linearize a curve (or an already linearized curve).

        Args:
            shape (Union[SubdividedCurve, \
            ~bezier.hazmat.geometric_intersection.Linearization]): A curve or
                an already linearized curve.

        Returns:
            Union[SubdividedCurve, \
            ~bezier.hazmat.geometric_intersection.Linearization]: The
            (potentially linearized) curve.
        """
        # NOTE: In the below we replace ``isinstance(a, B)`` with
        #       ``a.__class__ is B``, which is a 3-3.5x speedup.
        if shape.__class__ is cls:
            return shape

        else:
            error = linearization_error(shape.nodes)
            if error < _ERROR_VAL:
                linearized = cls(shape, error)
                return linearized

            else:
                return shape


def self_intersections(nodes):
    r"""Determine the self-intersections of a planar B |eacute| zier curve.

    See: https://doi.org/10.1016/0166-3615(89)90072-9

    Args:
        nodes (numpy.ndarray): The nodes in the curve.

    Returns:
        numpy.ndarray: Pairs of parameters along the curve where
        self-intersections occur (as a 2D array).
    """
    angle = curve_helpers.discrete_turning_angle(nodes)
    if angle < np.pi:
        return np.empty((2, 0), order="F")

    left_nodes, right_nodes = curve_helpers.subdivide_nodes(nodes)
    #  left_nodes: [0.0, 0.5] <- [0.0, 1.0]
    left_self = 0.5 * self_intersections(left_nodes)
    # right_nodes: [0.5, 1.0] <- [0.0, 1.0]
    right_self = 0.5 + 0.5 * self_intersections(right_nodes)

    left_right_intersections, _ = all_intersections(left_nodes, right_nodes)
    left_right_intersections[0, :] *= 0.5
    left_right_intersections[1, :] *= 0.5
    left_right_intersections[1, :] += 0.5

    # ``left_right_intersections`` should contain a redundant intersection
    # for the split point: left(1) = right(0).
    compare_split = left_right_intersections == _SPLIT_POINT
    keep_columns = ~np.all(compare_split, axis=0)
    left_right_intersections = left_right_intersections[:, keep_columns]

    result = np.hstack([left_self, left_right_intersections, right_self])
    return np.asfortranarray(result)
