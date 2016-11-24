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


import enum
import itertools

import numpy as np
import six

from bezier import _curve_helpers


_FREXP = np.frexp  # pylint: disable=no-member
# Set the threshold for exponetn at half the bits available,
# this way one round of Newton's method can finish the job
# by squaring the error.
_ERROR_EXPONENT = -26
_MAX_INTERSECT_SUBDIVISIONS = 20
_EPS = 2.0**(-40)
_TOO_MANY_TEMPLATE = (
    'The number of candidate intersections is too high.\n'
    '{:d} accepted pairs gives {:d} candidate pairs.')
# Allow wiggle room for ``s`` and ``t`` computed during segment
# intersection. Any over- or under-shooting will (hopefully) be
# resolved in the Newton refinement step. If it isn't resolved, the
# call to _check_parameters() will fail the intersection.
_WIGGLE_START = -2.0**(-16)
_WIGGLE_END = 1.0 - _WIGGLE_START


def _in_interval(value, start, end):
    """Checks if a ``value`` is an interval (inclusive).

    .. note::

       The current implementation does the most basic check,
       however, in the future, a more generic check may be desired
       that allows wiggle room around the endpoints to account
       for round-off.

    Args:
        value (float): The value to check.
        start (float): The (inclusive) start of the interval.
        end (float): The (inclusive) end of the interval.

    Returns:
        bool: Indicating if the value is in the interval.
    """
    return start <= value <= end


def _vector_close(vec1, vec2):
    r"""Checks that two vectors are equal to some threshold.

    Does so by computing :math:`s_1 = \|v_1\|_2` and
    :math:`s_2 = \|v_2\|_2` and then checking if

    .. math::

       \|v_1 - v_2\|_2 \leq \varepsilon \min(s_1, s_2)

    where :math:`\varepsilon = 2^{-40} \approx 10^{-12}` is a fixed
    threshold. In the rare case that one of ``vec1`` or ``vec2`` is
    the zero vector (i.e. when :math:`\min(s_1, s_2) = 0`) instead
    checks that the other vector is close enough to zero:

    .. math::

       \|v_1\|_2 = 0 \Longrightarrow \|v_2\|_2 \leq \varepsilon

    .. note::

       This function assumes that both vectors have finite numbers,
       i.e. that no NaN or infinite numbers occur. NumPy provides
       :func:`np.allclose` for coverage of **all** cases.

    Args:
        vec1 (numpy.ndarray): First vector for comparison.
        vec2 (numpy.ndarray): Second vector for comparison.

    Returns:
        bool: Flag indicating if they are close to precision.
    """
    size1 = np.linalg.norm(vec1, ord=2)
    size2 = np.linalg.norm(vec2, ord=2)
    if size1 == 0:
        return size2 <= _EPS
    elif size2 == 0:
        return size1 <= _EPS
    else:
        upper_bound = _EPS * min(size1, size2)
        return np.linalg.norm(vec1 - vec2, ord=2) <= upper_bound


def _check_close(s, curve1, t, curve2):
    r"""Checks that two curves intersect to some threshold.

    Verifies :math:`B_1(s) \approx B_2(t)` and then returns
    :math:`B_1(s)` if they are sufficiently close.

    Args:
        s (float): Parameter of a near-intersection along ``curve1``.
        curve1 (.Curve): First curve forming intersection.
        t (float): Parameter of a near-intersection along ``curve2``.
        curve2 (.Curve): Second curve forming intersection.

    Raises:
        ValueError: If :math:`B_1(s)` is not sufficiently close to
            :math:`B_2(t)`.

    Returns:
        numpy.ndarray: The value of :math:`B_1(s)`.
    """
    vec1 = curve1.evaluate(s)
    vec2 = curve2.evaluate(t)
    if not _vector_close(vec1, vec2):
        raise ValueError('B_1(s) and B_2(t) are not sufficiently close')

    return vec1


def _check_parameters(s, t):
    r"""Check if ``s``, ``t`` are in :math:`\left(0, 1\right)`.

    Args:
        s (float): Parameter on a curve.
        t (float): Parameter on a curve.

    Raises:
        ValueError: If one of the values falls outside the unit interval.
    """
    if not _in_interval(s, 0.0, 1.0):
        raise ValueError('s outside of unit interval', s)
    if not _in_interval(t, 0.0, 1.0):
        raise ValueError('t outside of unit interval', t)


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
        BoxIntersectionType: Enum indicating the type of bounding
        box intersection.
    """
    left1, bottom1 = np.min(nodes1, axis=0)
    right1, top1 = np.max(nodes1, axis=0)
    left2, bottom2 = np.min(nodes2, axis=0)
    right2, top2 = np.max(nodes2, axis=0)

    if (right2 < left1 or right1 < left2 or
            top2 < bottom1 or top1 < bottom2):
        return BoxIntersectionType.disjoint

    if (right2 == left1 or right1 == left2 or
            top2 == bottom1 or top1 == bottom2):
        return BoxIntersectionType.tangent
    else:
        return BoxIntersectionType.intersection


def linearization_error(curve):
    r"""Compute the maximum error of a linear approximation.

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
           + \left[\begin{array}{c} 1 \\ 3 \end{array}\right] 2s(1 - s)
           + \left[\begin{array}{c} -2 \\ 9 \end{array}\right] s^2

    has
    :math:`B''(s) \equiv \left[\begin{array}{c} -8 \\ 6 \end{array}\right]`
    which has norm :math:`10` everywhere, hence the maximum error is

    .. math::

       \left.\frac{s(1 - s)}{2!} \cdot 10\right|_{s = \frac{1}{2}}
          = \frac{5}{4}.

    .. testsetup::

       import numpy as np
       import bezier
       from bezier._intersection_helpers import linearization_error

    .. doctest::

       >>> nodes = np.array([
       ...     [ 0.0, 0.0],
       ...     [ 1.0, 3.0],
       ...     [-2.0, 9.0],
       ... ])
       >>> curve = bezier.Curve(nodes)
       >>> linearization_error(curve)
       1.25

    Args:
        curve (~bezier.curve.Curve): A curve to be approximated by a line.

    Returns:
        float: The maximum error between the curve and the
        linear approximation.
    """
    if curve.degree == 1:
        return 0.0

    nodes = curve._nodes  # pylint: disable=protected-access
    second_deriv = nodes[:-2, :] - 2.0 * nodes[1:-1, :] + nodes[2:, :]
    worst_case = np.max(np.abs(second_deriv), axis=0)

    # max_{0 <= s <= 1} s(1 - s)/2 = 1/8 = 0.125
    multiplier = 0.125 * curve.degree * (curve.degree - 1)
    return multiplier * np.linalg.norm(worst_case, ord=2)


def _evaluate_hodograph(nodes, degree, s):
    r"""Evaluate the Hodograph curve at a point :math:`s`.

    The Hodograph (first derivative) of a B |eacute| zier curve
    degree :math:`d = n - 1` and is given by

    .. math::

       B'(s) = n \sum_{j = 0}^{d} \binom{d}{j} s^j
       (1 - s)^{d - j} \cdot \Delta v_j

    where each forward difference is given by
    :math:`\Delta v_j = v_{j + 1} - v_j`.

    Args:
        nodes (numpy.ndarray): The nodes of a curve.
        degree (int): The degree of the curve (assumed to be one less than
            the number of ``nodes``.
        s (float): A parameter along the curve at which the Hodograph
            is to be evaluated.

    Returns:
        numpy.ndarray: The point on the Hodograph curve (as a one
        dimensional NumPy array).
    """
    first_deriv = nodes[1:, :] - nodes[:-1, :]
    # NOTE: Taking the derivative drops the degree by 1.
    return degree * _curve_helpers.evaluate_multi(
        first_deriv, degree - 1, np.array([s])).flatten()


def newton_refine(s, curve1, t, curve2):
    r"""Apply one step of 2D Newton's method.

    We want to use Newton's method on the function

    .. math::

       F(s, t) = B_1(s) - B_2(t)

    to refine :math:`\left(s_{\ast}, t_{\ast}\right)`. Using this,
    and the Jacobian :math:`DF`, we "solve"

    .. math::

       \left[\begin{array}{c}
           0 \\ 0 \end{array}\right] \approx
           F\left(s_{\ast} + \Delta s, t_{\ast} + \Delta t\right) \approx
           F\left(s_{\ast}, t_{\ast}\right) +
           \left[\begin{array}{c c}
               B_1'\left(s_{\ast}\right) &
               - B_2'\left(t_{\ast}\right) \end{array}\right]
           \left[\begin{array}{c}
               \Delta s \\ \Delta t \end{array}\right]

    and refine with the component updates :math:`\Delta s` and
    :math:`\Delta t`.

    .. note::

       This implementation assumes ``curve1`` and ``curve2`` live in
       :math:`\mathbf{R}^2`.

    For example, the curves

    .. math::

        \begin{align*}
        B_1(s) &= \left[\begin{array}{c} 0 \\ 0 \end{array}\right] (1 - s)^2
            + \left[\begin{array}{c} 2 \\ 4 \end{array}\right] 2s(1 - s)
            + \left[\begin{array}{c} 4 \\ 0 \end{array}\right] s^2 \\
        B_2(t) &= \left[\begin{array}{c} 2 \\ 0 \end{array}\right] (1 - t)
            + \left[\begin{array}{c} 0 \\ 3 \end{array}\right] t
        \end{align*}

    intersect at the point
    :math:`B_1\left(\frac{1}{4}\right) = B_2\left(\frac{1}{2}\right) =
    \frac{1}{2} \left[\begin{array}{c} 2 \\ 3 \end{array}\right]`.
    However, starting from the wrong point we have

    .. math::

        \begin{align*}
        F\left(\frac{3}{8}, \frac{1}{4}\right) &= \frac{1}{8}
            \left[\begin{array}{c} 0 \\ 9 \end{array}\right] \\
        DF\left(\frac{3}{8}, \frac{1}{4}\right) &=
            \left[\begin{array}{c c}
            4 & 2 \\ 2 & -3 \end{array}\right] \\
        \Longrightarrow \left[\begin{array}{c} \Delta s \\ \Delta t
            \end{array}\right] &= \frac{9}{64} \left[\begin{array}{c}
            -1 \\ 2 \end{array}\right].
        \end{align*}

    .. testsetup::

       import numpy as np
       import bezier
       from bezier._intersection_helpers import newton_refine

    .. doctest::

       >>> nodes1 = np.array([
       ...     [0.0, 0.0],
       ...     [2.0, 4.0],
       ...     [4.0, 0.0],
       ... ])
       >>> curve1 = bezier.Curve(nodes1)
       >>> nodes2 = np.array([
       ...     [2.0, 0.0],
       ...     [0.0, 3.0],
       ... ])
       >>> curve2 = bezier.Curve(nodes2)
       >>> s, t = 0.375, 0.25
       >>> new_s, new_t = newton_refine(s, curve1, t, curve2)
       >>> 64.0 * (new_s - s)
       -9.0
       >>> 64.0 * (new_t - t)
       18.0

    Args:
        s (float): Parameter of a near-intersection along ``curve1``.
        curve1 (.Curve): First curve forming intersection.
        t (float): Parameter of a near-intersection along ``curve2``.
        curve2 (.Curve): Second curve forming intersection.

    Returns:
        Tuple[float, float]: The refined parameters from a single Newton
        step.
    """
    # NOTE: We form -F(s, t) since we want to solve -DF^{-1} F(s, t).
    func_val = curve2.evaluate(t) - curve1.evaluate(s)
    if np.all(func_val == 0.0):
        # No refinement is needed.
        return s, t

    # NOTE: This assumes the curves are 2D.
    jac_mat = np.zeros((2, 2))
    # In curve.evaluate() and evaluate_hodograph() the roles of
    # columns and rows are swapped.
    nodes1 = curve1._nodes  # pylint: disable=protected-access
    jac_mat[0, :] = _evaluate_hodograph(nodes1, curve1.degree, s)
    nodes2 = curve2._nodes  # pylint: disable=protected-access
    jac_mat[1, :] = - _evaluate_hodograph(nodes2, curve2.degree, t)

    # Solve the system (via the transposes, since, as above, the roles
    # of columns and rows are reversed). Note that since ``func_val``
    # is a  1D object, then ``delta_s``, ``delta_t`` can be unpacked
    # without worry of them being vectors.
    delta_s, delta_t = np.linalg.solve(jac_mat.T, func_val.T)
    return s + delta_s, t + delta_t


def _cross_product(vec0, vec1):
    r"""Compute the cross-product of vectors in :math:`\mathbf{R}^2`.

    Utilizes the fact that

    .. math::

       \left[\begin{array}{c} A \\ B \\ 0 \end{array}\right] \times
           \left[\begin{array}{c} C \\ D \\ 0 \end{array}\right] =
           \left[\begin{array}{c} 0 \\ 0 \\ AD - BC \end{array}\right]

    and just returns the :math:`z` component.

    Args:
        vec0 (numpy.ndarray): A vector as a 1x2 NumPy array.
        vec1 (numpy.ndarray): A vector as a 1x2 NumPy array.

    Returns:
        float: The cross-product (or rather, its :math:`z` component).
    """
    return vec0[0, 0] * vec1[0, 1] - vec0[0, 1] * vec1[0, 0]


def segment_intersection(start0, end0, start1, end1, _fail=True):
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
       cross-product in :math:`\mathbf{R}^3` will always point in the
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

    .. testsetup::

       import numpy as np
       import bezier
       from bezier._intersection_helpers import segment_intersection

    .. doctest::
       :options: +NORMALIZE_WHITESPACE

       >>> start0 = np.array([[0.0, 0.0]])
       >>> end0 = np.array([[2.0, 2.0]])
       >>> start1 = np.array([[-1.0, 2.0]])
       >>> end1 = np.array([[1.0, 0.0]])
       >>> s, t = segment_intersection(start0, end0, start1, end1)
       >>> s
       0.25
       >>> t
       0.75

    Args:
        start0 (numpy.ndarray): A 1x2 NumPy array that is the start
            vector :math:`S_0` of the parametric line :math:`L_0(s)`.
        end0 (numpy.ndarray): A 1x2 NumPy array that is the end
            vector :math:`E_0` of the parametric line :math:`L_0(s)`.
        start1 (numpy.ndarray): A 1x2 NumPy array that is the start
            vector :math:`S_1` of the parametric line :math:`L_1(s)`.
        end1 (numpy.ndarray): A 1x2 NumPy array that is the end
            vector :math:`E_1` of the parametric line :math:`L_1(s)`.
        _fail (bool): Flag indicating if an exception should be raised
            when an unimplemented situation is encountered. If
            :data:`False`, then a non-sense answer will be returned.

    Returns:
        Tuple[float, float]: Pair of :math:`s_{\ast}` and :math:`t_{\ast}`
        such that the lines intersect:
        :math:`L_0\left(s_{\ast}\right) = L_1\left(t_{\ast}\right)`.

    Raises:
        NotImplementedError: If the lines are parallel (or one of the lines
            is degenerate). This manifests via
            :math:`\Delta_0 \times \Delta_1 = 0`.
    """
    delta0 = end0 - start0
    delta1 = end1 - start1
    cross_d0_d1 = _cross_product(delta0, delta1)
    if cross_d0_d1 == 0.0:
        if _fail:
            raise NotImplementedError('Delta_0 x Delta_1 = 0 not supported')
        else:
            return np.nan, np.nan
    else:
        start_delta = start1 - start0
        s = _cross_product(start_delta, delta1) / cross_d0_d1
        t = _cross_product(start_delta, delta0) / cross_d0_d1
        return s, t


def from_linearized(left, right, intersections):
    """Determine curve-curve intersection from pair of linearizations.

    If there is an intersection along the segments, adds that intersection
    to ``intersections``. Otherwise, returns without doing anything.

    Args:
        left (Linearization): First curve being intersected.
        right (Linearization): Second curve being intersected.
        intersections (list): A list of existing intersections.
    """
    s, t = segment_intersection(
        left.start_node, left.end_node,
        right.start_node, right.end_node)
    if not _in_interval(s, _WIGGLE_START, _WIGGLE_END):
        return
    if not _in_interval(t, _WIGGLE_START, _WIGGLE_END):
        return

    # Now, promote `s` and `t` onto the original curves.
    orig_s = (1 - s) * left.curve.start + s * left.curve.end
    orig_left = left.curve.root
    orig_t = (1 - t) * right.curve.start + t * right.curve.end
    orig_right = right.curve.root
    # Perform one step of Newton iteration to refine the computed
    # values of s and t.
    refined_s, refined_t = newton_refine(
        orig_s, orig_left, orig_t, orig_right)
    _check_parameters(refined_s, refined_t)
    at_endpoint = (refined_s in (left.curve.start, left.curve.end) or
                   refined_t in (right.curve.start, right.curve.end))
    intersection = Intersection(
        orig_left, refined_s, orig_right, refined_t,
        at_endpoint=at_endpoint)
    _add_intersection(intersection, intersections)


def _add_intersection(intersection, intersections):
    """Adds an intersection at bounding box tangency.

    If the intersection has already been found, does nothing.

    Args:
        intersection (Intersection): A new intersection to add.
        intersections (list): List of existing intersections.
    """
    # pylint: disable=protected-access
    for existing in intersections:
        if (existing.at_endpoint and
                existing.left is intersection.left and
                existing.right is intersection.right and
                existing._s_val == intersection._s_val and
                existing._t_val == intersection._t_val):
            return
    # pylint: enable=protected-access

    intersections.append(intersection)


def _tangent_bbox_intersection(left, right, intersections):
    r"""Check if two curves with tangent bounding boxes intersect.

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

       This method assumes callers will not pass curves that can be
       linearized / are linear. In :func:`all_intersections`, curves
       are pre-processed to do any linearization before the
       subdivision / intersection process begins.

    Args:
        left (.Curve): First curve being intersected.
        right (.Curve): Second curve being intersected.
        intersections (list): A list of already encountered
            intersections. If these curves intersect at their tangeny,
            then those intersections will be added to this list.
    """
    # pylint: disable=protected-access
    left_nodes = left._nodes
    right_nodes = right._nodes
    # pylint: enable=protected-access
    for i, s in ((0, 0.0), (-1, 1.0)):
        node_left = left_nodes[i, :]
        for j, t in ((0, 0.0), (-1, 1.0)):
            node_right = right_nodes[j, :]
            if _vector_close(node_left, node_right):
                orig_s = (1 - s) * left.start + s * left.end
                orig_t = (1 - t) * right.start + t * right.end
                intersection = Intersection(
                    left.root, orig_s, right.root, orig_t,
                    point=node_left, at_endpoint=True)
                _add_intersection(intersection, intersections)


# pylint: disable=too-many-return-statements
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
        nodes (numpy.ndarray): Point (Nx2) that determine a bounding box.
        line_start (numpy.ndarray): Beginning of a line segment (1x2).
        line_end (numpy.ndarray): End of a line segment (1x2).

    Returns:
        BoxIntersectionType: Enum indicating the type of bounding
        box intersection.
    """
    left, bottom = np.min(nodes, axis=0)
    right, top = np.max(nodes, axis=0)

    if (_in_interval(line_start[0, 0], left, right) and
            _in_interval(line_start[0, 1], bottom, top)):
        return BoxIntersectionType.intersection
    if (_in_interval(line_end[0, 0], left, right) and
            _in_interval(line_end[0, 1], bottom, top)):
        return BoxIntersectionType.intersection

    # NOTE: We pass ``_fail=False`` to ``segment_intersection`` below.
    #       At first, this may appear to "ignore" some potential
    #       intersections of parallel lines. However, no intersections
    #       will be missed. If parallel lines don't overlap, then
    #       there is nothing to miss. If they do overlap, then either
    #       the segment will have endpoints on the box (already covered)
    #       or the segment will contain an entire side of the box, which
    #       will force it to intersect the 3 edges that meet at the
    #       two ends of those sides. The parallel edge will be skipped,
    #       but the other two will be covered.

    # Bottom Edge
    s_bottom, t_bottom = segment_intersection(
        np.array([[left, bottom]]), np.array([[right, bottom]]),
        line_start, line_end, _fail=False)
    if _in_interval(s_bottom, 0.0, 1.0) and _in_interval(t_bottom, 0.0, 1.0):
        return BoxIntersectionType.intersection
    # Right Edge
    s_right, t_right = segment_intersection(
        np.array([[right, bottom]]), np.array([[right, top]]),
        line_start, line_end, _fail=False)
    if _in_interval(s_right, 0.0, 1.0) and _in_interval(t_right, 0.0, 1.0):
        return BoxIntersectionType.intersection
    # Top Edge
    s_top, t_top = segment_intersection(
        np.array([[right, top]]), np.array([[left, top]]),
        line_start, line_end, _fail=False)
    if _in_interval(s_top, 0.0, 1.0) and _in_interval(t_top, 0.0, 1.0):
        return BoxIntersectionType.intersection
    # NOTE: We skip the "last" edge. This is because any curve
    #       that doesn't have an endpoint on a curve must cross
    #       at least two, so we will already covered such curves
    #       in one of the branches above.

    return BoxIntersectionType.disjoint
# pylint: enable=too-many-return-statements


def intersect_one_round(candidates, intersections):
    """Perform one step of the intersection process.

    Checks if the bounding boxes of each pair in ``candidates``
    intersect. If the bounding boxes do not intersect, the pair
    is discarded. Otherwise, the pair is "accepted". Then we
    attempt to linearize each curve in an "accepted" pair and
    track the overall Linearization error for every curve
    encountered.

    Args:
        candidates (Union[list, itertools.chain]): An iterable of
            pairs of curves (or linearized curves).
        intersections (list): A list of already encountered
            intersections. If any intersections can be readily determined
            during this round of subdivision, then they will be added
            to this list.

    Returns:
        list: Returns a list of ``accepted`` pairs (among ``candidates``).
    """
    accepted = []

    for left, right in candidates:
        if isinstance(left, Linearization):
            if isinstance(right, Linearization):
                # Just jump right to line-line intersection.
                bbox_int = None
            else:
                right_nodes = right._nodes  # pylint: disable=protected-access
                bbox_int = bbox_line_intersect(
                    right_nodes, left.start_node, left.end_node)
        elif isinstance(right, Linearization):
            # pylint: disable=protected-access
            left_nodes = left._nodes
            # pylint: enable=protected-access
            bbox_int = bbox_line_intersect(
                left_nodes, right.start_node, right.end_node)
        else:
            left_nodes = left._nodes  # pylint: disable=protected-access
            right_nodes = right._nodes  # pylint: disable=protected-access
            bbox_int = bbox_intersect(left_nodes, right_nodes)

        if bbox_int is BoxIntersectionType.disjoint:
            continue
        elif bbox_int is BoxIntersectionType.tangent:
            _tangent_bbox_intersection(left, right, intersections)
            continue

        # Attempt to replace the curves with linearizations
        # if they are close enough to lines.
        # NOTE: This may be a wasted computation, e.g. if ``left``
        #       occurs in multiple accepted pairs. However, in practice
        #       the number of such pairs will be small so this cost
        #       will be low.
        left = Linearization.from_shape(left)
        # Now do the same for the right.
        right = Linearization.from_shape(right)

        # If both ``left`` and ``right`` are linearizations, then we can
        # intersect them immediately.
        if (isinstance(left, Linearization) and
                isinstance(right, Linearization)):
            from_linearized(left, right, intersections)
        else:
            # Add the accepted pair.
            accepted.append((left, right))

    return accepted


def all_intersections(candidates):
    r"""Find the points of intersection among pairs of curves.

    .. note::

       This assumes all curves in a candidate pair are in
       :math:`\mathbf{R}^2`, but does not **explicitly** check this.
       However, functions used here will fail if that assumption
       fails, e.g. :func:`bbox_intersect` and :func:`newton_refine`.

    Args:
        candidates (iterable): Iterable of pairs of curves that may
            intersect.

    Returns:
        list: List of all :class:`Intersection`s (possibly empty).

    Raises:
        ValueError: If the subdivision iteration does not terminate
            before exhausting the maximum number of subdivisions.
        NotImplementedError: If the subdivision process picks up too
            many candidate pairs. This typically indicates tangent
            curves or coincident curves.
    """
    # First make sure any curves that are linear / near-linear are
    # linearized (to avoid unnecessary checks, e.g. bbox intersect check).
    candidates = [
        (Linearization.from_shape(left), Linearization.from_shape(right))
        for left, right in candidates
    ]

    intersections = []
    for _ in six.moves.xrange(_MAX_INTERSECT_SUBDIVISIONS):
        accepted = intersect_one_round(candidates, intersections)
        if len(accepted) > 16:
            msg = _TOO_MANY_TEMPLATE.format(
                len(accepted), 4 * len(accepted))
            raise NotImplementedError(msg)

        # If none of the pairs have been accepted, then there is
        # no intersection.
        if not accepted:
            return intersections

        # If we **do** require more subdivisions, we need to update
        # the list of candidates.
        candidates = itertools.chain(*[
            itertools.product(left.subdivide(), right.subdivide())
            for left, right in accepted])

    raise ValueError(
        'Curve intersection failed to converge to approximately '
        'linear subdivisions after max iterations.',
        _MAX_INTERSECT_SUBDIVISIONS)


class BoxIntersectionType(enum.Enum):
    """Enum representing all possible bounding box intersections."""
    intersection = 'intersection'
    tangent = 'tangent'
    disjoint = 'disjoint'


class Linearization(object):
    """A linearization of a curve.

    This class is provided as a stand-in for a curve, so it
    provides a similar interface.

    Args:
        curve (.Curve): A curve that is linearized.
        error (Optional[float]): The linearization error. If not
            provided, this value is computed on the fly via
            :func:`linearization_error`.
    """

    def __init__(self, curve, error=None):
        self._curve = curve
        self._error = error

    def subdivide(self):
        """Do-nothing method to match the :class:`.Curve` interface.

        Returns:
            Tuple[Linearization]: List of all subdivided parts, which is
            just the current object.
        """
        return self,

    @property
    def curve(self):
        """.Curve: The curve that this linearization approximates."""
        return self._curve

    @property
    def error(self):
        """float: The linearization error for the linearized curve."""
        if self._error is None:
            self._error = linearization_error(self._curve)
        return self._error

    @property
    def start_node(self):
        """numpy.ndarray: The start vector of this linearization."""
        return self._curve._nodes[[0], :]  # pylint: disable=protected-access

    @property
    def end_node(self):
        """numpy.ndarray: The end vector of this linearization."""
        return self._curve._nodes[[-1], :]  # pylint: disable=protected-access

    @classmethod
    def from_shape(cls, shape):
        """Try to linearize a curve (or an already linearized curve).

        Args:
            shape (Union[.Curve, Linearization]): A curve or an already
                linearized curve.

        Returns:
            Union[.Curve, Linearization]: The (potentially linearized) curve.
        """
        if isinstance(shape, cls):
            return shape
        else:
            error = linearization_error(shape)
            if error == 0.0:
                err_exp = -np.inf
            else:
                _, err_exp = _FREXP(error)

            if err_exp <= _ERROR_EXPONENT:
                linearized = cls(shape, error=error)
                return linearized
            else:
                return shape


class Intersection(object):
    """Representation of a curve-curve intersection.

    Args:
        left (.Curve): The "left" curve in the intersection.
        s (float): The parameter along ``left`` where the
            intersection occurs.
        right (.Curve): The "right" curve in the intersection.
        t (float): The parameter along ``right`` where the
            intersection occurs.
        point (Optional[numpy.ndarray]): The point where the two
            curves actually intersect. If not provided, will be
            computed on the fly when first accessed.
        at_endpoint (bool): Indicates if the intersection happened at
            the endpoint of one (or both) of the two subdivided curves.
            This can be used to help with de-duplication of encountered
            intersections.
    """

    def __init__(self, left, s, right, t, point=None,
                 at_endpoint=False):
        self._left = left
        self._s_val = s
        self._right = right
        self._t_val = t
        self._point = point
        self.at_endpoint = at_endpoint

    @property
    def left(self):
        """numpy.ndarray: The "left" curve in the intersection."""
        return self._left

    @property
    def right(self):
        """numpy.ndarray: The "right" curve in the intersection."""
        return self._right

    @property
    def point(self):
        """numpy.ndarray: The point where the intersection occurs."""
        if self._point is None:
            self._point = _check_close(
                self._s_val, self._left, self._t_val, self._right)
        return self._point
