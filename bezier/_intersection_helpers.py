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
        NotImplementedError: If one of the values is equal to, or very
            close to, one of the endpoints of the unit interval.
    """
    if np.any(np.isclose([[s, t]], [[0.0], [1.0]])):
        raise NotImplementedError(
            'Intersection **at** endpoints not currently supported',
            'One parameter is very near 0 or 1', s, t)
    if s < 0.0 or s > 1.0:
        raise ValueError('s outside of unit interval', s)
    if t < 0.0 or t > 1.0:
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
        raise NotImplementedError('Delta_0 x Delta_1 = 0 not supported')
    else:
        start_delta = start1 - start0
        s = _cross_product(start_delta, delta1) / cross_d0_d1
        t = _cross_product(start_delta, delta0) / cross_d0_d1
        return s, t


def from_linearized(linearized_pairs):
    """Determine curve-curve intersections from pairs of linearizations.

    Args:
        linearized_pairs (list): List of pairs of :class:`Linearization`
            objects.

    Returns:
        list: List of all :class:`Intersection`s.
    """
    intersections = []
    for left, right in linearized_pairs:
        s, t = segment_intersection(
            left.start, left.end, right.start, right.end)
        _check_parameters(s, t)
        # Now, promote `s` and `t` onto the original curves.
        orig_s = (1 - s) * left.curve.start + s * left.curve.end
        orig_left = left.curve.root
        orig_t = (1 - t) * right.curve.start + t * right.curve.end
        orig_right = right.curve.root
        # Perform one step of Newton iteration to refine the computed
        # values of s and t.
        refined_s, refined_t = newton_refine(
            orig_s, orig_left, orig_t, orig_right)
        intersection = Intersection(
            orig_left, refined_s, orig_right, refined_t)
        intersections.append(intersection)

    return intersections


def intersect_one_round(candidates):
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

    Returns:
        Tuple[list, float]: Returns a list of ``accepted`` pairs
        (among ``candidates``) and the maximum linearization error
        among all curves in the list of accepted pairs.
    """
    accepted = []
    max_err = 0.0

    for left, right in candidates:
        # pylint: disable=protected-access
        left_nodes = left._nodes
        right_nodes = right._nodes
        # pylint: enable=protected-access
        if not bbox_intersect(left_nodes, right_nodes):
            continue

        # Attempt to replace the curves with linearizations
        # if they are close enough to lines.
        # NOTE: This may be a wasted computation, e.g. if ``left``
        #       occurs in multiple accepted pairs. However, in practice
        #       the number of such pairs will be small so this cost
        #       will be low.
        left, err_left = Linearization.from_shape(left)
        max_err = max(max_err, err_left)
        # Now do the same for the right.
        right, err_right = Linearization.from_shape(right)
        max_err = max(max_err, err_right)
        # Add the accepted pair.
        accepted.append((left, right))

    return accepted, max_err


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
    """
    for _ in six.moves.xrange(_MAX_INTERSECT_SUBDIVISIONS):
        accepted, max_err = intersect_one_round(candidates)

        # If none of the pairs have been accepted, then there is
        # no intersection.
        if not accepted:
            return []

        # In the case of ``accepted`` pairs, if the pairs are
        # sufficiently close to their linearizations, we can stop
        # the subdivisions and move on to the next step.
        _, max_exp = _FREXP(max_err)
        if max_exp <= _ERROR_EXPONENT:
            return from_linearized(accepted)

        # If we **do** require more subdivisions, we need to update
        # the list of candidates.
        candidates = itertools.chain(*[
            itertools.product(left.subdivide(), right.subdivide())
            for left, right in accepted])

    raise ValueError(
        'Curve intersection failed to converge to approximately '
        'linear subdivisions after max iterations.',
        _MAX_INTERSECT_SUBDIVISIONS)


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
    def _nodes(self):
        """numpy.ndarray: The nodes defining the linearized curve."""
        # NOTE: It's unclear if self._curve._nodes is appropriate here
        #       or if self._curve._nodes[[0, -1], :] is.
        return self._curve._nodes  # pylint: disable=protected-access

    @property
    def start(self):
        """numpy.ndarray: The start vector of this linearization."""
        return self._curve._nodes[[0], :]  # pylint: disable=protected-access

    @property
    def end(self):
        """numpy.ndarray: The end vector of this linearization."""
        return self._curve._nodes[[-1], :]  # pylint: disable=protected-access

    @classmethod
    def from_shape(cls, shape):
        """Try to linearize a curve (or an already linearized curve).

        Args:
            shape (Union[.Curve, Linearization]): A curve or an already
                linearized curve.

        Returns:
            Tuple[Union[.Curve, Linearization], float]: A pair of the
            (potentially linearized) curve and the linearization error.
        """
        if isinstance(shape, cls):
            return shape, shape.error
        else:
            error = linearization_error(shape)
            if error == 0.0:
                err_exp = -np.inf
            else:
                _, err_exp = _FREXP(error)

            if err_exp <= _ERROR_EXPONENT:
                linearized = cls(shape, error=error)
                return linearized, error
            else:
                return shape, error


class Intersection(object):
    """Representation of a curve-curve intersection.

    Args:
        left (.Curve): The "left" curve in the intersection.
        s (float): The parameter along ``left`` where the
            intersection occurs.
        right (.Curve): The "right" curve in the intersection.
        t (float): The parameter along ``right`` where the
            intersection occurs.
    """

    def __init__(self, left, s, right, t):
        self._left = left
        self._s_val = s
        self._right = right
        self._t_val = t
        self._point = None

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
