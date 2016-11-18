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
    r"""Refine a near-intersection using Newton's method.

    We assume we have an "almost" solution

    .. math::

       B_1\left(s_{\ast}\right) \approx B_2\left(t_{\ast}\right)

    and want to use Newton's method on the function

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

        This assumes ``curve1`` and ``curve2`` live in
        :math:`\mathbf{R}^2`.

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
        L_1(s) &= S_1 (1 - t) + E_1 t &&= S_1 + t \Delta_1.
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
            is degenerate. This manifests via
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
        numpy.ndarray: Array of all intersections.
    """
    intersections = []
    for left, right in linearized_pairs:
        s, t = segment_intersection(
            left.start, left.end, right.start, right.end)
        # TODO: Check if s, t are in [0, 1].
        # Now, promote `s` and `t` onto the original curves.
        orig_s = (1 - s) * left.curve.start + s * left.curve.end
        orig_left = left.curve.root
        orig_t = (1 - t) * right.curve.start + t * right.curve.end
        orig_right = right.curve.root
        # Perform one step of Newton iteration to refine the computed
        # values of s and t.
        refined_s, _ = newton_refine(
            orig_s, orig_left, orig_t, orig_right)
        # TODO: Check that (orig_left.evaluate(refined_s) ~=
        #                   orig_right.evaluate(refined_t))
        intersections.append(orig_left.evaluate(refined_s))

    return np.vstack(intersections)


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
        numpy.ndarray: Array of intersection points (possibly empty).

    Raises:
        ValueError: If the subdivision iteration does not terminate
            before exhausting the maximum number of subdivisions.
    """
    for _ in six.moves.xrange(_MAX_INTERSECT_SUBDIVISIONS):
        accepted, max_err = intersect_one_round(candidates)

        # If none of the pairs have been accepted, then there is
        # no intersection.
        if not accepted:
            return np.zeros((0, 2))

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

    return ValueError(
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
            _, err_exp = _FREXP(error)
            if err_exp <= _ERROR_EXPONENT:
                linearized = cls(shape, error=error)
                return linearized, error
            else:
                return shape, error
