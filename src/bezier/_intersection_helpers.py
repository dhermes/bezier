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
from bezier import _helpers
try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


# Set the threshold for exponent at half the bits available,
# this way one round of Newton's method can finish the job
# by squaring the error.
_ERROR_VAL = 0.5**26
_MAX_INTERSECT_SUBDIVISIONS = 20
_MAX_CANDIDATES = 16
_TOO_MANY_TEMPLATE = (
    'The number of candidate intersections is too high.\n'
    '{:d} accepted pairs gives {:d} candidate pairs.')
# Allow wiggle room for ``s`` and ``t`` computed during segment
# intersection. Any over- or under-shooting will (hopefully) be
# resolved in the Newton refinement step. If it isn't resolved, the
# call to _wiggle_interval() will fail the intersection.
_WIGGLE_START = -2.0**(-16)
_WIGGLE_END = 1.0 - _WIGGLE_START


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
    if not _helpers.vector_close(vec1, vec2):
        raise ValueError('B_1(s) and B_2(t) are not sufficiently close')

    return vec1


def _bbox_intersect(nodes1, nodes2):
    r"""Bounding box intersection predicate.

    Determines if the bounding box of two sets of control points
    intersects in :math:`\mathbf{R}^2` with non-trivial
    intersection (i.e. tangent bounding boxes are insufficient).

    .. note::

       Though we assume (and the code relies on this fact) that
       the nodes are two-dimensional, we don't check it.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes1 (numpy.ndarray): Set of control points for a
            B |eacute| zier shape.
        nodes2 (numpy.ndarray): Set of control points for a
            B |eacute| zier shape.

    Returns:
        int: Enum from ``BoxIntersectionType`` indicating the type of
        bounding box intersection.
    """
    left1, right1, bottom1, top1 = _helpers.bbox(nodes1)
    left2, right2, bottom2, top2 = _helpers.bbox(nodes2)

    if (right2 < left1 or right1 < left2 or
            top2 < bottom1 or top1 < bottom2):
        return BoxIntersectionType.DISJOINT

    if (right2 == left1 or right1 == left2 or
            top2 == bottom1 or top1 == bottom2):
        return BoxIntersectionType.TANGENT
    else:
        return BoxIntersectionType.INTERSECTION


def _linearization_error(nodes, degree):
    r"""Compute the maximum error of a linear approximation.

    Helper for :class:`.Linearization`, which is used during the
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

    .. image:: images/linearization_error.png
       :align: center

    .. testsetup:: linearization-error, linearization-error-fail

       import numpy as np
       import bezier
       from bezier._intersection_helpers import linearization_error

    .. doctest:: linearization-error

       >>> nodes = np.asfortranarray([
       ...     [0.0,  0.0],
       ...     [3.0,  1.0],
       ...     [9.0, -2.0],
       ... ])
       >>> linearization_error(nodes, 2)
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
       ...     [ 0.0,  0.0],
       ...     [ 5.0, 12.0],
       ...     [10.0, 24.0],
       ...     [30.0, 72.0],
       ... ])
       >>> linearization_error(nodes, 3)
       29.25

    Though it may seem that ``0`` is a more appropriate answer, consider
    the **goal** of this function. We seek to linearize curves and then
    intersect the linear approximations. Then the :math:`s`-values from
    the line-line intersection is lifted back to the curves. Thus
    the error :math:`\|B(s) - L(s)\|_2` is more relevant than the
    underyling algebraic curve containing :math:`B(s)`.

    .. note::

       It may be more appropriate to use a **relative** linearization error
       rather than the **absolute** error provided here. It's unclear if
       the domain :math:`\left[0, 1\right]` means the error is **already**
       adequately scaled or if the error should be scaled by the arc
       length of the curve or the (easier-to-compute) length of the line.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Nodes of a curve.
        degree (int): The degree of the curve

    Returns:
        float: The maximum error between the curve and the
        linear approximation.
    """
    if degree == 1:
        return 0.0

    second_deriv = nodes[:-2, :] - 2.0 * nodes[1:-1, :] + nodes[2:, :]
    worst_case = np.max(np.abs(second_deriv), axis=0)

    # max_{0 <= s <= 1} s(1 - s)/2 = 1/8 = 0.125
    multiplier = 0.125 * degree * (degree - 1)
    # NOTE: worst_case is 1D due to np.max(), so this is the vector norm.
    return multiplier * np.linalg.norm(worst_case, ord=2)


def _newton_refine(s, nodes1, t, nodes2):
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

       This implementation assumes the curves live in
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

    .. image:: images/newton_refine1.png
       :align: center

    .. testsetup:: newton-refine1, newton-refine2, newton-refine3

       import numpy as np
       import bezier
       from bezier._intersection_helpers import newton_refine

    .. doctest:: newton-refine1

       >>> nodes1 = np.asfortranarray([
       ...     [0.0, 0.0],
       ...     [2.0, 4.0],
       ...     [4.0, 0.0],
       ... ])
       >>> nodes2 = np.asfortranarray([
       ...     [2.0, 0.0],
       ...     [0.0, 3.0],
       ... ])
       >>> s, t = 0.375, 0.25
       >>> new_s, new_t = newton_refine(s, nodes1, t, nodes2)
       >>> 64.0 * (new_s - s)
       -9.0
       >>> 64.0 * (new_t - t)
       18.0

    .. testcleanup:: newton-refine1

       import make_images
       curve1 = bezier.Curve(nodes1, degree=2)
       curve2 = bezier.Curve(nodes2, degree=1)
       make_images.newton_refine1(s, new_s, curve1, t, new_t, curve2)

    For "typical" curves, we converge to a solution quadratically.
    This means that the number of correct digits doubles every
    iteration (until machine precision is reached).

    .. image:: images/newton_refine2.png
       :align: center

    .. doctest:: newton-refine2

       >>> nodes1 = np.asfortranarray([
       ...     [0.0 ,  0.0],
       ...     [0.25,  2.0],
       ...     [0.5 , -2.0],
       ...     [0.75,  2.0],
       ...     [1.0 ,  0.0],
       ... ])
       >>> nodes2 = np.asfortranarray([
       ...     [0.0 , 1.0],
       ...     [0.25, 0.5],
       ...     [0.5 , 0.5],
       ...     [0.75, 0.5],
       ...     [1.0 , 0.0],
       ... ])
       >>> # The expected intersection is the only real root of
       >>> # 28 s^3 - 30 s^2 + 9 s - 1.
       >>> omega = (28.0 * np.sqrt(17.0) + 132.0)**(1.0 / 3.0) / 28.0
       >>> expected = 5.0 / 14.0 + omega + 1 / (49.0 * omega)
       >>> s_vals = [0.625, None, None, None, None]
       >>> t = 0.625
       >>> np.log2(abs(expected - s_vals[0]))
       -4.399...
       >>> s_vals[1], t = newton_refine(s_vals[0], nodes1, t, nodes2)
       >>> np.log2(abs(expected - s_vals[1]))
       -7.901...
       >>> s_vals[2], t = newton_refine(s_vals[1], nodes1, t, nodes2)
       >>> np.log2(abs(expected - s_vals[2]))
       -16.010...
       >>> s_vals[3], t = newton_refine(s_vals[2], nodes1, t, nodes2)
       >>> np.log2(abs(expected - s_vals[3]))
       -32.110...
       >>> s_vals[4], t = newton_refine(s_vals[3], nodes1, t, nodes2)
       >>> s_vals[4] == expected
       True

    .. testcleanup:: newton-refine2

       import make_images
       curve1 = bezier.Curve(nodes1, degree=4)
       curve2 = bezier.Curve(nodes2, degree=4)
       make_images.newton_refine2(s_vals, curve1, curve2)

    However, when the intersection occurs at a point of tangency,
    the convergence becomes linear. This means that the number of
    correct digits added each iteration is roughly constant.

    .. image:: images/newton_refine3.png
       :align: center

    .. doctest:: newton-refine3

       >>> nodes1 = np.asfortranarray([
       ...     [0.0, 0.0],
       ...     [0.5, 1.0],
       ...     [1.0, 0.0],
       ... ])
       >>> nodes2 = np.asfortranarray([
       ...     [0.0, 0.5],
       ...     [1.0, 0.5],
       ... ])
       >>> expected = 0.5
       >>> s_vals = [0.375, None, None, None, None, None]
       >>> t = 0.375
       >>> np.log2(abs(expected - s_vals[0]))
       -3.0
       >>> s_vals[1], t = newton_refine(s_vals[0], nodes1, t, nodes2)
       >>> np.log2(abs(expected - s_vals[1]))
       -4.0
       >>> s_vals[2], t = newton_refine(s_vals[1], nodes1, t, nodes2)
       >>> np.log2(abs(expected - s_vals[2]))
       -5.0
       >>> s_vals[3], t = newton_refine(s_vals[2], nodes1, t, nodes2)
       >>> np.log2(abs(expected - s_vals[3]))
       -6.0
       >>> s_vals[4], t = newton_refine(s_vals[3], nodes1, t, nodes2)
       >>> np.log2(abs(expected - s_vals[4]))
       -7.0
       >>> s_vals[5], t = newton_refine(s_vals[4], nodes1, t, nodes2)
       >>> np.log2(abs(expected - s_vals[5]))
       -8.0

    .. testcleanup:: newton-refine3

       import make_images
       curve1 = bezier.Curve(nodes1, degree=2)
       curve2 = bezier.Curve(nodes2, degree=1)
       make_images.newton_refine3(s_vals, curve1, curve2)

    Unfortunately, the process terminates with an error that is not close
    to machine precision :math:`\varepsilon` when
    :math:`\Delta s = \Delta t = 0`.

    .. testsetup:: newton-refine3-continued

       import numpy as np
       import bezier
       from bezier._intersection_helpers import newton_refine

       nodes1 = np.asfortranarray([
           [0.0, 0.0],
           [0.5, 1.0],
           [1.0, 0.0],
       ])
       nodes2 = np.asfortranarray([
           [0.0, 0.5],
           [1.0, 0.5],
       ])

    .. doctest:: newton-refine3-continued

       >>> s1 = t1 = 0.5 - 0.5**27
       >>> np.log2(0.5 - s1)
       -27.0
       >>> s2, t2 = newton_refine(s1, nodes1, t1, nodes2)
       >>> s2 == t2
       True
       >>> np.log2(0.5 - s2)
       -28.0
       >>> s3, t3 = newton_refine(s2, nodes1, t2, nodes2)
       >>> s3 == t3 == s2
       True

    Due to round-off near the point of tangency, the final error
    resembles :math:`\sqrt{\varepsilon}` rather than machine
    precision as expected.

    .. note::

       The following is not implemented in this function. It's just
       an exploration on how the shortcomings might be addressed.

    However, this can be overcome. At the point of tangency, we want
    :math:`B_1'(s) \parallel B_2'(t)`. This can be checked numerically via

    .. math::

        B_1'(s) \times B_2'(t) = 0.

    For the last example (the one that converges linearly), this is

    .. math::

        0 = \left[\begin{array}{c} 1 \\ 2 - 4s \end{array}\right] \times
            \left[\begin{array}{c} 1 \\ 0 \end{array}\right] = 4 s - 2.

    With this, we can modify Newton's method to find a zero of the
    over-determined system

    .. math::

        G(s, t) = \left[\begin{array}{c} B_0(s) - B_1(t) \\
            B_1'(s) \times B_2'(t) \end{array}\right] =
            \left[\begin{array}{c} s - t \\ 2 s (1 - s) - \frac{1}{2} \\
            4 s - 2\end{array}\right].

    Since :math:`DG` is :math:`3 \times 2`, we can't invert it. However,
    we can find a least-squares solution:

    .. math::

        \left(DG^T DG\right) \left[\begin{array}{c}
            \Delta s \\ \Delta t \end{array}\right] = -DG^T G.

    This only works if :math:`DG` has full rank. In this case, it does
    since the submatrix containing the first and last rows has rank two:

    .. math::

        DG = \left[\begin{array}{c c} 1 & -1 \\
            2 - 4 s & 0 \\
            4 & 0 \end{array}\right].

    Though this avoids a singular system, the normal equations have a
    condition number that is the square of the condition number of the matrix.

    Starting from :math:`s = t = \frac{3}{8}` as above:

    .. testsetup:: newton-refine4

       import numpy as np
       from bezier import _helpers

       def modified_update(s, t):
           minus_G = np.asfortranarray([
               [t - s],
               [0.5 - 2.0 * s * (1.0 - s)],
               [2.0 - 4.0 * s],
           ])
           DG = np.asfortranarray([
               [1.0, -1.0],
               [2.0 - 4.0 * s, 0.0],
               [4.0, 0.0],
           ])
           DG_t = np.asfortranarray(DG.T)

           LHS = _helpers.matrix_product(DG_t, DG)
           RHS = _helpers.matrix_product(DG_t, minus_G)
           delta_params = np.linalg.solve(LHS, RHS)
           delta_s, delta_t = delta_params.flatten()
           return s + delta_s, t + delta_t

    .. doctest:: newton-refine4

       >>> s0, t0 = 0.375, 0.375
       >>> np.log2(0.5 - s0)
       -3.0
       >>> s1, t1 = modified_update(s0, t0)
       >>> s1 == t1
       True
       >>> 1040.0 * s1
       519.0
       >>> np.log2(0.5 - s1)
       -10.022...
       >>> s2, t2 = modified_update(s1, t1)
       >>> s2 == t2
       True
       >>> np.log2(0.5 - s2)
       -31.067...
       >>> s3, t3 = modified_update(s2, t2)
       >>> s3 == t3 == 0.5
       True

    Args:
        s (float): Parameter of a near-intersection along the first curve.
        nodes1 (numpy.ndarray): Nodes of first curve forming intersection.
        t (float): Parameter of a near-intersection along the second curve.
        nodes2 (numpy.ndarray): Nodes of second curve forming intersection.

    Returns:
        Tuple[float, float]: The refined parameters from a single Newton
        step.
    """
    # NOTE: We form -F(s, t) since we want to solve -DF^{-1} F(s, t).
    func_val = (
        _curve_helpers.evaluate_multi(nodes2, np.asfortranarray([t])) -
        _curve_helpers.evaluate_multi(nodes1, np.asfortranarray([s])))
    if np.all(func_val == 0.0):
        # No refinement is needed.
        return s, t

    # NOTE: This assumes nodes are (d + 1)x2.
    degree1 = nodes1.shape[0] - 1
    degree2 = nodes2.shape[0] - 1

    # NOTE: This assumes the curves are 2D.
    jac_mat = np.empty((2, 2), order='F')
    # In curve.evaluate() and evaluate_hodograph() the roles of
    # columns and rows are swapped.
    jac_mat[0, :] = _curve_helpers.evaluate_hodograph(s, nodes1, degree1)
    jac_mat[1, :] = - _curve_helpers.evaluate_hodograph(t, nodes2, degree2)

    # Solve the system (via the transposes, since, as above, the roles
    # of columns and rows are reversed).
    result = np.linalg.solve(jac_mat.T, func_val.T)
    # Convert to row-vector and unpack (also makes assertion on shape).
    (delta_s, delta_t), = result.T
    return s + delta_s, t + delta_t


def _segment_intersection(start0, end0, start1, end1):
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

    .. image:: images/segment_intersection1.png
       :align: center

    .. testsetup:: segment-intersection1, segment-intersection2

       import numpy as np
       from bezier._intersection_helpers import segment_intersection

    .. doctest:: segment-intersection1
       :options: +NORMALIZE_WHITESPACE

       >>> start0 = np.asfortranarray([[0.0, 0.0]])
       >>> end0 = np.asfortranarray([[2.0, 2.0]])
       >>> start1 = np.asfortranarray([[-1.0, 2.0]])
       >>> end1 = np.asfortranarray([[1.0, 0.0]])
       >>> s, t, _ = segment_intersection(start0, end0, start1, end1)
       >>> s
       0.25
       >>> t
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

    .. image:: images/segment_intersection2.png
       :align: center

    .. doctest:: segment-intersection2
       :options: +NORMALIZE_WHITESPACE

       >>> start0 = np.asfortranarray([[1.0, 0.0]])
       >>> end0 = np.asfortranarray([[0.0, 1.0]])
       >>> start1 = np.asfortranarray([[-1.0, 3.0]])
       >>> end1 = np.asfortranarray([[3.0, -1.0]])
       >>> _, _, success = segment_intersection(start0, end0, start1, end1)
       >>> bool(success)
       False

    .. testcleanup:: segment-intersection2

       import make_images
       make_images.segment_intersection2(start0, end0, start1, end1)

    Instead, we use :func:`.parallel_different`:

    .. testsetup:: segment-intersection2-continued

       import numpy as np
       from bezier._intersection_helpers import parallel_different

       start0 = np.asfortranarray([[1.0, 0.0]])
       end0 = np.asfortranarray([[0.0, 1.0]])
       start1 = np.asfortranarray([[-1.0, 3.0]])
       end1 = np.asfortranarray([[3.0, -1.0]])

    .. doctest:: segment-intersection2-continued

       >>> bool(parallel_different(start0, end0, start1, end1))
       True

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

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
        Tuple[float, float, bool]: Pair of :math:`s_{\ast}` and
        :math:`t_{\ast}` such that the lines intersect:
        :math:`L_0\left(s_{\ast}\right) = L_1\left(t_{\ast}\right)` and then
        a boolean indicating if an intersection was found.
    """
    delta0 = end0 - start0
    delta1 = end1 - start1
    cross_d0_d1 = _helpers.cross_product(delta0, delta1)
    if cross_d0_d1 == 0.0:
        return None, None, False
    else:
        start_delta = start1 - start0
        s = _helpers.cross_product(start_delta, delta1) / cross_d0_d1
        t = _helpers.cross_product(start_delta, delta0) / cross_d0_d1
        return s, t, True


def _parallel_different(start0, end0, start1, end1):
    r"""Checks if two parallel lines ever meet.

    Meant as a back-up when :func:`.segment_intersection` fails.

    .. note::

       This function assumes but never verifies that the lines
       are parallel.

    In the case that the segments are parallel and lie on the **exact**
    same line, finding a unique intersection is not possible. However, if
    they are parallel but on **different** lines, then there is a
    **guarantee** of no intersection.

    In :func:`.segment_intersection`, we utilized the normal form of the
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

    .. image:: images/parallel_different1.png
       :align: center

    .. testsetup:: parallel-different1, parallel-different2

       import numpy as np
       from bezier._intersection_helpers import parallel_different

    .. doctest:: parallel-different1

       >>> # Line: y = 1
       >>> start0 = np.asfortranarray([[0.0, 1.0]])
       >>> end0 = np.asfortranarray([[1.0, 1.0]])
       >>> # Vertical shift up: y = 2
       >>> start1 = np.asfortranarray([[-1.0, 2.0]])
       >>> end1 = np.asfortranarray([[3.0, 2.0]])
       >>> bool(parallel_different(start0, end0, start1, end1))
       True

    .. testcleanup:: parallel-different1

       import make_images
       make_images.helper_parallel_different(
           start0, end0, start1, end1, 'parallel_different1.png')

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

    .. image:: images/parallel_different2.png
       :align: center

    .. doctest:: parallel-different2

       >>> start0 = np.asfortranarray([[1.0, 0.0]])
       >>> delta0 = np.asfortranarray([[2.0, -1.0]])
       >>> end0 = start0 + 1.0 * delta0
       >>> start1 = start0 + 1.5 * delta0
       >>> end1 = start0 + 2.0 * delta0
       >>> bool(parallel_different(start0, end0, start1, end1))
       True

    .. testcleanup:: parallel-different2

       import make_images
       make_images.helper_parallel_different(
           start0, end0, start1, end1, 'parallel_different2.png')

    but if the intervals overlap, like :math:`\left[0, 1\right]` and
    :math:`\left[-1, \frac{1}{2}\right]`, the segments meet:

    .. image:: images/parallel_different3.png
       :align: center

    .. testsetup:: parallel-different3, parallel-different4

       import numpy as np
       from bezier._intersection_helpers import parallel_different

       start0 = np.asfortranarray([[1.0, 0.0]])
       delta0 = np.asfortranarray([[2.0, -1.0]])
       end0 = start0 + 1.0 * delta0

    .. doctest:: parallel-different3

       >>> start1 = start0 - 1.0 * delta0
       >>> end1 = start0 + 0.5 * delta0
       >>> bool(parallel_different(start0, end0, start1, end1))
       False

    .. testcleanup:: parallel-different3

       import make_images
       make_images.helper_parallel_different(
           start0, end0, start1, end1, 'parallel_different3.png')

    Similarly, if the second interval completely contains the first,
    the segments meet:

    .. image:: images/parallel_different4.png
       :align: center

    .. doctest:: parallel-different4

       >>> start1 = start0 + 3.0 * delta0
       >>> end1 = start0 - 2.0 * delta0
       >>> bool(parallel_different(start0, end0, start1, end1))
       False

    .. testcleanup:: parallel-different4

       import make_images
       make_images.helper_parallel_different(
           start0, end0, start1, end1, 'parallel_different4.png')

    .. note::

       This function doesn't currently allow wiggle room around the
       desired value, i.e. the two values must be bitwise identical.
       However, the most "correct" version of this function likely
       should allow for some round off.

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
        bool: Indicating if the lines are different.
    """
    delta0 = end0 - start0
    line0_const = _helpers.cross_product(start0, delta0)
    start1_against = _helpers.cross_product(start1, delta0)
    if line0_const != start1_against:
        return True

    # Each array is 1x2 (i.e. a row vector), we want the vector dot product.
    norm0_sq = np.vdot(delta0[0, :], delta0[0, :])
    start_numer = np.vdot((start1 - start0)[0, :], delta0[0, :])
    #      0 <= start_numer / norm0_sq <= 1
    # <==> 0 <= start_numer            <= norm0_sq
    if _helpers.in_interval(start_numer, 0.0, norm0_sq):
        return False

    end_numer = np.vdot((end1 - start0)[0, :], delta0[0, :])
    #      0 <= end_numer / norm0_sq <= 1
    # <==> 0 <= end_numer            <= norm0_sq
    if _helpers.in_interval(end_numer, 0.0, norm0_sq):
        return False

    # We know neither the start or end parameters are in [0, 1], but
    # they may contain [0, 1] between them.
    min_val, max_val = sorted([start_numer, end_numer])
    # So we make sure that 0 isn't between them.
    return not _helpers.in_interval(0.0, min_val, max_val)


# pylint: disable=too-many-arguments
def _from_linearized_low_level_py(
        error1, start1, end1, start_node1, end_node1, nodes1,
        error2, start2, end2, start_node2, end_node2, nodes2):
    """Determine curve-curve intersection from pair of linearizations.

    The inputs are the "fully-unpacked" values from two
    :class:`.Linearization` objects: ``first`` and ``second``.

    Args:
        error1 (float): The linearization error for the first curve.
        start1 (float): The start parameter of the (potentially
            subdivided) first curve.
        end1 (float): The end parameter of the (potentially
            subdivided) first curve.
        start_node1 (numpy.ndarray): The (``1x2``) start node in the
            first curve (also in the linearization).
        end_node1 (numpy.ndarray): The (``1x2``) end node in the
            first curve (also in the linearization).
        nodes1 (numpy.ndarray): The (``(D + 1)x2``) nodes of the "root"
            curve that contains the first curve.
        error2 (float): The linearization error for the second curve.
        start2 (float): The start parameter of the (potentially
            subdivided) second curve.
        end2 (float): The end parameter of the (potentially
            subdivided) second curve.
        start_node2 (numpy.ndarray): The (``1x2``) start node in the
            second curve (also in the linearization).
        end_node2 (numpy.ndarray): The (``1x2``) end node in the
            second curve (also in the linearization).
        nodes2 (numpy.ndarray): The (``(D2 + 1)x2``) nodes of the "root"
            curve that contains the second curve.

    Returns:
        Tuple[float, float, bool]: Triple of:

        * ``s``-coordinate of intersection
        * ``t``-coordinate of intersection
        * Flag indicating of the linearizations actually intersect

    Raises:
        NotImplementedError: If the segment intersection fails.
    """
    # pylint: disable=too-many-locals
    s, t, success = segment_intersection(
        start_node1, end_node1, start_node2, end_node2)
    if success:
        if error1 == 0.0 and not _helpers.in_interval(s, 0.0, 1.0):
            return None, None, False
        if error2 == 0.0 and not _helpers.in_interval(t, 0.0, 1.0):
            return None, None, False
        if not _helpers.in_interval(s, _WIGGLE_START, _WIGGLE_END):
            return None, None, False
        if not _helpers.in_interval(t, _WIGGLE_START, _WIGGLE_END):
            return None, None, False
    else:
        # Handle special case where the curves are actually lines.
        if error1 == 0.0 and error2 == 0.0:
            if parallel_different(start_node1, end_node1,
                                  start_node2, end_node2):
                return None, None, False

        raise NotImplementedError('Line segments parallel.')

    # Now, promote `s` and `t` onto the original curves.
    orig_s = (1 - s) * start1 + s * end1
    orig_t = (1 - t) * start2 + t * end2
    # Perform one step of Newton iteration to refine the computed
    # values of s and t.
    refined_s, refined_t = newton_refine(orig_s, nodes1, orig_t, nodes2)
    refined_s = _helpers.wiggle_interval(refined_s)
    refined_t = _helpers.wiggle_interval(refined_t)
    return refined_s, refined_t, True
    # pylint: enable=too-many-locals
# pylint: enable=too-many-arguments


def from_linearized(first, second, intersections):
    """Determine curve-curve intersection from pair of linearizations.

    If there is an intersection along the segments, adds that intersection
    to ``intersections``. Otherwise, returns without doing anything.

    Args:
        first (Linearization): First curve being intersected.
        second (Linearization): Second curve being intersected.
        intersections (list): A list of existing intersections.
    """
    curve1 = first.curve
    curve2 = second.curve
    # pylint: disable=protected-access
    orig_first = curve1._root
    orig_second = curve2._root
    refined_s, refined_t, success = _from_linearized_low_level(
        first.error, curve1._start, curve1._end, first.start_node,
        first.end_node, orig_first._nodes,
        second.error, curve2._start, curve2._end, second.start_node,
        second.end_node, orig_second._nodes)
    # pylint: enable=protected-access
    if success:
        intersection = Intersection(
            orig_first, refined_s, orig_second, refined_t)
        _add_intersection(intersection, intersections)


def _add_intersection(intersection, intersections):
    """Adds an intersection to list of ``intersections``.

    Accounts for repeated points at curve endpoints. If the
    intersection has already been found, does nothing.

    Args:
        intersection (Intersection): A new intersection to add.
        intersections (list): List of existing intersections.
    """
    for existing in intersections:
        if (existing.first is intersection.first and
                existing.second is intersection.second and
                _helpers.ulps_away(existing.s, intersection.s) and
                _helpers.ulps_away(existing.t, intersection.t)):
            return

    intersections.append(intersection)


def _endpoint_check(first, node_first, s,
                    second, node_second, t, intersections):
    r"""Check if curve endpoints are identical.

    Helper for :func:`_tangent_bbox_intersection`.

    Args:
        first (.Curve): First curve being intersected (assumed in
            :math:\mathbf{R}^2`).
        node_first (numpy.ndarray): ``1x2`` array, one of the endpoints
            of ``first``.
        s (float): The parameter corresponding to ``node_first``, so
             expected to be one of ``0.0`` or ``1.0``.
        second (.Curve): Second curve being intersected (assumed in
            :math:\mathbf{R}^2`).
        node_second (numpy.ndarray): ``1x2`` array, one of the endpoints
            of ``second``.
        t (float): The parameter corresponding to ``node_second``, so
             expected to be one of ``0.0`` or ``1.0``.
        intersections (list): A list of already encountered
            intersections. If these curves intersect at their tangeny,
            then those intersections will be added to this list.
    """
    if _helpers.vector_close(node_first, node_second):
        # pylint: disable=protected-access
        orig_s = (1 - s) * first._start + s * first._end
        orig_t = (1 - t) * second._start + t * second._end
        # ``node_first`` is a row in a Fortran-contiguous array, so
        # the values aren't contiguous.
        intersection = Intersection(
            first._root, orig_s, second._root, orig_t,
            point=np.asfortranarray(node_first))
        # pylint: enable=protected-access
        _add_intersection(intersection, intersections)


def _tangent_bbox_intersection(first, second, intersections):
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

       This function assumes callers will not pass curves that can be
       linearized / are linear. In :func:`all_intersections`, curves
       are pre-processed to do any linearization before the
       subdivision / intersection process begins.

    Args:
        first (.Curve): First curve being intersected (assumed in
            :math:\mathbf{R}^2`).
        second (.Curve): Second curve being intersected (assumed in
            :math:\mathbf{R}^2`).
        intersections (list): A list of already encountered
            intersections. If these curves intersect at their tangeny,
            then those intersections will be added to this list.
    """
    # NOTE: We want the nodes to be 1x2 but accessing
    #       ``first_nodes[[index], :]`` makes a copy while the
    #       accesses below **do not** copy. See
    #       (https://docs.scipy.org/doc/numpy-1.6.0/reference/
    #        arrays.indexing.html#advanced-indexing)
    #       We don't use order='F' when reshaping since the data isn't
    #       contiguous (since we assumes `_nodes` is Fortran contiguous).
    node_first1 = first._nodes[0, :].reshape((1, 2))
    node_first2 = first._nodes[-1, :].reshape((1, 2))
    node_second1 = second._nodes[0, :].reshape((1, 2))
    node_second2 = second._nodes[-1, :].reshape((1, 2))

    _endpoint_check(
        first, node_first1, 0.0, second, node_second1, 0.0, intersections)
    _endpoint_check(
        first, node_first1, 0.0, second, node_second2, 1.0, intersections)
    _endpoint_check(
        first, node_first2, 1.0, second, node_second1, 0.0, intersections)
    _endpoint_check(
        first, node_first2, 1.0, second, node_second2, 1.0, intersections)


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
        nodes (numpy.ndarray): Points (Nx2) that determine a bounding box.
        line_start (numpy.ndarray): Beginning of a line segment (1x2).
        line_end (numpy.ndarray): End of a line segment (1x2).

    Returns:
        int: Enum from ``BoxIntersectionType`` indicating the type of
        bounding box intersection.
    """
    left, right, bottom, top = _helpers.bbox(nodes)

    if (_helpers.in_interval(line_start[0, 0], left, right) and
            _helpers.in_interval(line_start[0, 1], bottom, top)):
        return BoxIntersectionType.INTERSECTION
    if (_helpers.in_interval(line_end[0, 0], left, right) and
            _helpers.in_interval(line_end[0, 1], bottom, top)):
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
        np.asfortranarray([[left, bottom]]),
        np.asfortranarray([[right, bottom]]),
        line_start, line_end)
    if (success and _helpers.in_interval(s_bottom, 0.0, 1.0) and
            _helpers.in_interval(t_bottom, 0.0, 1.0)):
        return BoxIntersectionType.INTERSECTION
    # Right Edge
    s_right, t_right, success = segment_intersection(
        np.asfortranarray([[right, bottom]]),
        np.asfortranarray([[right, top]]),
        line_start, line_end)
    if (success and _helpers.in_interval(s_right, 0.0, 1.0) and
            _helpers.in_interval(t_right, 0.0, 1.0)):
        return BoxIntersectionType.INTERSECTION
    # Top Edge
    s_top, t_top, success = segment_intersection(
        np.asfortranarray([[right, top]]),
        np.asfortranarray([[left, top]]),
        line_start, line_end)
    if (success and _helpers.in_interval(s_top, 0.0, 1.0) and
            _helpers.in_interval(t_top, 0.0, 1.0)):
        return BoxIntersectionType.INTERSECTION
    # NOTE: We skip the "last" edge. This is because any curve
    #       that doesn't have an endpoint on a curve must cross
    #       at least two, so we will already covered such curves
    #       in one of the branches above.

    return BoxIntersectionType.DISJOINT


def intersect_one_round(candidates, intersections):
    """Perform one step of the intersection process.

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
        list: Returns a list of ``accepted`` pairs (among ``candidates``).
    """
    accepted = []

    # NOTE: In the below we replace ``isinstance(a, B)`` with
    #       ``a.__class__ is B``, which is a 3-3.5x speedup.
    for first, second in candidates:
        if first.__class__ is Linearization:
            if second.__class__ is Linearization:
                # If both ``first`` and ``second`` are linearizations, then
                # we can intersect them immediately.
                from_linearized(first, second, intersections)
                continue
            else:
                bbox_int = bbox_line_intersect(
                    second._nodes, first.start_node, first.end_node)
        elif second.__class__ is Linearization:
            bbox_int = bbox_line_intersect(
                first._nodes, second.start_node, second.end_node)
        else:
            bbox_int = bbox_intersect(first._nodes, second._nodes)

        if bbox_int == BoxIntersectionType.DISJOINT:
            continue
        elif bbox_int == BoxIntersectionType.TANGENT:
            _tangent_bbox_intersection(first, second, intersections)
            continue

        # If we haven't ``continue``-d, add the accepted pair.
        accepted.append((first, second))

    return accepted


def _next_candidates(first, second):
    """Take a pair of "accepted" curves and subdivide them.

    Attempts to replace the subdivided curves with linearizations
    if they are close enough to lines.

    Args:
        first (Union[.Curve, Linearization]): First curve in pair.
        second (Union[.Curve, Linearization]): Second curve in pair.

    Returns:
        itertools.product: Iterator of pairs of first and second curves after
        subdivision, some of which may be linearized.
    """
    # NOTE: This may be a wasted computation, e.g. if ``first``
    #       or ``second`` occur in multiple accepted pairs (the caller
    #       only passes one pair at a time). However, in practice
    #       the number of such pairs will be small so this cost
    #       will be low.
    lin1 = six.moves.map(Linearization.from_shape, first.subdivide())
    lin2 = six.moves.map(Linearization.from_shape, second.subdivide())
    return itertools.product(lin1, lin2)


class IntersectionStrategy(enum.Enum):
    """Enum determining if the type of intersection algorithm to use."""

    geometric = 'geometric'
    """Geometric approach to intersection (via subdivision)."""

    algebraic = 'algebraic'
    """Algebraic approach to intersection (via implicitization)."""


def _all_intersections_geometric(candidates):
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
        (Linearization.from_shape(first), Linearization.from_shape(second))
        for first, second in candidates
    ]

    intersections = []
    for _ in six.moves.xrange(_MAX_INTERSECT_SUBDIVISIONS):
        accepted = intersect_one_round(candidates, intersections)
        if len(accepted) > _MAX_CANDIDATES:
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
            _next_candidates(first, second)
            for first, second in accepted])

    raise ValueError(
        'Curve intersection failed to converge to approximately '
        'linear subdivisions after max iterations.',
        _MAX_INTERSECT_SUBDIVISIONS)


def _all_intersections_algebraic(candidates):
    r"""Find the points of intersection among pairs of curves.

    .. note::

       This assumes all curves in a candidate pair are in
       :math:`\mathbf{R}^2`, but does not **explicitly** check this.
       However, functions used here will fail if that assumption
       fails.

    Args:
        candidates (iterable): Iterable of pairs of curves that may
            intersect.

    Returns:
        list: List of all :class:`Intersection`s (possibly empty).
    """
    # NOTE: This cyclic import need be resolved.
    from bezier import _implicitization

    result = []
    for first, second in candidates:
        # NOTE: In the below we replace ``isinstance(a, B)`` with
        #       ``a.__class__ is B``, which is a 3-3.5x speedup.
        curve1 = first.curve if first.__class__ is Linearization else first
        nodes1 = curve1._nodes
        curve2 = second.curve if second.__class__ is Linearization else second
        nodes2 = curve2._nodes
        # Only attempt this if the bounding boxes intersect.
        if bbox_intersect(nodes1, nodes2) == BoxIntersectionType.DISJOINT:
            continue

        st_vals = _implicitization.intersect_curves(nodes1, nodes2)
        for s, t in st_vals:
            intersection = Intersection(curve1, s, curve2, t)
            result.append(intersection)

    return result


def all_intersections(candidates, strategy=IntersectionStrategy.geometric):
    r"""Find the points of intersection among pairs of curves.

    .. note::

       This assumes all curves in a candidate pair are in
       :math:`\mathbf{R}^2`, but does not **explicitly** check this.
       However, functions used here will fail if that assumption
       fails.

    Args:
        candidates (iterable): Iterable of pairs of curves that may
            intersect.
        strategy (Optional[~bezier.curve.IntersectionStrategy]): The
            intersection algorithm to use. Defaults to geometric.

    Returns:
        list: List of all :class:`Intersection`s (possibly empty).

    Raises:
        ValueError: If the strategy is not known.
    """
    if strategy is IntersectionStrategy.geometric:
        return _all_intersections_geometric(candidates)
    elif strategy is IntersectionStrategy.algebraic:
        return _all_intersections_algebraic(candidates)
    else:
        raise ValueError('Unexpected strategy.')


class BoxIntersectionType(object):  # pylint: disable=too-few-public-methods
    """Enum representing all possible bounding box intersections.

    .. note::

       This class would be more "correct" as an ``enum.Enum``, but it we keep
       the values integers to make interfacing with Fortran easier.
    """
    INTERSECTION = 0
    TANGENT = 1
    DISJOINT = 2


class Linearization(object):
    """A linearization of a curve.

    This class is provided as a stand-in for a curve, so it
    provides a similar interface.

    Args:
        curve (.Curve): A curve that is linearized.
        error (float): The linearization error. Expected to have been
            computed via :func:`.linearization_error`.
    """

    __slots__ = ('curve', 'error', 'start_node', 'end_node')

    def __init__(self, curve, error):
        self.curve = curve
        """Curve: The curve that this linearization approximates."""
        self.error = error
        """float: The linearization error for the linearized curve."""
        # NOTE: We want the nodes to be 1x2 but accessing
        #       ``curve._nodes[[0], :]`` makes a copy while the access
        #       below **does not** copy. See
        #       (https://docs.scipy.org/doc/numpy-1.6.0/reference/
        #        arrays.indexing.html#advanced-indexing)
        self.start_node = curve._nodes[0, :].reshape((1, 2), order='F')
        """numpy.ndarray: The start vector of this linearization."""
        self.end_node = curve._nodes[-1, :].reshape((1, 2), order='F')
        """numpy.ndarray: The end vector of this linearization."""

    def subdivide(self):
        """Do-nothing method to match the :class:`.Curve` interface.

        Returns:
            Tuple[~bezier._intersection_helpers.Linearization]: List of all
            subdivided parts, which is just the current object.
        """
        return self,

    @classmethod
    def from_shape(cls, shape):
        """Try to linearize a curve (or an already linearized curve).

        Args:
            shape (Union[~bezier.curve.Curve, \
            ~bezier._intersection_helpers.Linearization]): A curve or an
                already linearized curve.

        Returns:
            Union[~bezier.curve.Curve, \
            ~bezier._intersection_helpers.Linearization]: The
            (potentially linearized) curve.
        """
        # NOTE: In the below we replace ``isinstance(a, B)`` with
        #       ``a.__class__ is B``, which is a 3-3.5x speedup.
        if shape.__class__ is cls:
            return shape
        else:
            error = linearization_error(shape._nodes, shape._degree)
            if error < _ERROR_VAL:
                linearized = cls(shape, error)
                return linearized
            else:
                return shape


class Intersection(object):  # pylint: disable=too-few-public-methods
    """Representation of a curve-curve intersection.

    Args:
        first (.Curve): The "first" curve in the intersection.
        s (float): The parameter along ``first`` where the
            intersection occurs.
        second (.Curve): The "second" curve in the intersection.
        t (float): The parameter along ``second`` where the
            intersection occurs.
        point (Optional[numpy.ndarray]): The point where the two
            curves actually intersect.
        interior_curve (Optional[ \
            ~bezier._surface_helpers.IntersectionClassification]): The
            classification of the intersection.
    """

    __slots__ = ('first', 's', 'second', 't',
                 'point', 'interior_curve')

    def __init__(self, first, s, second, t, point=None, interior_curve=None):
        self.first = first
        """Curve: The "first" curve in the intersection."""
        self.s = s
        """float: The intersection parameter for the ``first`` curve."""
        self.second = second
        """Curve: The "second" curve in the intersection."""
        self.t = t
        """float: The intersection parameter for the ``second`` curve."""
        self.point = point
        """numpy.ndarray: The point where the intersection occurs."""
        self.interior_curve = interior_curve
        """IntersectionClassification: Which of the curves is on the interior.

        See :func:`.classify_intersection` for more details.
        """

    def get_point(self):
        """The point where the intersection occurs.

        This exists primarily for :meth:`.Curve.intersect`.

        Returns:
            numpy.ndarray: The point where the intersection occurs.
            Returns ``point`` if stored on the current value, otherwise
            computes the value on the fly.
        """
        if self.point is None:
            return _check_close(
                self.s, self.first, self.t, self.second)
        else:
            return self.point


# pylint: disable=invalid-name
if _speedup is None:  # pragma: NO COVER
    linearization_error = _linearization_error
    segment_intersection = _segment_intersection
    newton_refine = _newton_refine
    bbox_intersect = _bbox_intersect
    parallel_different = _parallel_different
    _from_linearized_low_level = _from_linearized_low_level_py
else:
    linearization_error = _speedup.speedup.linearization_error
    segment_intersection = _speedup.speedup.segment_intersection
    newton_refine = _speedup.speedup.newton_refine_intersect
    bbox_intersect = _speedup.speedup.bbox_intersect
    parallel_different = _speedup.speedup.parallel_different
    _from_linearized_low_level = _speedup.speedup.from_linearized
# pylint: enable=invalid-name
