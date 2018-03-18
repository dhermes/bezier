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
"""Private helper methods for intersecting B |eacute| zier shapes.

As a convention, the functions defined here with a leading underscore
(e.g. :func:`_newton_refine`) have a special meaning.

Each of these functions have a Cython speedup with the exact same
interface which calls out to a Fortran implementation. The speedup
will be used if the extension can be built. The name **without** the
leading underscore will be surfaced as the actual interface (e.g.
``newton_refine``) whether that is the pure Python implementation
or the speedup.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import enum

import numpy as np
import six

from bezier import _curve_helpers
from bezier import _helpers

try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None
# For ``full_newton()``.
ZERO_THRESHOLD = 0.5 ** 10  # ~1e-3
MAX_NEWTON_ITERATIONS = 10
NEWTON_ERROR_RATIO = 0.5 ** 36
NEWTON_NO_CONVERGE = """\
Unsupported multiplicity.

Newton's method failed to converge to a solution under the
following assumptions:

- The starting ``s-t`` values were already near a solution
- The root / solution has multiplicity 1 or 2
  - 1: The root is "simple", i.e. the curves are not tangent
       and have no self-intersections at the point of intersection.
  - 2: The root is a double root, i.e. the curves are tangent
       but have different curvatures at the point of intersection.

The failure to converge may have been caused by one of:

- The root was of multiplicity greater than 2
- The curves don't actually intersect, though they come very close
- Numerical issues caused the iteration to leave the region
  of convergence
"""


def _newton_refine(s, nodes1, t, nodes2):
    r"""Apply one step of 2D Newton's method.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

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

       machine_eps = np.finfo(np.float64).eps

       def cuberoot(value):
           return np.cbrt(value)

    .. doctest:: newton-refine1

       >>> nodes1 = np.asfortranarray([
       ...     [0.0, 2.0, 4.0],
       ...     [0.0, 4.0, 0.0],
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
       ...     [0.0, 0.25,  0.5, 0.75, 1.0],
       ...     [0.0, 2.0 , -2.0, 2.0 , 0.0],
       ... ])
       >>> nodes2 = np.asfortranarray([
       ...     [0.0, 0.25, 0.5, 0.75, 1.0],
       ...     [1.0, 0.5 , 0.5, 0.5 , 0.0],
       ... ])
       >>> # The expected intersection is the only real root of
       >>> # 28 s^3 - 30 s^2 + 9 s - 1.
       >>> omega = cuberoot(28.0 * np.sqrt(17.0) + 132.0) / 28.0
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
       >>> np.allclose(s_vals[4], expected, rtol=machine_eps, atol=0.0)
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
       ...     [0.0, 0.5, 1.0],
       ...     [0.0, 1.0, 0.0],
       ... ])
       >>> nodes2 = np.asfortranarray([
       ...     [0.0, 1.0],
       ...     [0.5, 0.5],
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
           [0.0, 0.5, 1.0],
           [0.0, 1.0, 0.0],
       ])
       nodes2 = np.asfortranarray([
           [0.0, 1.0],
           [0.5, 0.5],
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

    Raises:
        ValueError: If the Jacobian is singular at ``(s, t)``.
    """
    # NOTE: We form -F(s, t) since we want to solve -DF^{-1} F(s, t).
    func_val = (
        _curve_helpers.evaluate_multi(nodes2, np.asfortranarray([t])) -
        _curve_helpers.evaluate_multi(nodes1, np.asfortranarray([s]))
    )
    if np.all(func_val == 0.0):
        # No refinement is needed.
        return s, t

    # NOTE: This assumes the curves are 2D.
    jac_mat = np.empty((2, 2), order='F')
    jac_mat[:, :1] = _curve_helpers.evaluate_hodograph(s, nodes1)
    jac_mat[:, 1:] = -_curve_helpers.evaluate_hodograph(t, nodes2)
    # Solve the system.
    singular, delta_s, delta_t = _helpers.solve2x2(jac_mat, func_val[:, 0])
    if singular:
        raise ValueError('Jacobian is singular.')

    else:
        return s + delta_s, t + delta_t


class NewtonSimpleRoot(object):  # pylint: disable=too-few-public-methods
    r"""Callable object that facilitates Newton's method.

    This is meant to be used to compute the Newton update via:

    .. math::

       DF(s, t) \left[\begin{array}{c}
           \Delta s \\ \Delta t \end{array}\right] = -F(s, t).

    Args:
        nodes1 (numpy.ndarray): Control points of the first curve.
        first_deriv1 (numpy.ndarray): Control points of the curve
            :math:`B_1'(s)`.
        nodes2 (numpy.ndarray): Control points of the second curve.
        first_deriv2 (numpy.ndarray): Control points of the curve
            :math:`B_2'(t)`.
    """

    def __init__(self, nodes1, first_deriv1, nodes2, first_deriv2):
        self.nodes1 = nodes1
        self.first_deriv1 = first_deriv1
        self.nodes2 = nodes2
        self.first_deriv2 = first_deriv2

    def __call__(self, s, t):
        r"""This computes :math:`F = B_1(s) - B_2(t)` and :math:`DF(s, t)`.

        .. note::

           There is **almost** identical code in :func:`._newton_refine`, but
           that code can avoid computing the ``first_deriv1`` and
           ``first_deriv2`` nodes in cases that :math:`F(s, t) = 0` whereas
           this function assumes they have been given.

        In the case that :math:`DF(s, t)` is singular, the assumption is that
        the intersection has a multiplicity higher than one (i.e. the root is
        non-simple). **Near** a simple root, it must be the case that
        :math:`DF(s, t)` has non-zero determinant, so due to continuity, we
        assume the Jacobian will be invertible nearby.

        Args:
            s (float): The parameter where we'll compute :math:`B_1(s)` and
                :math:`DF(s, t)`.
            t (float): The parameter where we'll compute :math:`B_2(t)` and
                :math:`DF(s, t)`.

        Returns:
            Tuple[Optional[numpy.ndarray], numpy.ndarray]: Pair of

            * The LHS matrix ``DF``, a ``2 x 2`` array. If ``F == 0`` then
              this matrix won't be computed and :data:`None` will be returned.
            * The RHS vector ``F``, a ``2 x 1`` array.
        """
        s_vals = np.asfortranarray([s])
        b1_s = _curve_helpers.evaluate_multi(self.nodes1, s_vals)
        t_vals = np.asfortranarray([t])
        b2_t = _curve_helpers.evaluate_multi(self.nodes2, t_vals)
        func_val = b1_s - b2_t
        if np.all(func_val == 0.0):
            return None, func_val

        else:
            jacobian = np.empty((2, 2), order='F')
            jacobian[:, :1] = _curve_helpers.evaluate_multi(
                self.first_deriv1, s_vals
            )
            jacobian[:, 1:] = -_curve_helpers.evaluate_multi(
                self.first_deriv2, t_vals
            )
            return jacobian, func_val


class NewtonDoubleRoot(object):  # pylint: disable=too-few-public-methods
    r"""Callable object that facilitates Newton's method for double roots.

    This is an augmented version of :class:`NewtonSimpleRoot`.

    For non-simple intersections (i.e. multiplicity greater than 1),
    the curves will be tangent, which forces :math:`B_1'(s) \times B_2'(t)`
    to be zero. Unfortunately, that quantity is also equal to the
    determinant of the Jacobian, so :math:`DF` will not be full rank.

    In order to produce a system that **can** be solved, an
    an augmented function is computed:

    .. math::

        G(s, t) = \left[\begin{array}{c}
            F(s, t) \\ \hline
            B_1'(s) \times B_2'(t)
            \end{array}\right]

    The use of :math:`B_1'(s) \times B_2'(t)` (with lowered degree in
    :math:`s` and :math:`t`) means that the rank deficiency in
    :math:`DF` can be fixed in:

    .. math::

        DG(s, t) = \left[\begin{array}{c | c}
            B_1'(s) & -B_2'(t) \\ \hline
            B_1''(s) \times B_2'(t) & B_1'(s) \times B_2''(t)
            \end{array}\right]

    (This may not always be full rank, but in the double root / multiplicity
    2 case it will be full rank near a solution.)

    Rather than finding a least squares solution to the overdetermined system

    .. math::

       DG(s, t) \left[\begin{array}{c}
           \Delta s \\ \Delta t \end{array}\right] = -G(s, t)

    we find a solution to the square (and hopefully full rank) system:

    .. math::

       DG^T DG \left[\begin{array}{c}
           \Delta s \\ \Delta t \end{array}\right] = -DG^T G.

    Forming ``DG^T DG`` squares the condition number, so it would be "better"
    to use :func:`~numpy.linalg.lstsq` (which wraps the LAPACK routine
    ``dgelsd``). However, using :func:`.solve2x2` is **much** more
    straightforward and in practice this is just as accurate.

    Args:
        nodes1 (numpy.ndarray): Control points of the first curve.
        first_deriv1 (numpy.ndarray): Control points of the curve
            :math:`B_1'(s)`.
        second_deriv1 (numpy.ndarray): Control points of the curve
            :math:`B_1''(s)`.
        nodes2 (numpy.ndarray): Control points of the second curve.
        first_deriv2 (numpy.ndarray): Control points of the curve
            :math:`B_2'(t)`.
        second_deriv2 (numpy.ndarray): Control points of the curve
            :math:`B_2''(t)`.
    """

    def __init__(
        self,
        nodes1,
        first_deriv1,
        second_deriv1,
        nodes2,
        first_deriv2,
        second_deriv2,
    ):
        self.nodes1 = nodes1
        self.first_deriv1 = first_deriv1
        self.second_deriv1 = second_deriv1
        self.nodes2 = nodes2
        self.first_deriv2 = first_deriv2
        self.second_deriv2 = second_deriv2

    def __call__(self, s, t):
        r"""This computes :math:`DG^T G` and :math:`DG^T DG`.

        If :math:`DG^T DG` is not full rank, this means either :math:`DG`
        was not full rank or that it was, but with a relatively high condition
        number. So, in the case that :math:`DG^T DG` is singular, the
        assumption is that the intersection has a multiplicity higher than two.

        Args:
            s (float): The parameter where we'll compute :math:`G(s, t)` and
                :math:`DG(s, t)`.
            t (float): The parameter where we'll compute :math:`G(s, t)` and
                :math:`DG(s, t)`.

        Returns:
            Tuple[Optional[numpy.ndarray], Optional[numpy.ndarray]]: Pair of

            * The LHS matrix ``DG^T DG``, a ``2 x 2`` array. If ``G == 0`` then
              this matrix won't be computed and :data:`None` will be returned.
            * The RHS vector ``DG^T G``, a ``2 x 1`` array.
        """
        s_vals = np.asfortranarray([s])
        b1_s = _curve_helpers.evaluate_multi(self.nodes1, s_vals)
        b1_ds = _curve_helpers.evaluate_multi(self.first_deriv1, s_vals)
        t_vals = np.asfortranarray([t])
        b2_t = _curve_helpers.evaluate_multi(self.nodes2, t_vals)
        b2_dt = _curve_helpers.evaluate_multi(self.first_deriv2, t_vals)
        func_val = np.empty((3, 1), order='F')
        func_val[:2, :] = b1_s - b2_t
        func_val[2, :] = _helpers.cross_product(b1_ds[:, 0], b2_dt[:, 0])
        if np.all(func_val == 0.0):
            return None, func_val[:2, :]

        else:
            jacobian = np.empty((3, 2), order='F')
            jacobian[:2, :1] = b1_ds
            jacobian[:2, 1:] = -b2_dt
            if self.second_deriv1.size == 0:
                jacobian[2, 0] = 0.0
            else:
                jacobian[2, 0] = _helpers.cross_product(
                    _curve_helpers.evaluate_multi(self.second_deriv1, s_vals)[
                        :, 0
                    ],
                    b2_dt[:, 0],
                )
            if self.second_deriv2.size == 0:
                jacobian[2, 1] = 0.0
            else:
                jacobian[2, 1] = _helpers.cross_product(
                    b1_ds[:, 0],
                    _curve_helpers.evaluate_multi(self.second_deriv2, t_vals)[
                        :, 0
                    ],
                )
            modified_lhs = _helpers.matrix_product(jacobian.T, jacobian)
            modified_rhs = _helpers.matrix_product(jacobian.T, func_val)
            return modified_lhs, modified_rhs


def newton_iterate(evaluate_fn, s, t):
    r"""Perform a Newton iteration.

    In this function, we assume that :math:`s` and :math:`t` are nonzero,
    this makes convergence easier to detect since "relative error" at
    ``0.0`` is not a useful measure.

    There are several tolerance / threshold quantities used below:

    * :math:`10` (:attr:`MAX_NEWTON_ITERATIONS`) iterations will be done before
      "giving up". This is based on the assumption that we are already starting
      near a root, so quadratic convergence should terminate quickly.
    * :math:`\tau = \frac{1}{4}` is used as the boundary between linear
      and superlinear convergence. So if the current error
      :math:`\|p_{n + 1} - p_n\|` is not smaller than :math:`\tau` times
      the previous error :math:`\|p_n - p_{n - 1}\|`, then convergence
      is considered to be linear at that point.
    * :math:`\frac{2}{3}` of all iterations must be converging linearly
      for convergence to be stopped (and moved to the next regime). This
      will only be checked after 4 or more updates have occurred.
    * :math:`\tau = 2^{-42}` (:attr:`NEWTON_ERROR_RATIO`) is used to
      determine that an update is sufficiently small to stop iterating. So if
      the error :math:`\|p_{n + 1} - p_n\|` smaller than :math:`\tau` times
      size of the term being updated :math:`\|p_n\|`, then we
      exit with the "correct" answer.

    It is assumed that ``evaluate_fn`` will use a Jacobian return value of
    :data:`None` to indicate that :math:`F(s, t)` is exactly ``0.0``. We
    **assume** that if the function evaluates to exactly ``0.0``, then we are
    at a solution. It is possible however, that badly parameterized curves
    can evaluate to exactly ``0.0`` for inputs that are relatively far away
    from a solution (see issue #21).

    Args:
        evaluate_fn (Callable[Tuple[float, float], tuple]): A callable
            which takes :math:`s` and :math:`t` and produces an evaluated
            function value and the Jacobian matrix.
        s (float): The (first) parameter where the iteration will start.
        t (float): The (second) parameter where the iteration will start.

    Returns:
        Tuple[bool, float, float]: The triple of

        * Flag indicating if the iteration converged.
        * The current :math:`s` value when the iteration stopped.
        * The current :math:`t` value when the iteration stopped.
    """
    # Several quantities will be tracked throughout the iteration:
    # * norm_update_prev: ||p{n}   - p{n-1}|| = ||dp{n-1}||
    # * norm_update     : ||p{n+1} - p{n}  || = ||dp{n}  ||
    # * linear_updates  : This is a count on the number of times that
    #                     ``dp{n}`` "looks like" ``dp{n-1}`` (i.e.
    #                     is within a constant factor of it).
    norm_update_prev = None
    norm_update = None
    linear_updates = 0  # Track the number of "linear" updates.
    current_s = s
    current_t = t
    for index in six.moves.xrange(MAX_NEWTON_ITERATIONS):
        jacobian, func_val = evaluate_fn(current_s, current_t)
        if jacobian is None:
            return True, current_s, current_t

        singular, delta_s, delta_t = _helpers.solve2x2(
            jacobian, func_val[:, 0]
        )
        if singular:
            break

        norm_update_prev = norm_update
        norm_update = np.linalg.norm([delta_s, delta_t], ord=2)
        # If ||p{n} - p{n-1}|| > 0.25 ||p{n-1} - p{n-2}||, then that means
        # our convergence is acting linear at the current step.
        if index > 0 and norm_update > 0.25 * norm_update_prev:
            linear_updates += 1
        # If ``>=2/3`` of the updates have been linear, we are near a
        # non-simple root. (Make sure at least 5 updates have occurred.)
        if index >= 4 and 3 * linear_updates >= 2 * index:
            break

        # Determine the norm of the "old" solution before updating.
        norm_soln = np.linalg.norm([current_s, current_t], ord=2)
        current_s -= delta_s
        current_t -= delta_t
        if norm_update < NEWTON_ERROR_RATIO * norm_soln:
            return True, current_s, current_t

    return False, current_s, current_t


def full_newton_nonzero(s, nodes1, t, nodes2):
    r"""Perform a Newton iteration until convergence to a solution.

    This is the "implementation" for :func:`full_newton`. In this
    function, we assume that :math:`s` and :math:`t` are nonzero.

    Args:
        s (float): The parameter along the first curve where the iteration
            will start.
        nodes1 (numpy.ndarray): Control points of the first curve.
        t (float): The parameter along the second curve where the iteration
            will start.
        nodes2 (numpy.ndarray): Control points of the second curve.

    Returns:
        Tuple[float, float]: The pair of :math:`s` and :math:`t` values that
        Newton's method converged to.

    Raises:
        NotImplementedError: If Newton's method doesn't converge in either the
            multiplicity 1 or 2 cases.
    """
    # NOTE: We somewhat replicate code in ``evaluate_hodograph()``
    #       here. This is so we don't re-compute the nodes for the first
    #       (and possibly second) derivatives every time they are evaluated.
    _, num_nodes1 = np.shape(nodes1)
    first_deriv1 = (num_nodes1 - 1) * (nodes1[:, 1:] - nodes1[:, :-1])
    _, num_nodes2 = np.shape(nodes2)
    first_deriv2 = (num_nodes2 - 1) * (nodes2[:, 1:] - nodes2[:, :-1])
    evaluate_fn = NewtonSimpleRoot(nodes1, first_deriv1, nodes2, first_deriv2)
    converged, current_s, current_t = newton_iterate(evaluate_fn, s, t)
    if converged:
        return current_s, current_t

    # If Newton's method did not converge, then assume the root is not simple.
    second_deriv1 = (num_nodes1 - 2) * (
        first_deriv1[:, 1:] - first_deriv1[:, :-1]
    )
    second_deriv2 = (num_nodes2 - 2) * (
        first_deriv2[:, 1:] - first_deriv2[:, :-1]
    )
    evaluate_fn = NewtonDoubleRoot(
        nodes1,
        first_deriv1,
        second_deriv1,
        nodes2,
        first_deriv2,
        second_deriv2,
    )
    converged, current_s, current_t = newton_iterate(
        evaluate_fn, current_s, current_t
    )
    if converged:
        return current_s, current_t

    raise NotImplementedError(NEWTON_NO_CONVERGE)


def full_newton(s, nodes1, t, nodes2):
    r"""Perform a Newton iteration until convergence to a solution.

    This assumes :math:`s` and :math:`t` are sufficiently close to an
    intersection. It **does not** govern the maximum distance away
    that the solution can lie, though the subdivided intervals that contain
    :math:`s` and :math:`t` could be used.

    To avoid round-off issues near ``0.0``, this reverses the direction
    of a curve and replaces the parameter value :math:`\nu` with
    :math:`1 - \nu` whenever :math:`\nu < \tau` (here we use a threshold
    :math:`\tau` equal to :math:`2^{-10}`, i.e. ``ZERO_THRESHOLD``).

    Args:
        s (float): The parameter along the first curve where the iteration
            will start.
        nodes1 (numpy.ndarray): Control points of the first curve.
        t (float): The parameter along the second curve where the iteration
            will start.
        nodes2 (numpy.ndarray): Control points of the second curve.

    Returns:
        Tuple[float, float]: The pair of :math:`s` and :math:`t` values that
        Newton's method converged to.
    """
    if s < ZERO_THRESHOLD:
        reversed1 = np.asfortranarray(nodes1[:, ::-1])
        if t < ZERO_THRESHOLD:
            reversed2 = np.asfortranarray(nodes2[:, ::-1])
            refined_s, refined_t = full_newton_nonzero(
                1.0 - s, reversed1, 1.0 - t, reversed2
            )
            return 1.0 - refined_s, 1.0 - refined_t

        else:
            refined_s, refined_t = full_newton_nonzero(
                1.0 - s, reversed1, t, nodes2
            )
            return 1.0 - refined_s, refined_t

    else:
        if t < ZERO_THRESHOLD:
            reversed2 = np.asfortranarray(nodes2[:, ::-1])
            refined_s, refined_t = full_newton_nonzero(
                s, nodes1, 1.0 - t, reversed2
            )
            return refined_s, 1.0 - refined_t

        else:
            return full_newton_nonzero(s, nodes1, t, nodes2)


class IntersectionClassification(enum.Enum):
    """Enum classifying the "interior" curve in an intersection.

    Provided as the output values for :func:`.classify_intersection`.
    """
    FIRST = 0
    """The first curve is on the interior."""
    SECOND = 1
    """The second curve is on the interior."""
    OPPOSED = 2
    """Tangent intersection with opposed interiors."""
    TANGENT_FIRST = 3
    """Tangent intersection, first curve is on the interior."""
    TANGENT_SECOND = 4
    """Tangent intersection, second curve is on the interior."""
    IGNORED_CORNER = 5
    """Intersection at a corner, interiors don't intersect."""
    TANGENT_BOTH = 6
    """Tangent intersection, both curves are interior from some perspective."""
    COINCIDENT = 7
    """Intersection is actually an endpoint of a coincident segment."""
    COINCIDENT_UNUSED = 8
    """Unused because the edges are moving in opposite directions."""


class Intersection(object):  # pylint: disable=too-few-public-methods
    """Representation of a curve-curve intersection.

    Args:
        index_first (int): The index of the first curve within a list of
            curves. Expected to be used to index within the three edges of
            a surface.
        s (float): The parameter along the first curve where the
            intersection occurs.
        index_second (int): The index of the second curve within a list of
            curves. Expected to be used to index within the three edges of
            a surface.
        t (float): The parameter along the second curve where the
            intersection occurs.
        interior_curve (Optional[ \
            ~bezier._intersection_helpers.IntersectionClassification]): The
            classification of the intersection.
    """
    __slots__ = ('index_first', 's', 'index_second', 't', 'interior_curve')

    def __init__(self, index_first, s, index_second, t, interior_curve=None):
        self.index_first = index_first
        """int: Index of the first curve within a list of edges."""
        self.s = s
        """float: The intersection parameter for the first curve."""
        self.index_second = index_second
        """int: Index of the second curve within a list of edges."""
        self.t = t
        """float: The intersection parameter for the second curve."""
        self.interior_curve = interior_curve
        """IntersectionClassification: Which of the curves is on the interior.

        See :func:`.classify_intersection` for more details.
        """

    @property
    def __dict__(self):
        """dict: Dictionary of current intersection's property namespace.

        This is just a stand-in property for the usual ``__dict__``. This
        class defines ``__slots__`` so by default would not provide a
        ``__dict__``.

        This also means that the current object can't be modified by the
        returned dictionary.
        """
        return {
            'index_first': self.index_first,
            's': self.s,
            'index_second': self.index_second,
            't': self.t,
            'interior_curve': self.interior_curve,
        }


class IntersectionStrategy(enum.Enum):
    """Enum determining the type of intersection algorithm to use."""
    GEOMETRIC = 0
    """Geometric approach to intersection (via subdivision)."""
    ALGEBRAIC = 1
    """Algebraic approach to intersection (via implicitization)."""


# pylint: disable=invalid-name
if _speedup is None:  # pragma: NO COVER
    newton_refine = _newton_refine
else:
    newton_refine = _speedup.newton_refine_curve_intersect
# pylint: enable=invalid-name
