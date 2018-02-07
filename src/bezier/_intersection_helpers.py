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

from bezier import _curve_helpers
try:
    from bezier import _speedup
except ImportError:  # pragma: NO COVER
    _speedup = None


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
    """
    # NOTE: We form -F(s, t) since we want to solve -DF^{-1} F(s, t).
    func_val = (
        _curve_helpers.evaluate_multi(nodes2, np.asfortranarray([t])) -
        _curve_helpers.evaluate_multi(nodes1, np.asfortranarray([s])))
    if np.all(func_val == 0.0):
        # No refinement is needed.
        return s, t

    # NOTE: This assumes the curves are 2D.
    jac_mat = np.empty((2, 2), order='F')
    jac_mat[:, 0] = _curve_helpers.evaluate_hodograph(s, nodes1)[:, 0]
    jac_mat[:, 1] = - _curve_helpers.evaluate_hodograph(t, nodes2)[:, 0]

    # Solve the system.
    result = np.linalg.solve(jac_mat, func_val)
    # Convert to row-vector and unpack (also makes assertion on shape).
    (delta_s, delta_t), = result.T
    return s + delta_s, t + delta_t


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
            ~bezier._surface_helpers.IntersectionClassification]): The
            classification of the intersection.
    """

    __slots__ = (
        'index_first', 's',
        'index_second', 't',
        'interior_curve',
    )

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
