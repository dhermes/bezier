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

"""Private helper methods for :mod:`bezier.surface`.

As a convention, the functions defined here with a leading underscore
(e.g. :func:`_newton_refine`) have a special meaning.

Each of these functions have a Cython speedup with the exact same
interface which calls out to a Fortran implementation. The speedup
will be used if the extension can be built. The name **without** the
leading underscore will be surfaced as the actual interface (e.g.
``newton_refine``) whether that is the pure Python implementation
or the speedup.
"""

try:
    from bezier import _surface_intersection_speedup
except ImportError:  # pragma: NO COVER
    _surface_intersection_speedup = None


def newton_refine_solve(jac_both, x_val, surf_x, y_val, surf_y):
    r"""Helper for :func:`newton_refine`.

    We have a system:

    .. code-block:: rest

       [A C][ds] = [E]
       [B D][dt]   [F]

    This is not a typo, ``A->B->C->D`` matches the data in ``jac_both``.
    We solve directly rather than using a linear algebra utility:

    .. code-block:: rest

       ds = (D E - C F) / (A D - B C)
       dt = (A F - B E) / (A D - B C)

    Args:
        jac_both (numpy.ndarray): A 1x4 matrix of entries in a Jacobian.
        x_val (float): An ``x``-value we are trying to reach.
        surf_x (float): The actual ``x``-value we are currently at.
        y_val (float): An ``y``-value we are trying to reach.
        surf_y (float): The actual ``x``-value we are currently at.

    Returns:
        Tuple[float, float]: The pair of values the solve the
         linear system.
    """
    a_val, b_val, c_val, d_val = jac_both[0, :]
    #       and
    e_val = x_val - surf_x
    f_val = y_val - surf_y
    # Now solve:
    denom = a_val * d_val - b_val * c_val
    delta_s = (d_val * e_val - c_val * f_val) / denom
    delta_t = (a_val * f_val - b_val * e_val) / denom
    return delta_s, delta_t


def _newton_refine(nodes, degree, x_val, y_val, s, t):
    r"""Refine a solution to :math:`B(s, t) = p` using Newton's method.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Computes updates via

    .. math::

       \left[\begin{array}{c}
           0 \\ 0 \end{array}\right] \approx
           \left(B\left(s_{\ast}, t_{\ast}\right) -
           \left[\begin{array}{c} x \\ y \end{array}\right]\right) +
           \left[\begin{array}{c c}
               B_s\left(s_{\ast}, t_{\ast}\right) &
               B_t\left(s_{\ast}, t_{\ast}\right) \end{array}\right]
           \left[\begin{array}{c}
               \Delta s \\ \Delta t \end{array}\right]

    For example, (with weights
    :math:`\lambda_1 = 1 - s - t, \lambda_2 = s, \lambda_3 = t`)
    consider the surface

    .. math::

       B(s, t) =
           \left[\begin{array}{c} 0 \\ 0 \end{array}\right] \lambda_1^2 +
           \left[\begin{array}{c} 1 \\ 0 \end{array}\right]
               2 \lambda_1 \lambda_2 +
           \left[\begin{array}{c} 2 \\ 0 \end{array}\right] \lambda_2^2 +
           \left[\begin{array}{c} 2 \\ 1 \end{array}\right]
               2 \lambda_1 \lambda_3 +
           \left[\begin{array}{c} 2 \\ 2 \end{array}\right]
               2 \lambda_2 \lambda_1 +
           \left[\begin{array}{c} 0 \\ 2 \end{array}\right] \lambda_3^2

    and the point
    :math:`B\left(\frac{1}{4}, \frac{1}{2}\right) =
    \frac{1}{4} \left[\begin{array}{c} 5 \\ 5 \end{array}\right]`.

    Starting from the **wrong** point
    :math:`s = \frac{1}{2}, t = \frac{1}{4}`, we have

    .. math::

       \begin{align*}
       \left[\begin{array}{c} x \\ y \end{array}\right] -
          B\left(\frac{1}{2}, \frac{1}{4}\right) &= \frac{1}{4}
          \left[\begin{array}{c} -1 \\ 2 \end{array}\right] \\
       DB\left(\frac{1}{2}, \frac{1}{4}\right) &= \frac{1}{2}
           \left[\begin{array}{c c} 3 & 2 \\ 1 & 6 \end{array}\right] \\
       \Longrightarrow \left[\begin{array}{c}
           \Delta s \\ \Delta t \end{array}\right] &= \frac{1}{32}
           \left[\begin{array}{c} -10 \\ 7 \end{array}\right]
       \end{align*}

    .. image:: images/newton_refine_surface.png
       :align: center

    .. testsetup:: newton-refine-surface

       import numpy as np
       import bezier
       from bezier._surface_intersection import newton_refine

    .. doctest:: newton-refine-surface

       >>> nodes = np.asfortranarray([
       ...     [0.0, 0.0],
       ...     [1.0, 0.0],
       ...     [2.0, 0.0],
       ...     [2.0, 1.0],
       ...     [2.0, 2.0],
       ...     [0.0, 2.0],
       ... ])
       >>> surface = bezier.Surface(nodes, degree=2)
       >>> surface.is_valid
       True
       >>> (x_val, y_val), = surface.evaluate_cartesian(0.25, 0.5)
       >>> x_val, y_val
       (1.25, 1.25)
       >>> s, t = 0.5, 0.25
       >>> new_s, new_t = newton_refine(nodes, 2, x_val, y_val, s, t)
       >>> 32 * (new_s - s)
       -10.0
       >>> 32 * (new_t - t)
       7.0

    .. testcleanup:: newton-refine-surface

       import make_images
       make_images.newton_refine_surface(
           surface, x_val, y_val, s, t, new_s, new_t)

    Args:
        nodes (numpy.ndarray): Array of nodes in a surface.
        degree (int): The degree of the surface.
        x_val (float): The :math:`x`-coordinate of a point
            on the surface.
        y_val (float): The :math:`y`-coordinate of a point
            on the surface.
        s (float): Approximate :math:`s`-value to be refined.
        t (float): Approximate :math:`t`-value to be refined.

    Returns:
        Tuple[float, float]: The refined :math:`s` and :math:`t` values.
    """
    # NOTE: This circular import is temporary (will be resolved) once
    #       the surface intersection helpers are excised from
    #       ``_surface_helpers.py``.
    from bezier import _surface_helpers

    lambda1 = 1.0 - s - t
    (surf_x, surf_y), = _surface_helpers.evaluate_barycentric(
        nodes, degree, lambda1, s, t)
    if surf_x == x_val and surf_y == y_val:
        # No refinement is needed.
        return s, t

    # NOTE: This function assumes ``dimension==2`` (i.e. since ``x, y``).
    jac_nodes = _surface_helpers.jacobian_both(nodes, degree, 2)

    # The degree of the jacobian is one less.
    jac_both = _surface_helpers.evaluate_barycentric(
        jac_nodes, degree - 1, lambda1, s, t)

    # The first column of the jacobian matrix is B_s (i.e. the
    # left-most values in ``jac_both``).
    delta_s, delta_t = newton_refine_solve(
        jac_both, x_val, surf_x, y_val, surf_y)
    return s + delta_s, t + delta_t


# pylint: disable=invalid-name
if _surface_intersection_speedup is None:  # pragma: NO COVER
    newton_refine = _newton_refine
else:
    newton_refine = _surface_intersection_speedup.newton_refine
# pylint: enable=invalid-name
