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

"""Pure Python implementations of helper methods for B |eacute| zier curves.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import functools

import numpy as np

from bezier.hazmat import helpers as _py_helpers


_MAX_LOCATE_SUBDIVISIONS = 20
_LOCATE_STD_CAP = 0.5**20
_FLOAT64 = np.float64  # pylint: disable=no-member
_REDUCE_THRESHOLD = 0.5**26  # sqrt(machine precision)
# Projections onto the space of degree-elevated nodes.
# If v --> vE is the (right) elevation map, then P = E^T (E E^T)^{-1} E
# is the (right) projection.
_PROJECTION0 = np.asfortranarray([[0.5, 0.5], [0.5, 0.5]])
_PROJ_DENOM0 = 1.0
_PROJECTION1 = np.asfortranarray(
    [[2.5, 1.0, -0.5], [1.0, 1.0, 1.0], [-0.5, 1.0, 2.5]]
)
_PROJ_DENOM1 = 3.0
_PROJECTION2 = np.asfortranarray(
    [
        [4.75, 0.75, -0.75, 0.25],
        [0.75, 2.75, 2.25, -0.75],
        [-0.75, 2.25, 2.75, 0.75],
        [0.25, -0.75, 0.75, 4.75],
    ]
)
_PROJ_DENOM2 = 5.0
_PROJECTION3 = np.asfortranarray(
    [
        [34.5, 2.0, -3.0, 2.0, -0.5],
        [2.0, 27.0, 12.0, -8.0, 2.0],
        [-3.0, 12.0, 17.0, 12.0, -3.0],
        [2.0, -8.0, 12.0, 27.0, 2.0],
        [-0.5, 2.0, -3.0, 2.0, 34.5],
    ]
)
_PROJ_DENOM3 = 35.0
# Reductions for a set of degree-elevated nodes.
# If v --> vE is the (right) elevation map, then R = E^T (E E^T)^{-1} -- the
# (right) pseudo-inverse of E -- actually reduces a set of nodes.
_REDUCTION0 = np.asfortranarray([[0.5], [0.5]])
_REDUCTION_DENOM0 = 1.0
_REDUCTION1 = np.asfortranarray([[2.5, -0.5], [1.0, 1.0], [-0.5, 2.5]])
_REDUCTION_DENOM1 = 3.0
_REDUCTION2 = np.asfortranarray(
    [
        [4.75, -1.25, 0.25],
        [0.75, 3.75, -0.75],
        [-0.75, 3.75, 0.75],
        [0.25, -1.25, 4.75],
    ]
)
_REDUCTION_DENOM2 = 5.0
_REDUCTION3 = np.asfortranarray(
    [
        [103.5, -26.5, 8.5, -1.5],
        [6.0, 106.0, -34.0, 6.0],
        [-9.0, 51.0, 51.0, -9.0],
        [6.0, -34.0, 106.0, 6.0],
        [-1.5, 8.5, -26.5, 103.5],
    ]
)
_REDUCTION_DENOM3 = 105.0
_LINEAR_SUBDIVIDE_LEFT = np.asfortranarray([[1.0, 0.5], [0.0, 0.5]])
_LINEAR_SUBDIVIDE_RIGHT = np.asfortranarray([[0.5, 0.0], [0.5, 1.0]])
_QUADRATIC_SUBDIVIDE_LEFT = np.asfortranarray(
    [[1.0, 0.5, 0.25], [0.0, 0.5, 0.5], [0.0, 0.0, 0.25]]
)
_QUADRATIC_SUBDIVIDE_RIGHT = np.asfortranarray(
    [[0.25, 0.0, 0.0], [0.5, 0.5, 0.0], [0.25, 0.5, 1.0]]
)
_CUBIC_SUBDIVIDE_LEFT = np.asfortranarray(
    [
        [1.0, 0.5, 0.25, 0.125],
        [0.0, 0.5, 0.5, 0.375],
        [0.0, 0.0, 0.25, 0.375],
        [0.0, 0.0, 0.0, 0.125],
    ]
)
_CUBIC_SUBDIVIDE_RIGHT = np.asfortranarray(
    [
        [0.125, 0.0, 0.0, 0.0],
        [0.375, 0.25, 0.0, 0.0],
        [0.375, 0.5, 0.5, 0.0],
        [0.125, 0.25, 0.5, 1.0],
    ]
)


def make_subdivision_matrices(degree):
    """Make the matrix used to subdivide a curve.

    .. note::

        This is a helper for :func:`subdivide_nodes`. It does not have a
        Fortran speedup because it is **only** used by a function which has
        a Fortran speedup.

    Args:
        degree (int): The degree of the curve.

    Returns:
        Tuple[numpy.ndarray, numpy.ndarray]: The matrices used to convert
        the nodes into left and right nodes, respectively.
    """
    left = np.zeros((degree + 1, degree + 1), order="F")
    right = np.zeros((degree + 1, degree + 1), order="F")
    left[0, 0] = 1.0
    right[-1, -1] = 1.0
    for col in range(1, degree + 1):
        half_prev = 0.5 * left[:col, col - 1]
        left[:col, col] = half_prev
        left[1 : col + 1, col] += half_prev  # noqa: E203
        # Populate the complement col (in right) as well.
        complement = degree - col
        # NOTE: We "should" reverse the results when using
        #       the complement, but they are symmetric so
        #       that would be a waste.
        right[-(col + 1) :, complement] = left[: col + 1, col]  # noqa: E203
    return left, right


def subdivide_nodes(nodes):
    """Subdivide a curve into two sub-curves.

    Does so by taking the unit interval (i.e. the domain of the curve) and
    splitting it into two sub-intervals by splitting down the middle.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a B |eacute| zier curve.

    Returns:
        Tuple[numpy.ndarray, numpy.ndarray]: The nodes for the two sub-curves.
    """
    _, num_nodes = np.shape(nodes)
    if num_nodes == 2:
        left_nodes = _py_helpers.matrix_product(nodes, _LINEAR_SUBDIVIDE_LEFT)
        right_nodes = _py_helpers.matrix_product(
            nodes, _LINEAR_SUBDIVIDE_RIGHT
        )
    elif num_nodes == 3:
        left_nodes = _py_helpers.matrix_product(
            nodes, _QUADRATIC_SUBDIVIDE_LEFT
        )
        right_nodes = _py_helpers.matrix_product(
            nodes, _QUADRATIC_SUBDIVIDE_RIGHT
        )
    elif num_nodes == 4:
        left_nodes = _py_helpers.matrix_product(nodes, _CUBIC_SUBDIVIDE_LEFT)
        right_nodes = _py_helpers.matrix_product(nodes, _CUBIC_SUBDIVIDE_RIGHT)
    else:
        left_mat, right_mat = make_subdivision_matrices(num_nodes - 1)
        left_nodes = _py_helpers.matrix_product(nodes, left_mat)
        right_nodes = _py_helpers.matrix_product(nodes, right_mat)
    return left_nodes, right_nodes


def evaluate_multi(nodes, s_vals):
    r"""Computes multiple points along a curve.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        s_vals (numpy.ndarray): Parameters along the curve (as a
            1D array).

    Returns:
        numpy.ndarray: The evaluated points on the curve as a two dimensional
        NumPy array, with the columns corresponding to each ``s``
        value and the rows to the dimension.
    """
    one_less = 1.0 - s_vals
    return evaluate_multi_barycentric(nodes, one_less, s_vals)


def evaluate_multi_vs(nodes, lambda1, lambda2):
    r"""Evaluates a B |eacute| zier type-function.

    .. _VS Algorithm: https://doi.org/10.1016/0167-8396(86)90018-X

    Of the form

    .. math::

       B(\lambda_1, \lambda_2) = \sum_j \binom{n}{j}
           \lambda_1^{n - j} \lambda_2^j \cdot v_j

    for some set of vectors :math:`v_j` given by ``nodes``.

    Does so via a modified Horner's method (the `VS Algorithm`_) for each
    pair of values in ``lambda1`` and ``lambda2``.

    .. math::

       \begin{align*}
       w_0 &= \lambda_1 v_0 \\
       w_j &= \lambda_1 \left[w_{j - 1} +
           \binom{n}{j} \lambda_2^j v_j\right] \\
       w_n &= w_{n - 1} + \lambda_2^n v_n \\
       B(\lambda_1, \lambda_2) &= w_n
       \end{align*}

    Additionally, binomial coefficients are computed by utilizing the fact that
    :math:`\binom{n}{j} = \binom{n}{j - 1} \frac{n - j + 1}{j}`.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        lambda1 (numpy.ndarray): Parameters along the curve (as a
            1D array).
        lambda2 (numpy.ndarray): Parameters along the curve (as a
            1D array). Typically we have ``lambda1 + lambda2 == 1``.

    Returns:
        numpy.ndarray: The evaluated points as a two dimensional
        NumPy array, with the columns corresponding to each pair of parameter
        values and the rows to the dimension.
    """
    # NOTE: We assume but don't check that lambda2 has the same shape.
    (num_vals,) = lambda1.shape
    dimension, num_nodes = nodes.shape
    degree = num_nodes - 1
    # Resize as row vectors for broadcast multiplying with
    # columns of ``nodes``.
    lambda1 = lambda1[np.newaxis, :]
    lambda2 = lambda2[np.newaxis, :]
    result = np.zeros((dimension, num_vals), order="F")
    result += lambda1 * nodes[:, [0]]
    binom_val = 1.0
    lambda2_pow = np.ones((1, num_vals), order="F")
    for index in range(1, degree):
        lambda2_pow *= lambda2
        binom_val = (binom_val * (degree - index + 1)) / index
        result += binom_val * lambda2_pow * nodes[:, [index]]
        result *= lambda1
    result += lambda2 * lambda2_pow * nodes[:, [degree]]
    return result


def evaluate_multi_de_casteljau(nodes, lambda1, lambda2):
    r"""Evaluates a B |eacute| zier type-function.

    Of the form

    .. math::

       B(\lambda_1, \lambda_2) = \sum_j \binom{n}{j}
           \lambda_1^{n - j} \lambda_2^j \cdot v_j

    for some set of vectors :math:`v_j` given by ``nodes``.

    Does so via the de Castljau algorithm:

    .. math::

       \begin{align*}
       v_j^{(n)} &= v_j \\
       v_j^{(k)} &= \lambda_1 \cdot v_j^{(k + 1)} +
           \lambda_2 \cdot v_{j + 1}^{(k + 1)} \\
       B(\lambda_1, \lambda_2) &= v_0^{(0)}
       \end{align*}

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        lambda1 (numpy.ndarray): Parameters along the curve (as a
            1D array).
        lambda2 (numpy.ndarray): Parameters along the curve (as a
            1D array). Typically we have ``lambda1 + lambda2 == 1``.

    Returns:
        numpy.ndarray: The evaluated points as a two dimensional
        NumPy array, with the columns corresponding to each pair of parameter
        values and the rows to the dimension.
    """
    # NOTE: We assume but don't check that lambda2 has the same shape.
    (num_vals,) = lambda1.shape
    dimension, num_nodes = nodes.shape
    degree = num_nodes - 1

    lambda1_wide = np.empty((dimension, num_vals, degree), order="F")
    lambda2_wide = np.empty((dimension, num_vals, degree), order="F")
    workspace = np.empty((dimension, num_vals, degree), order="F")
    for index in range(num_vals):
        lambda1_wide[:, index, :] = lambda1[index]
        lambda2_wide[:, index, :] = lambda2[index]
        workspace[:, index, :] = (
            lambda1[index] * nodes[:, :degree] + lambda2[index] * nodes[:, 1:]
        )

    for index in range(degree - 1, 0, -1):
        workspace[:, :, :index] = (
            lambda1_wide[:, :, :index] * workspace[:, :, :index]
            + lambda2_wide[:, :, :index]
            * workspace[:, :, 1 : (index + 1)]  # noqa: E203
        )

    # NOTE: This returns an array with `evaluated.flags.owndata` false, though
    #       it is Fortran contiguous.
    return workspace[:, :, 0]


def evaluate_multi_barycentric(nodes, lambda1, lambda2):
    r"""Evaluates a B |eacute| zier type-function.

    Of the form

    .. math::

       B(\lambda_1, \lambda_2) = \sum_j \binom{n}{j}
           \lambda_1^{n - j} \lambda_2^j \cdot v_j

    for some set of vectors :math:`v_j` given by ``nodes``. This uses the
    more efficient :func:`.evaluate_multi_vs` until degree 55, at which point
    :math:`\binom{55}{26}` and other coefficients cannot be computed exactly.
    For degree 55 and higher, the classical de Casteljau algorithm will be
    used via :func:`.evaluate_multi_de_casteljau`.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        lambda1 (numpy.ndarray): Parameters along the curve (as a
            1D array).
        lambda2 (numpy.ndarray): Parameters along the curve (as a
            1D array). Typically we have ``lambda1 + lambda2 == 1``.

    Returns:
        numpy.ndarray: The evaluated points as a two dimensional
        NumPy array, with the columns corresponding to each pair of parameter
        values and the rows to the dimension.
    """
    _, num_nodes = nodes.shape
    # NOTE: The computation of (degree C k) values in ``evaluate_multi_vs``
    #       starts to introduce round-off when computing (55 C 26). For very
    #       large degree, we ditch the VS algorithm and use de Casteljau
    #       (which has quadratic runtime and cubic space usage).
    if num_nodes > 55:
        return evaluate_multi_de_casteljau(nodes, lambda1, lambda2)

    return evaluate_multi_vs(nodes, lambda1, lambda2)


def vec_size(nodes, s_val):
    r"""Compute :math:`\|B(s)\|_2`.

    .. note::

        This is a helper for :func:`compute_length` and does not have
        a Fortran speedup.

    Intended to be used with :func:`functools.partial` to fill in the
    value of ``nodes`` and create a callable that only accepts ``s_val``.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        s_val (float): Parameter to compute :math:`B(s)`.

    Returns:
        float: The norm of :math:`B(s)`.
    """
    result_vec = evaluate_multi(nodes, np.asfortranarray([s_val]))
    # NOTE: We convert to 1D to make sure NumPy uses vector norm.
    return np.linalg.norm(result_vec[:, 0], ord=2)


def compute_length(nodes):
    r"""Approximately compute the length of a curve.

    .. _QUADPACK: https://en.wikipedia.org/wiki/QUADPACK

    If ``degree`` is :math:`n`, then the Hodograph curve
    :math:`B'(s)` is degree :math:`d = n - 1`. Using this curve, we
    approximate the integral:

    .. math::

       \int_{B\left(\left[0, 1\right]\right)} 1 \, d\mathbf{x} =
       \int_0^1 \left\lVert B'(s) \right\rVert_2 \, ds

    using `QUADPACK`_ (via SciPy).

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.

    Returns:
        float: The length of the curve.

    Raises:
        ValueError: If ``nodes`` has zero columns.
    """
    _, num_nodes = np.shape(nodes)
    # NOTE: We somewhat replicate code in ``evaluate_hodograph()``
    #       here. This is so we don't re-compute the nodes for the first
    #       derivative every time it is evaluated.
    first_deriv = (num_nodes - 1) * (nodes[:, 1:] - nodes[:, :-1])
    if num_nodes == 0:
        raise ValueError("Curve should have at least one node.")

    if num_nodes == 1:
        return 0.0

    if num_nodes == 2:
        # NOTE: We convert to 1D to make sure NumPy uses vector norm.
        return np.linalg.norm(first_deriv[:, 0], ord=2)

    # NOTE: We import SciPy at runtime to avoid the import-time cost for users
    #       that don't need pure Python curve helpers (e.g. if the ``_speedup``
    #       module is available). The ``scipy`` import is a tad expensive.
    import scipy.integrate  # pylint: disable=import-outside-toplevel

    size_func = functools.partial(vec_size, first_deriv)
    length, _ = scipy.integrate.quad(size_func, 0.0, 1.0)
    return length


def elevate_nodes(nodes):
    r"""Degree-elevate a B |eacute| zier curve.

    Does this by converting the current nodes :math:`v_0, \ldots, v_n`
    to new nodes :math:`w_0, \ldots, w_{n + 1}` where

    .. math::

       \begin{align*}
       w_0 &= v_0 \\
       w_j &= \frac{j}{n + 1} v_{j - 1} + \frac{n + 1 - j}{n + 1} v_j \\
       w_{n + 1} &= v_n
       \end{align*}

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.

    Returns:
        numpy.ndarray: The nodes of the degree-elevated curve.
    """
    dimension, num_nodes = np.shape(nodes)
    new_nodes = np.empty((dimension, num_nodes + 1), order="F")
    multipliers = np.arange(1, num_nodes, dtype=_FLOAT64)[np.newaxis, :]
    denominator = float(num_nodes)
    new_nodes[:, 1:-1] = (
        multipliers * nodes[:, :-1]
        + (denominator - multipliers) * nodes[:, 1:]
    )
    # Hold off on division until the end, to (attempt to) avoid round-off.
    new_nodes /= denominator
    # After setting the internal nodes (which require division), set the
    # boundary nodes.
    new_nodes[:, 0] = nodes[:, 0]
    new_nodes[:, -1] = nodes[:, -1]
    return new_nodes


def de_casteljau_one_round(nodes, lambda1, lambda2):
    """Perform one round of de Casteljau's algorithm.

    .. note::

        This is a helper for :func:`specialize_curve`. It does not have a
        Fortran speedup because it is **only** used by a function which has
        a Fortran speedup.

    The weights are assumed to sum to one.

    Args:
        nodes (numpy.ndarray): Control points for a curve.
        lambda1 (float): First barycentric weight on interval.
        lambda2 (float): Second barycentric weight on interval.

    Returns:
        numpy.ndarray: The nodes for a "blended" curve one degree
        lower.
    """
    return np.asfortranarray(lambda1 * nodes[:, :-1] + lambda2 * nodes[:, 1:])


def specialize_curve(nodes, start, end):
    """Specialize a curve to a re-parameterization

    .. note::

       This assumes the curve is degree 1 or greater but doesn't check.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): Control points for a curve.
        start (float): The start point of the interval we are specializing to.
        end (float): The end point of the interval we are specializing to.

    Returns:
        numpy.ndarray: The control points for the specialized curve.
    """
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-locals
    _, num_nodes = np.shape(nodes)
    # Uses start-->0, end-->1 to represent the specialization used.
    weights = ((1.0 - start, start), (1.0 - end, end))
    partial_vals = {
        (0,): de_casteljau_one_round(nodes, *weights[0]),
        (1,): de_casteljau_one_round(nodes, *weights[1]),
    }
    for _ in range(num_nodes - 2, 0, -1):
        new_partial = {}
        for key, sub_nodes in partial_vals.items():
            # Our keys are ascending so we increment from the last value.
            for next_id in range(key[-1], 1 + 1):
                new_key = key + (next_id,)
                new_partial[new_key] = de_casteljau_one_round(
                    sub_nodes, *weights[next_id]
                )
        partial_vals = new_partial
    result = np.empty(nodes.shape, order="F")
    for index in range(num_nodes):
        key = (0,) * (num_nodes - index - 1) + (1,) * index
        result[:, [index]] = partial_vals[key]
    return result


def evaluate_hodograph(s, nodes):
    r"""Evaluate the Hodograph curve at a point :math:`s`.

    The Hodograph (first derivative) of a B |eacute| zier curve
    degree :math:`d = n - 1` and is given by

    .. math::

       B'(s) = n \sum_{j = 0}^{d} \binom{d}{j} s^j
       (1 - s)^{d - j} \cdot \Delta v_j

    where each forward difference is given by
    :math:`\Delta v_j = v_{j + 1} - v_j`.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes of a curve.
        s (float): A parameter along the curve at which the Hodograph
            is to be evaluated.

    Returns:
        numpy.ndarray: The point on the Hodograph curve (as a two
        dimensional NumPy array with a single row).
    """
    _, num_nodes = np.shape(nodes)
    first_deriv = nodes[:, 1:] - nodes[:, :-1]
    return (num_nodes - 1) * evaluate_multi(
        first_deriv, np.asfortranarray([s])
    )


def get_curvature(nodes, tangent_vec, s):
    r"""Compute the signed curvature of a curve at :math:`s`.

    Computed via

    .. math::

       \frac{B'(s) \times B''(s)}{\left\lVert B'(s) \right\rVert_2^3}

    .. image:: ../../images/get_curvature.png
       :align: center

    .. testsetup:: get-curvature

       import numpy as np
       import bezier
       from bezier.hazmat.curve_helpers import evaluate_hodograph
       from bezier.hazmat.curve_helpers import get_curvature

    .. doctest:: get-curvature
       :options: +NORMALIZE_WHITESPACE

       >>> import numpy as np
       >>> nodes = np.asfortranarray([
       ...     [1.0, 0.75,  0.5, 0.25, 0.0],
       ...     [0.0, 2.0 , -2.0, 2.0 , 0.0],
       ... ])
       >>> s = 0.5
       >>> tangent_vec = evaluate_hodograph(s, nodes)
       >>> tangent_vec
       array([[-1.],
              [ 0.]])
       >>> curvature = get_curvature(nodes, tangent_vec, s)
       >>> curvature
       -12.0

    .. testcleanup:: get-curvature

       import make_images
       make_images.get_curvature(nodes, s, tangent_vec, curvature)

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes of a curve.
        tangent_vec (numpy.ndarray): The already computed value of
            :math:`B'(s)`
        s (float): The parameter value along the curve.

    Returns:
        float: The signed curvature.
    """
    _, num_nodes = np.shape(nodes)
    if num_nodes == 2:  # Lines have no curvature.
        return 0.0

    # NOTE: We somewhat replicate code in ``evaluate_hodograph()`` here.
    first_deriv = nodes[:, 1:] - nodes[:, :-1]
    second_deriv = first_deriv[:, 1:] - first_deriv[:, :-1]
    concavity = (
        (num_nodes - 1)
        * (num_nodes - 2)
        * evaluate_multi(second_deriv, np.asfortranarray([s]))
    )
    curvature = _py_helpers.cross_product(
        tangent_vec.ravel(order="F"), concavity.ravel(order="F")
    )
    # NOTE: We convert to 1D to make sure NumPy uses vector norm.
    curvature /= np.linalg.norm(tangent_vec[:, 0], ord=2) ** 3
    return curvature


def newton_refine(nodes, point, s):
    r"""Refine a solution to :math:`B(s) = p` using Newton's method.

    Computes updates via

    .. math::

       \mathbf{0} \approx
           \left(B\left(s_{\ast}\right) - p\right) +
           B'\left(s_{\ast}\right) \Delta s

    For example, consider the curve

    .. math::

       B(s) =
           \left[\begin{array}{c} 0 \\ 0 \end{array}\right] (1 - s)^2 +
           \left[\begin{array}{c} 1 \\ 2 \end{array}\right]
               2 (1 - s) s +
           \left[\begin{array}{c} 3 \\ 1 \end{array}\right] s^2

    and the point :math:`B\left(\frac{1}{4}\right) =
    \frac{1}{16} \left[\begin{array}{c} 9 \\ 13 \end{array}\right]`.

    Starting from the **wrong** point :math:`s = \frac{3}{4}`, we have

    .. math::

       \begin{align*}
       p - B\left(\frac{1}{2}\right) &= -\frac{1}{2}
           \left[\begin{array}{c} 3 \\ 1 \end{array}\right] \\
       B'\left(\frac{1}{2}\right) &= \frac{1}{2}
           \left[\begin{array}{c} 7 \\ -1 \end{array}\right] \\
       \Longrightarrow \frac{1}{4} \left[\begin{array}{c c}
           7 & -1 \end{array}\right] \left[\begin{array}{c}
           7 \\ -1 \end{array}\right] \Delta s &= -\frac{1}{4}
           \left[\begin{array}{c c} 7 & -1 \end{array}\right]
           \left[\begin{array}{c} 3 \\ 1 \end{array}\right] \\
       \Longrightarrow \Delta s &= -\frac{2}{5}
       \end{align*}

    .. image:: ../../images/newton_refine_curve.png
       :align: center

    .. testsetup:: newton-refine-curve, newton-refine-curve-cusp

       import numpy as np
       import bezier
       from bezier.hazmat.curve_helpers import newton_refine

    .. doctest:: newton-refine-curve
       :options: +NORMALIZE_WHITESPACE

       >>> import bezier
       >>> nodes = np.asfortranarray([
       ...     [0.0, 1.0, 3.0],
       ...     [0.0, 2.0, 1.0],
       ... ])
       >>> curve = bezier.Curve(nodes, degree=2)
       >>> point = curve.evaluate(0.25)
       >>> point
       array([[0.5625],
              [0.8125]])
       >>> s = 0.75
       >>> new_s = newton_refine(nodes, point, s)
       >>> 5 * (new_s - s)
       -2.0

    .. testcleanup:: newton-refine-curve

       import make_images
       make_images.newton_refine_curve(curve, point, s, new_s)

    On curves that are not "valid" (i.e. :math:`B(s)` is not
    injective with non-zero gradient), Newton's method may
    break down and converge linearly:

    .. image:: ../../images/newton_refine_curve_cusp.png
       :align: center

    .. doctest:: newton-refine-curve-cusp
       :options: +NORMALIZE_WHITESPACE

       >>> nodes = np.asfortranarray([
       ...     [ 6.0, -2.0, -2.0, 6.0],
       ...     [-3.0,  3.0, -3.0, 3.0],
       ... ])
       >>> curve = bezier.Curve(nodes, degree=3)
       >>> expected = 0.5
       >>> point = curve.evaluate(expected)
       >>> point
       array([[0.],
              [0.]])
       >>> s_vals = [0.625, None, None, None, None, None]
       >>> np.log2(abs(expected - s_vals[0]))
       -3.0
       >>> s_vals[1] = newton_refine(nodes, point, s_vals[0])
       >>> np.log2(abs(expected - s_vals[1]))
       -3.983...
       >>> s_vals[2] = newton_refine(nodes, point, s_vals[1])
       >>> np.log2(abs(expected - s_vals[2]))
       -4.979...
       >>> s_vals[3] = newton_refine(nodes, point, s_vals[2])
       >>> np.log2(abs(expected - s_vals[3]))
       -5.978...
       >>> s_vals[4] = newton_refine(nodes, point, s_vals[3])
       >>> np.log2(abs(expected - s_vals[4]))
       -6.978...
       >>> s_vals[5] = newton_refine(nodes, point, s_vals[4])
       >>> np.log2(abs(expected - s_vals[5]))
       -7.978...

    .. testcleanup:: newton-refine-curve-cusp

       import make_images
       make_images.newton_refine_curve_cusp(curve, s_vals)

    Due to round-off, the Newton process terminates with an error that is not
    close to machine precision :math:`\varepsilon` when :math:`\Delta s = 0`.

    .. testsetup:: newton-refine-curve-cusp-continued

       import numpy as np
       import bezier
       from bezier.hazmat.curve_helpers import newton_refine

       nodes = np.asfortranarray([
           [ 6.0, -2.0, -2.0, 6.0],
           [-3.0,  3.0, -3.0, 3.0],
       ])
       curve = bezier.Curve(nodes, degree=3)
       point = curve.evaluate(0.5)

    .. doctest:: newton-refine-curve-cusp-continued

       >>> s_vals = [0.625]
       >>> new_s = newton_refine(nodes, point, s_vals[-1])
       >>> while new_s not in s_vals:
       ...     s_vals.append(new_s)
       ...     new_s = newton_refine(nodes, point, s_vals[-1])
       ...
       >>> terminal_s = s_vals[-1]
       >>> terminal_s == newton_refine(nodes, point, terminal_s)
       True
       >>> 2.0**(-31) <= abs(terminal_s - 0.5) <= 2.0**(-28)
       True

    Due to round-off near the cusp, the final error resembles
    :math:`\sqrt{\varepsilon}` rather than machine precision
    as expected.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a B |eacute| zier curve.
        point (numpy.ndarray): A point on the curve.
        s (float): An "almost" solution to :math:`B(s) = p`.

    Returns:
        float: The updated value :math:`s + \Delta s`.
    """
    pt_delta = point - evaluate_multi(nodes, np.asfortranarray([s]))
    derivative = evaluate_hodograph(s, nodes)
    # Each array is 2 x 1 (i.e. a column vector), we want the vector
    # dot product.
    delta_s = np.vdot(pt_delta[:, 0], derivative[:, 0]) / np.vdot(
        derivative[:, 0], derivative[:, 0]
    )
    return s + delta_s


def locate_point(nodes, point):
    r"""Locate a point on a curve.

    Does so by recursively subdividing the curve and rejecting
    sub-curves with bounding boxes that don't contain the point.
    After the sub-curves are sufficiently small, uses Newton's
    method to zoom in on the parameter value.

    .. note::

       This assumes, but does not check, that ``point`` is ``D x 1``,
       where ``D`` is the dimension that ``curve`` is in.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes defining a B |eacute| zier curve.
        point (numpy.ndarray): The point to locate.

    Returns:
        Optional[float]: The parameter value (:math:`s`) corresponding
        to ``point`` or :data:`None` if the point is not on the ``curve``.

    Raises:
        ValueError: If the standard deviation of the remaining start / end
            parameters among the subdivided intervals exceeds a given
            threshold (e.g. :math:`2^{-20}`).
    """
    candidates = [(0.0, 1.0, nodes)]
    for _ in range(_MAX_LOCATE_SUBDIVISIONS + 1):
        next_candidates = []
        for start, end, candidate in candidates:
            if _py_helpers.contains_nd(candidate, point.ravel(order="F")):
                midpoint = 0.5 * (start + end)
                left, right = subdivide_nodes(candidate)
                next_candidates.extend(
                    ((start, midpoint, left), (midpoint, end, right))
                )
        candidates = next_candidates
    if not candidates:
        return None

    params = [(start, end) for start, end, _ in candidates]
    if np.std(params) > _LOCATE_STD_CAP:
        raise ValueError("Parameters not close enough to one another", params)

    s_approx = np.mean(params)
    s_approx = newton_refine(nodes, point, s_approx)
    # NOTE: Since ``np.mean(params)`` must be in ``[0, 1]`` it's
    #       "safe" to push the Newton-refined value back into the unit
    #       interval.
    if s_approx < 0.0:
        return 0.0

    elif s_approx > 1.0:
        return 1.0

    else:
        return s_approx


def reduce_pseudo_inverse(nodes):
    """Performs degree-reduction for a B |eacute| zier curve.

    Does so by using the pseudo-inverse of the degree elevation
    operator (which is overdetermined).

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes in the curve.

    Returns:
        numpy.ndarray: The reduced nodes.

    Raises:
        UnsupportedDegree: If the degree is not 1, 2, 3 or 4.
    """
    _, num_nodes = np.shape(nodes)
    if num_nodes == 2:
        reduction = _REDUCTION0
        denom = _REDUCTION_DENOM0
    elif num_nodes == 3:
        reduction = _REDUCTION1
        denom = _REDUCTION_DENOM1
    elif num_nodes == 4:
        reduction = _REDUCTION2
        denom = _REDUCTION_DENOM2
    elif num_nodes == 5:
        reduction = _REDUCTION3
        denom = _REDUCTION_DENOM3
    else:
        raise _py_helpers.UnsupportedDegree(
            num_nodes - 1, supported=(1, 2, 3, 4)
        )

    result = _py_helpers.matrix_product(nodes, reduction)
    result /= denom
    return result


def projection_error(nodes, projected):
    """Compute the error between ``nodes`` and the projected nodes.

    .. note::

        This is a helper for :func:`maybe_reduce`, which is in turn a helper
        for :func:`.full_reduce`. Hence there is no corresponding Fortran
        speedup.

    For now, just compute the relative error in the Frobenius norm. But,
    we may wish to consider the error per row / point instead.

    Args:
        nodes (numpy.ndarray): Nodes in a curve.
        projected (numpy.ndarray): The ``nodes`` projected into the
            space of degree-elevated nodes.

    Returns:
        float: The relative error.
    """
    relative_err = np.linalg.norm(nodes - projected, ord="fro")
    if relative_err != 0.0:
        relative_err /= np.linalg.norm(nodes, ord="fro")
    return relative_err


def maybe_reduce(nodes):
    r"""Reduce nodes in a curve if they are degree-elevated.

    .. note::

        This is a helper for :func:`.full_reduce`. Hence there is no
        corresponding Fortran speedup.

    We check if the nodes are degree-elevated by projecting onto the
    space of degree-elevated curves of the same degree, then comparing
    to the projection. We form the projection by taking the corresponding
    (right) elevation matrix :math:`E` (from one degree lower) and forming
    :math:`E^T \left(E E^T\right)^{-1} E`.

    Args:
        nodes (numpy.ndarray): The nodes in the curve.

    Returns:
        Tuple[bool, numpy.ndarray]: Pair of values. The first indicates
        if the ``nodes`` were reduced. The second is the resulting nodes,
        either the reduced ones or the original passed in.

    Raises:
        UnsupportedDegree: If the curve is degree 5 or higher.
    """
    _, num_nodes = nodes.shape
    if num_nodes < 2:
        return False, nodes

    elif num_nodes == 2:
        projection = _PROJECTION0
        denom = _PROJ_DENOM0
    elif num_nodes == 3:
        projection = _PROJECTION1
        denom = _PROJ_DENOM1
    elif num_nodes == 4:
        projection = _PROJECTION2
        denom = _PROJ_DENOM2
    elif num_nodes == 5:
        projection = _PROJECTION3
        denom = _PROJ_DENOM3
    else:
        raise _py_helpers.UnsupportedDegree(
            num_nodes - 1, supported=(0, 1, 2, 3, 4)
        )

    projected = _py_helpers.matrix_product(nodes, projection) / denom
    relative_err = projection_error(nodes, projected)
    if relative_err < _REDUCE_THRESHOLD:
        return True, reduce_pseudo_inverse(nodes)

    else:
        return False, nodes


def full_reduce(nodes):
    """Apply degree reduction to ``nodes`` until it can no longer be reduced.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        nodes (numpy.ndarray): The nodes in the curve.

    Returns:
        numpy.ndarray: The fully degree-reduced nodes.
    """
    was_reduced, nodes = maybe_reduce(nodes)
    while was_reduced:
        was_reduced, nodes = maybe_reduce(nodes)
    return nodes


def discrete_turning_angle(nodes):
    r"""Determine the absolute sum of B |eacute| zier node angles.

    .. note::

       This assumes, but does not check, that ``nodes`` is ``2 x N``.

    For the set of vectors :math:`v_j` given by ``nodes``, the discrete
    angles :math:`\theta_j` at each internal node is given by

    .. math::

       \left(v_{j} - v_{j - 1}\right) \cdot \left(v_{j + 1} - v_{j}\right) =
         \| v_{j} - v_{j - 1} \|_2 \| v_{j + 1} - v_{j} \|_2 \cos \theta_j

    and the discrete turning angle is :math:`\sum_{j} \left|\theta_j\right|`.
    This approximates the exact turning angle

    .. math::

       \int_0^1 \left|\theta'(s)\right| \, ds.

    This is done by considering how the angle of
    :math:`B'(s) = \left[x'(s), y'(s)\right]^T` changes in small intervals
    in the parameter space; where the angle is
    :math:`\theta(s) = \arctan(y'(s) / x'(s))`. Computing
    :math:`\theta(s + ds) - \theta(s)` as :math:`ds \longrightarrow 0` leaves
    us with the **signed** angle change

    .. math::

       \theta'(s) = \frac{y''(s) x'(s) - x''(s) y'(s)}{x'(s)^2 + y'(s)^2}.

    Args:
        nodes (numpy.ndarray): The nodes in the curve.

    Returns:
        float: The (discrete) turning angle.
    """
    _, num_nodes = nodes.shape
    if num_nodes < 3:
        return 0.0

    directed = nodes[:, 1:] - nodes[:, : (num_nodes - 1)]
    vector_theta = np.arctan2(directed[1, :], directed[0, :])
    # Two values in [-pi, pi] have a difference in [-2pi, 2pi]
    angle_theta = vector_theta[1:] - vector_theta[: (num_nodes - 2)]

    result = 0.0
    for angle in angle_theta:
        if angle > np.pi:
            # Convert value in [pi, 2pi] back into [-pi, pi] and take
            # absolute value
            result += 2 * np.pi - angle
        elif angle < -np.pi:
            # Convert value in [-2pi, -pi] back into [-pi, pi] and take
            # absolute value
            result += 2 * np.pi + angle
        else:
            # Absolute value of value in [-pi, pi]
            result += abs(angle)

    return result
