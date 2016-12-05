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

"""Private helper methods for B |eacute| zier curves.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


import functools

import numpy as np
import six

try:
    import scipy.integrate as _scipy_int
except ImportError:  # pragma: NO COVER
    _scipy_int = None

from bezier import _helpers


def make_subdivision_matrix(degree):
    """Make the matrix used to subdivide a curve.

    Args:
        degree (int): The degree of the curve.

    Returns:
        numpy.ndarray: The matrix used to convert the
           nodes into left and right nodes.
    """
    num_rows = 2 * degree + 1
    result = np.zeros((num_rows, degree + 1))
    result[0, 0] = 1.0
    result[-1, -1] = 1.0
    for row in six.moves.xrange(1, degree + 1):
        half_prev = 0.5 * result[row - 1, :row]
        result[row, :row] = half_prev
        result[row, 1:row + 1] += half_prev
        # Populate the complement row as well.
        complement = num_rows - row - 1
        # NOTE: We "should" reverse the results when using
        #       the complement, but they are symmetric so
        #       that would be a waste.
        result[complement, -(row + 1):] = result[row, :row + 1]
    return result


def evaluate_multi(nodes, degree, s_vals):
    r"""Computes multiple points along a curve.

    Does so by computing the Bernstein basis at each value in ``s_vals``
    rather than using the de Casteljau algorithm.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        degree (int): The degree of the curve (assumed to be one less than
            the number of ``nodes``.
        s_vals (numpy.ndarray): Parameters along the curve (as a
            1D array).

    Returns:
        numpy.ndarray: The evaluated points on the curve as a two dimensional
        NumPy array, with the rows corresponding to each ``s``
        value and the columns to the dimension.
    """
    num_vals, = s_vals.shape

    lambda2 = s_vals[:, np.newaxis]
    lambda1 = 1.0 - lambda2

    weights_next = np.zeros((num_vals, degree + 1))
    weights_curr = np.zeros((num_vals, degree + 1))
    weights_curr[:, 0] = 1.0

    # Increase from degree 0 to ``degree``.
    for curr_deg in six.moves.xrange(degree):
        weights_next[:, :curr_deg + 1] = (
            lambda1 * weights_curr[:, :curr_deg + 1])
        weights_next[:, 1:curr_deg + 2] += (
            lambda2 * weights_curr[:, :curr_deg + 1])
        weights_curr, weights_next = weights_next, weights_curr

    return weights_curr.dot(nodes)


def _vec_size(nodes, degree, s_val):
    r"""Compute :math:`\|B(s)\|_2`.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        degree (int): The degree of the curve (assumed to be one less than
            the number of ``nodes``.
        s_val (float): Parameter to compute :math:`B(s)`.

    Returns:
        float: The norm of :math:`B(s)`.
    """
    result_vec = evaluate_multi(nodes, degree, np.array([s_val]))
    return np.linalg.norm(result_vec, ord=2)


def compute_length(nodes, degree):
    r"""Approximately compute the length of a curve.

    .. _QUADPACK: https://en.wikipedia.org/wiki/QUADPACK

    If ``degree`` is :math:`n`, then the Hodograph curve
    :math:`B'(s)` is degree :math:`d = n - 1`. Using this curve, we
    approximate the integral:

    .. math::

       \ell\left(B\right) =
           \int_0^1 \| B'(s) \|_2 \, ds

    using `QUADPACK`_ (via SciPy).

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        degree (int): The degree of the curve (assumed to be one less than
            the number of ``nodes``.

    Returns:
        float: The length of the curve.

    Raises:
        OSError: If SciPy is not installed.
    """
    # NOTE: We somewhat replicate code in ``evaluate_hodograph()``
    #       here. This is so we don't re-compute the nodes for the first
    #       derivative every time it is evaluated.
    first_deriv = degree * (nodes[1:, :] - nodes[:-1, :])
    if degree == 1:
        return np.linalg.norm(first_deriv, ord=2)

    if _scipy_int is None:
        raise OSError('This function requires SciPy for quadrature.')

    size_func = functools.partial(_vec_size, first_deriv, degree - 1)
    length, _ = _scipy_int.quad(size_func, 0.0, 1.0)
    return length


def elevate_nodes(nodes, degree, dimension):
    r"""Degree-elevate a B |eacute| zier curves.

    Does this by converting the current nodes :math:`v_0, \ldots, v_n`
    to new nodes :math:`w_0, \ldots, w_{n + 1}` where

    .. math::

       \begin{align*}
       w_0 &= v_0 \\
       w_j &= \frac{j}{n + 1} v_{j - 1} + \frac{n + 1 - j}{n + 1} v_j \\
       w_{n + 1} &= v_n
       \end{align*}

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        degree (int): The degree of the curve (assumed to be one less than
            the number of ``nodes``.
        dimension (int): The dimension of the curve.

    Returns:
        numpy.ndarray: The nodes of the degree-elevated curve.
    """
    new_nodes = np.zeros((degree + 2, dimension))

    multipliers = np.arange(1, degree + 1, dtype=float)[:, np.newaxis]
    denominator = degree + 1.0
    new_nodes[1:-1, :] = (
        multipliers * nodes[:-1, :] +
        (denominator - multipliers) * nodes[1:, :])
    # Hold off on division until the end, to (attempt to) avoid round-off.
    new_nodes /= denominator

    # After setting the internal nodes (which require division), set the
    # boundary nodes.
    new_nodes[0, :] = nodes[0, :]
    new_nodes[-1, :] = nodes[-1, :]

    return new_nodes


def de_casteljau_one_round(nodes, lambda1, lambda2):
    """Perform one round of de Casteljau's algorithm.

    The weights are assumed to sum to one.

    Args:
        nodes (numpy.ndarray): Control points for a curve.
        lambda1 (float): First barycentric weight on interval.
        lambda2 (float): Second barycentric weight on interval.

    Returns:
        numpy.ndarray: The nodes for a "blended" curve one degree
        lower.
    """
    return lambda1 * nodes[:-1, :] + lambda2 * nodes[1:, :]


def specialize_curve(nodes, degree, start, end):
    """Specialize a curve to a re-parameterization

    Does so by taking two points along the number line and then
    reparameterizing the curve onto the interval formed by the
    start and end points.

    .. note::

       This assumes the curve is degree 1 or greater but doesn't check.

    Args:
        nodes (numpy.ndarray): Control points for a curve.
        degree (int): The degree of the curve.
        start (float): The start point of the interval we are specializing to.
        end (float): The end point of the interval we are specializing to.

    Returns:
        numpy.ndarray: The control points for the specialized curve.
    """
    # Uses start-->0, end-->1 to represent the specialization used.
    weights = (
        (1.0 - start, start),
        (1.0 - end, end),
    )
    partial_vals = {
        (0,): de_casteljau_one_round(nodes, *weights[0]),
        (1,): de_casteljau_one_round(nodes, *weights[1]),
    }

    for _ in six.moves.xrange(degree - 1, 0, -1):
        new_partial = {}
        for key, sub_nodes in six.iteritems(partial_vals):
            # Our keys are ascending so we increment from the last value.
            for next_id in six.moves.xrange(key[-1], 1 + 1):
                new_key = key + (next_id,)
                new_partial[new_key] = de_casteljau_one_round(
                    sub_nodes, *weights[next_id])

        partial_vals = new_partial

    result = np.empty(nodes.shape)
    for index in six.moves.xrange(degree + 1):
        key = (0,) * (degree - index) + (1,) * index
        result[index, :] = partial_vals[key]
    return result


def evaluate_hodograph(nodes, degree, s):
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
    return degree * evaluate_multi(
        first_deriv, degree - 1, np.array([s])).flatten()


def get_curvature(nodes, degree, tangent_vec, s):
    r"""Compute the signed curvature of a curve at :math:`s`.

    Computed via

    .. math::

       \frac{B'(s) \times B''(s)}{\left\lVert B'(s) \right\rVert^3}

    Args:
        nodes (numpy.ndarray): The nodes of a curve.
        degree (int): The degree of the curve.
        tangent_vec (numpy.ndarray): The already computed value of
            :math:`B'(s)`
        s (float): The parameter value along the curve.

    Returns:
        float: The signed curvature.
    """
    if degree == 1:
        return 0.0

    # NOTE: We somewhat replicate code in ``evaluate_hodograph()``
    #       here. It may be worthwhile to implement ``Curve.hodograph()``
    #       and ``Curve.concavity()`` to avoid re-computing the
    #       first and second node differences.
    first_deriv = nodes[1:, :] - nodes[:-1, :]
    second_deriv = first_deriv[1:, :] - first_deriv[:-1, :]
    concavity = degree * (degree - 1) * evaluate_multi(
        second_deriv, degree - 2, np.array([s]))
    # NOTE: This assumes ``tangent_vec`` was ``flatten()``-ed, but
    #       we intentionally don't flatten ``concavity``.
    curvature = _helpers.cross_product(
        np.array([tangent_vec]), concavity)
    curvature /= np.linalg.norm(tangent_vec, ord=2)**3
    return curvature
