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

r"""Helper for implicitizing B |eacute| zier curves.

.. _resultant: https://en.wikipedia.org/wiki/Resultant
.. _algebraic curve: https://en.wikipedia.org/wiki/Algebraic_curve
.. _Farouki and Rajan: http://dx.doi.org/10.1016/0167-8396(88)90016-7

Primarily uses the `resultant`_ to evaluate the implicitized
`algebraic curve`_. In order to do this on B |eacute| zier curves
without translating to a power basis, we utilize the work of
`Farouki and Rajan`_ to compute a modified Sylvester determinant.

Given two parametric curves :math:`(x_1(s), y_1(s))` and
:math:`(x_2(t), y_2(t))`, we can determine an "intersection polynomial"
for both :math:`s` and :math:`t`. For example, by implicitizing the
first curve, we determine :math:`f_1(x, y)` and plugging the second
curve into this we find

.. math::

   g(t) = f_1\left(x_2(t), y_2(t)\right) = 0

is the "intersection polynomial" for :math:`t`.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


import numpy as np

from bezier import _curve_helpers


def _evaluate3(nodes, x_val, y_val):
    """Helper for :func:`evaluate` when ``nodes`` is degree 3.

    Args:
        nodes (numpy.ndarray): ``4x2`` array of nodes in a curve.
        x_val (float): ``x``-coordinate for evaluation.
        y_val (float): ``y``-coordinate for evaluation.

    Returns:
        float: The computed value of :math:`f(x, y)`.
    """
    # NOTE: This may be (a) slower and (b) less precise than
    #       hard-coding the determinant.
    sylvester_mat = np.zeros((6, 6), order='F')
    delta = nodes - np.asfortranarray([[x_val, y_val]])
    delta[1:3, :] *= 3.0
    # Swap rows/columns so that x-y are right next to each other.
    # This will only change the determinant up to a sign.
    sylvester_mat[:4, :2] = delta
    sylvester_mat[1:5, 2:4] = delta
    sylvester_mat[2:, 4:] = delta
    return np.linalg.det(sylvester_mat)


def evaluate(nodes, x_val, y_val):
    r"""Evaluate the implicitized bivariate polynomial containing the curve.

    Assumes `algebraic curve`_ containing :math:`B(s, t)` is given by
    :math:`f(x, y) = 0`. This function evaluates :math:`f(x, y)`.

    .. note::

       This assumes, but doesn't check, that ``nodes`` has 2 columns.

    .. note::

       This assumes, but doesn't check, that ``nodes`` is not degree-elevated.
       If it were degree elevated, then the Sylvester matrix will always
       have zero determinant.

    Args:
        nodes (numpy.ndarray): ``Nx2`` array of nodes in a curve.
        x_val (float): ``x``-coordinate for evaluation.
        y_val (float): ``y``-coordinate for evaluation.

    Returns:
        float: The computed value of :math:`f(x, y)`.

    Raises:
        ValueError: If the curve is a point.
        NotImplementedError: If the curve is not degree 1 or 2.
    """
    num_nodes, _ = nodes.shape
    if num_nodes == 1:
        raise ValueError('A point cannot be implicitized')
    elif num_nodes == 2:
        # x(s) - x = (x0 - x) (1 - s) + (x1 - x) s
        # y(s) - y = (y0 - y) (1 - s) + (y1 - y) s
        # Modified Sylvester: [x0 - x, x1 - x]
        #                     [y0 - y, y1 - y]
        return (
            (nodes[0, 0] - x_val) * (nodes[1, 1] - y_val) -
            (nodes[1, 0] - x_val) * (nodes[0, 1] - y_val))
    elif num_nodes == 3:
        # x(s) - x = (x0 - x) (1 - s)^2 + 2 (x1 - x) s(1 - s) + (x2 - x) s^2
        # y(s) - y = (y0 - y) (1 - s)^2 + 2 (y1 - y) s(1 - s) + (y2 - y) s^2
        # Modified Sylvester: [x0 - x, 2(x1 - x),    x2 - x,      0] = A|B|C|0
        #                     [     0,    x0 - x, 2(x1 - x), x2 - x]   0|A|B|C
        #                     [y0 - y, 2(y1 - y),    y2 - y,      0]   D|E|F|0
        #                     [     0,    y0 - y, 2(y1 - y), y2 - y]   0|D|E|F
        val_a, val_b, val_c = nodes[:, 0] - x_val
        val_b *= 2
        val_d, val_e, val_f = nodes[:, 1] - y_val
        val_e *= 2
        #     [A, B, C]         [E, F, 0]
        # det [E, F, 0] = - det [A, B, C] = -E (BF - CE) + F(AF - CD)
        #     [D, E, F]         [D, E, F]
        sub1 = val_b * val_f - val_c * val_e
        sub2 = val_a * val_f - val_c * val_d
        sub_det_a = -val_e * sub1 + val_f * sub2
        #     [B, C, 0]
        # det [A, B, C] = B (BF - CE) - C (AF - CD)
        #     [D, E, F]
        sub_det_d = val_b * sub1 - val_c * sub2
        return val_a * sub_det_a + val_d * sub_det_d
    elif num_nodes == 4:
        return _evaluate3(nodes, x_val, y_val)
    else:
        raise NotImplementedError('Only degrees 1 and 2 supported')


def eval_intersection_polynomial(nodes1, nodes2, t):
    r"""Evaluates a parametric curve **on** an implicitized algebraic curve.

    Uses :func:`evaluate` to evaluate :math:`f_1(x, y)`, the implicitization
    of ``nodes1``. Then plugs ``t`` into the second parametric curve to
    get an ``x``- and ``y``-coordinate and evaluate the
    **intersection polynomial**:

    .. math::

       g(t) = f_1\left(x_2(t), y_2(t)right)

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.
        t (float): The parameter along ``nodes2`` where we evaluate
            the function.

    Returns:
        float: The computed value of :math:`f_1(x_2(t), y_2(t))`.
    """
    (x_val, y_val), = _curve_helpers.evaluate_multi(
        nodes2, np.asfortranarray([t]))
    return evaluate(nodes1, x_val, y_val)


def _to_power_basis11(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    .. _theorem: https://en.wikipedia.org/wiki/B%C3%A9zout's_theorem

    Helper for :func:`to_power_basis` in the case that each curve is
    degree one. In this case, B |eacute| zout's `theorem`_ tells us
    that the **intersection polynomial** is degree :math:`1 \cdot 1`
    hence we return two coefficients.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``2``-array of coefficients.
    """
    # We manually invert the Vandermonde matrix:
    # [1 0.0][c0] = [n0]
    # [1 1.0][c1]   [n1]
    n0 = eval_intersection_polynomial(nodes1, nodes2, 0.0)
    n1 = eval_intersection_polynomial(nodes1, nodes2, 1.0)
    # [c0] = [ 1 0][n0]
    # [c1] = [-1 1][n1]
    return np.array([n0, -n0 + n1])


def to_power_basis(nodes1, nodes2):
    """Compute the coefficients of an **intersection polynomial**.

    .. note::

       This assumes that the degree of the curve given by ``nodes1`` is
       less than or equal to the degree of that given by ``nodes2``.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``2``-array of coefficients.

    Raises:
        NotImplementedError: If the degree pair is not ``1-1``.
    """
    num_nodes1, _ = nodes1.shape
    num_nodes2, _ = nodes2.shape
    if num_nodes1 == 2:
        if num_nodes2 == 2:
            return _to_power_basis11(nodes1, nodes2)

    raise NotImplementedError(
        'Degree 1', num_nodes1 - 1, 'Degree2', num_nodes2 - 1,
        'Currently only supporting degree pair 1-1')
