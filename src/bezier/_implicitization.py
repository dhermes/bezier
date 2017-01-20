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

"""Helper for implicitizing B |eacute| zier curves.

.. _resultant: https://en.wikipedia.org/wiki/Resultant
.. _algebraic curve: https://en.wikipedia.org/wiki/Algebraic_curve
.. _Farouki and Rajan: http://dx.doi.org/10.1016/0167-8396(88)90016-7

Primarily uses the `resultant`_ to evaluate the implicitized
`algebraic curve`_. In order to do this on B |eacute| zier curves
without translating to a power basis, we utilize the work of
`Farouki and Rajan`_ to compute a modified Sylvester determinant.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


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
    else:
        raise NotImplementedError('Only degrees 1 and 2 supported')
