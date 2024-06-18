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

r"""Helpers for intersecting B |eacute| zier curves via algebraic methods.

Primarily helps implicitize B |eacute| zier curves.

.. _resultant: https://en.wikipedia.org/wiki/Resultant
.. _algebraic curve: https://en.wikipedia.org/wiki/Algebraic_curve
.. _Farouki and Rajan: https://dx.doi.org/10.1016/0167-8396(88)90016-7
.. _theorem: https://en.wikipedia.org/wiki/B%C3%A9zout's_theorem

Uses the `resultant`_ to evaluate the implicitized
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
from numpy.polynomial import polynomial

from bezier import _curve_helpers
from bezier import _geometric_intersection
from bezier import _helpers
from bezier import _intersection_helpers
from bezier.hazmat import geometric_intersection
from bezier.hazmat import helpers as _py_helpers


# NOTE: These are hardcoded from:
#         0.5 * (np.cos(np.pi * np.arange(1, 14, 2.0) / 14) + 1.0)
#         0.5 * (np.cos(np.pi * np.arange(1, 18, 2.0) / 18) + 1.0)
#         0.5 * (np.cos(np.pi * np.arange(1, 20, 2.0) / 20) + 1.0)
#       See 73f909805e971221cb7976cf69603b82f31a4a32 and
#       80907b1be03f5895f7132e0e920c5e3cdebba9ac.
_CHEB7 = np.asfortranarray(
    [
        float.fromhex("0x1.f994e02ac74b4p-1"),
        float.fromhex("0x1.c8261ba82ef26p-1"),
        float.fromhex("0x1.6f130135c6af0p-1"),
        float.fromhex("0x1.0000000000000p-1"),
        float.fromhex("0x1.21d9fd9472a20p-2"),
        float.fromhex("0x1.becf22be886e0p-4"),
        float.fromhex("0x1.9ac7f54e2d2c0p-7"),
    ]
)
_CHEB9 = np.asfortranarray(
    [
        float.fromhex("0x1.fc1c5c6408e0cp-1"),
        float.fromhex("0x1.ddb3d742c2656p-1"),
        float.fromhex("0x1.a48dba91e0b0ep-1"),
        float.fromhex("0x1.578ea1d2282fep-1"),
        float.fromhex("0x1.0000000000000p-1"),
        float.fromhex("0x1.50e2bc5bafa08p-2"),
        float.fromhex("0x1.6dc915b87d3c6p-3"),
        float.fromhex("0x1.126145e9ecd5cp-4"),
        float.fromhex("0x1.f1d1cdfb8fa40p-8"),
    ]
)
_CHEB10 = np.asfortranarray(
    [
        float.fromhex("0x1.fcd924a17f22ep-1"),
        float.fromhex("0x1.e41900e9e9636p-1"),
        float.fromhex("0x1.b504f333f9de6p-1"),
        float.fromhex("0x1.7438b8ad13780p-1"),
        float.fromhex("0x1.280c16cf50a6fp-1"),
        float.fromhex("0x1.afe7d2615eb25p-2"),
        float.fromhex("0x1.178e8ea5d9100p-2"),
        float.fromhex("0x1.2bec333018868p-3"),
        float.fromhex("0x1.be6ff16169ca0p-5"),
        float.fromhex("0x1.936daf406e940p-8"),
    ]
)
# Allow a buffer of sqrt(sqrt(machine precision)) for polynomial roots.
_IMAGINARY_WIGGLE = 0.5**13
_UNIT_INTERVAL_WIGGLE_START = -(0.5**13)
_UNIT_INTERVAL_WIGGLE_END = 1.0 + 0.5**13
_SIGMA_THRESHOLD = 0.5**20
_SINGULAR_EPS = 0.5**52
# Detect almost zero polynomials.
_L2_THRESHOLD = 0.5**40  # 4096 (machine precision)
_ZERO_THRESHOLD = 0.5**38  # 16384 (machine precision)
_COEFFICIENT_THRESHOLD = 0.5**26  # sqrt(machine precision)
_NON_SIMPLE_THRESHOLD = 0.5**48  # 16 (machine precision)
_COINCIDENT_ERR = "Coincident curves not currently supported"
_NON_SIMPLE_ERR = "Polynomial has non-simple roots"
_POWER_BASIS_ERR = (
    "Currently only supporting degree pairs "
    "1-1, 1-2, 1-3, 1-4, 2-2, 2-3, 2-4 and 3-3."
)
_DISJOINT = geometric_intersection.BoxIntersectionType.DISJOINT


def _evaluate3(nodes, x_val, y_val):
    """Helper for :func:`evaluate` when ``nodes`` is degree 3.

    Args:
        nodes (numpy.ndarray): ``2 x 4`` array of nodes in a curve.
        x_val (float): ``x``-coordinate for evaluation.
        y_val (float): ``y``-coordinate for evaluation.

    Returns:
        float: The computed value of :math:`f(x, y)`.
    """
    # NOTE: This may be (a) slower and (b) less precise than
    #       hard-coding the determinant.
    sylvester_mat = np.zeros((6, 6), order="F")
    delta = nodes - np.asfortranarray([[x_val], [y_val]])
    delta[:, 1:3] *= 3.0
    # Swap rows/columns so that x-y are right next to each other.
    # This will only change the determinant up to a sign.
    sylvester_mat[:2, :4] = delta
    sylvester_mat[2:4, 1:5] = delta
    sylvester_mat[4:, 2:] = delta
    return np.linalg.det(sylvester_mat)


def evaluate(nodes, x_val, y_val):
    r"""Evaluate the implicitized bivariate polynomial containing the curve.

    Assumes the `algebraic curve`_ containing :math:`B(s, t)` is given by
    :math:`f(x, y) = 0`. This function evaluates :math:`f(x, y)`.

    .. note::

       This assumes, but doesn't check, that ``nodes`` has 2 rows.

    .. note::

       This assumes, but doesn't check, that ``nodes`` is not degree-elevated.
       If it were degree-elevated, then the Sylvester matrix will always
       have zero determinant.

    Args:
        nodes (numpy.ndarray): ``2 x N`` array of nodes in a curve.
        x_val (float): ``x``-coordinate for evaluation.
        y_val (float): ``y``-coordinate for evaluation.

    Returns:
        float: The computed value of :math:`f(x, y)`.

    Raises:
        ValueError: If the curve is a point.
        UnsupportedDegree: If the degree is not 1, 2 or 3.
    """
    _, num_nodes = nodes.shape
    if num_nodes == 1:
        raise ValueError("A point cannot be implicitized")

    if num_nodes == 2:
        # x(s) - x = (x0 - x) (1 - s) + (x1 - x) s
        # y(s) - y = (y0 - y) (1 - s) + (y1 - y) s
        # Modified Sylvester: [x0 - x, x1 - x]
        #                     [y0 - y, y1 - y]
        return (nodes[0, 0] - x_val) * (nodes[1, 1] - y_val) - (
            nodes[0, 1] - x_val
        ) * (nodes[1, 0] - y_val)

    if num_nodes == 3:
        # x(s) - x = (x0 - x) (1 - s)^2 + 2 (x1 - x) s(1 - s) + (x2 - x) s^2
        # y(s) - y = (y0 - y) (1 - s)^2 + 2 (y1 - y) s(1 - s) + (y2 - y) s^2
        # Modified Sylvester: [x0 - x, 2(x1 - x),    x2 - x,      0] = A|B|C|0
        #                     [     0,    x0 - x, 2(x1 - x), x2 - x]   0|A|B|C
        #                     [y0 - y, 2(y1 - y),    y2 - y,      0]   D|E|F|0
        #                     [     0,    y0 - y, 2(y1 - y), y2 - y]   0|D|E|F
        val_a, val_b, val_c = nodes[0, :] - x_val
        val_b *= 2
        val_d, val_e, val_f = nodes[1, :] - y_val
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

    if num_nodes == 4:
        return _evaluate3(nodes, x_val, y_val)

    raise _py_helpers.UnsupportedDegree(num_nodes - 1, supported=(1, 2, 3))


def eval_intersection_polynomial(nodes1, nodes2, t):
    r"""Evaluates a parametric curve **on** an implicitized algebraic curve.

    Uses :func:`evaluate` to evaluate :math:`f_1(x, y)`, the implicitization
    of ``nodes1``. Then plugs ``t`` into the second parametric curve to
    get an ``x``- and ``y``-coordinate and evaluate the
    **intersection polynomial**:

    .. math::

       g(t) = f_1\left(x_2(t), y_2(t)\right)

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.
        t (float): The parameter along ``nodes2`` where we evaluate
            the function.

    Returns:
        float: The computed value of :math:`f_1(x_2(t), y_2(t))`.
    """
    (x_val,), (y_val,) = _curve_helpers.evaluate_multi(
        nodes2, np.asfortranarray([t])
    )
    return evaluate(nodes1, x_val, y_val)


def _to_power_basis11(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

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
    # [1 0][c0] = [n0]
    # [1 1][c1]   [n1]
    val0 = eval_intersection_polynomial(nodes1, nodes2, 0.0)
    val1 = eval_intersection_polynomial(nodes1, nodes2, 1.0)
    # [c0] = [ 1 0][n0]
    # [c1] = [-1 1][n1]
    return np.asfortranarray([val0, -val0 + val1])


def _to_power_basis12(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that the first curve is
    degree one and the second is degree two. In this case, B |eacute| zout's
    `theorem`_ tells us that the **intersection polynomial** is
    degree :math:`1 \cdot 2` hence we return three coefficients.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``3``-array of coefficients.
    """
    # We manually invert the Vandermonde matrix:
    # [1 0   0  ][c0] = [n0]
    # [1 1/2 1/4][c1]   [n1]
    # [1 1   1  ][c2]   [n2]
    val0 = eval_intersection_polynomial(nodes1, nodes2, 0.0)
    val1 = eval_intersection_polynomial(nodes1, nodes2, 0.5)
    val2 = eval_intersection_polynomial(nodes1, nodes2, 1.0)
    # [c0] = [ 1  0  0][n0]
    # [c1] = [-3  4 -1][n1]
    # [c2] = [ 2 -4  2][n2]
    return np.asfortranarray(
        [
            val0,
            -3.0 * val0 + 4.0 * val1 - val2,
            2.0 * val0 - 4.0 * val1 + 2.0 * val2,
        ]
    )


def _to_power_basis13(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that the first curve is
    degree one and the second is degree three. In this case, B |eacute| zout's
    `theorem`_ tells us that the **intersection polynomial** is
    degree :math:`1 \cdot 3` hence we return four coefficients.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``4``-array of coefficients.
    """
    # We manually invert the Vandermonde matrix:
    # Use exact f.p. numbers to avoid round-off wherever possible.
    # [1 0   0    0    ][c0] = [n0]
    # [1 1/4 1/16 1/64 ][c1]   [n1]
    # [1 3/4 9/16 27/64][c2]   [n2]
    # [1 1   1    1    ][c3]   [n3]
    val0 = eval_intersection_polynomial(nodes1, nodes2, 0.0)
    val1 = eval_intersection_polynomial(nodes1, nodes2, 0.25)
    val2 = eval_intersection_polynomial(nodes1, nodes2, 0.75)
    val3 = eval_intersection_polynomial(nodes1, nodes2, 1.0)
    # [c0] =       [  3   0   0   0][n0]
    # [c1] = 1 / 3 [-19  24  -8   3][n1]
    # [c2] =       [ 32 -56  40 -16][n2]
    # [c3] =       [-16  32 -32  16][n3]
    # Since polynomial coefficients, we don't need to divide by 3
    # to get the same polynomial. Avoid the division to avoid round-off.
    return np.asfortranarray(
        [
            3.0 * val0,
            -19.0 * val0 + 24.0 * val1 - 8.0 * val2 + 3.0 * val3,
            32.0 * val0 - 56.0 * val1 + 40.0 * val2 - 16.0 * val3,
            -16.0 * val0 + 32.0 * val1 - 32.0 * val2 + 16.0 * val3,
        ]
    )


def _to_power_basis_degree4(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that B |eacute| zout's
    `theorem`_ tells us the **intersection polynomial** is degree
    :math:`4`. This happens if the two curves have degrees two and two
    or have degrees one and four.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``5``-array of coefficients.
    """
    # We manually invert the Vandermonde matrix:
    # [1 0   0    0     0     ][c0] = [n0]
    # [1 1/4 1/16 1/64  1/256 ][c1]   [n1]
    # [1 1/2 1/4  1/8   1/16  ][c2]   [n2]
    # [1 3/4 9/16 27/64 81/256][c3]   [n3]
    # [1 1   1    1     1     ][c4]   [n4]
    val0 = eval_intersection_polynomial(nodes1, nodes2, 0.0)
    val1 = eval_intersection_polynomial(nodes1, nodes2, 0.25)
    val2 = eval_intersection_polynomial(nodes1, nodes2, 0.5)
    val3 = eval_intersection_polynomial(nodes1, nodes2, 0.75)
    val4 = eval_intersection_polynomial(nodes1, nodes2, 1.0)
    # [c0] =       [ 3   0    0    0    0 ][n0]
    # [c1] = 1 / 3 [-25  48  -36   16  -3 ][n1]
    # [c2] =       [ 70 -208  228 -112  22][n2]
    # [c3] =       [-80  288 -384  224 -48][n3]
    # [c4] =       [ 32 -128  192 -128  32][n4]
    # Since polynomial coefficients, we don't need to divide by 3
    # to get the same polynomial. Avoid the division to avoid round-off.
    return np.asfortranarray(
        [
            3.0 * val0,
            -25.0 * val0
            + 48.0 * val1
            - 36.0 * val2
            + 16.0 * val3
            - 3.0 * val4,
            70.0 * val0
            - 208.0 * val1
            + 228.0 * val2
            - 112.0 * val3
            + 22.0 * val4,
            (
                -80.0 * val0
                + 288.0 * val1
                - 384.0 * val2
                + 224.0 * val3
                - 48.0 * val4
            ),
            32.0 * val0
            - 128.0 * val1
            + 192.0 * val2
            - 128.0 * val3
            + 32.0 * val4,
        ]
    )


def _to_power_basis23(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that the first curve is
    degree two and the second is degree three. In this case, B |eacute| zout's
    `theorem`_ tells us that the **intersection polynomial** is
    degree :math:`2 \cdot 3` hence we return seven coefficients.

    .. note::

       This uses a least-squares fit to the function evaluated at the
       Chebyshev nodes (scaled and shifted onto ``[0, 1]``). Hence, the
       coefficients may be less stable than those produced for smaller
       degrees.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``7``-array of coefficients.
    """
    evaluated = [
        eval_intersection_polynomial(nodes1, nodes2, t_val) for t_val in _CHEB7
    ]
    return polynomial.polyfit(_CHEB7, evaluated, 6)


def _to_power_basis_degree8(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that B |eacute| zout's
    `theorem`_ tells us the **intersection polynomial** is degree
    :math:`8`. This happens if the two curves have degrees one and eight
    or have degrees two and four.

    .. note::

       This uses a least-squares fit to the function evaluated at the
       Chebyshev nodes (scaled and shifted onto ``[0, 1]``). Hence, the
       coefficients may be less stable than those produced for smaller
       degrees.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``9``-array of coefficients.
    """
    evaluated = [
        eval_intersection_polynomial(nodes1, nodes2, t_val) for t_val in _CHEB9
    ]
    return polynomial.polyfit(_CHEB9, evaluated, 8)


def _to_power_basis33(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that each curve is
    degree three. In this case, B |eacute| zout's `theorem`_ tells us
    that the **intersection polynomial** is degree :math:`3 \cdot 3`
    hence we return ten coefficients.

    .. note::

       This uses a least-squares fit to the function evaluated at the
       Chebyshev nodes (scaled and shifted onto ``[0, 1]``). Hence, the
       coefficients may be less stable than those produced for smaller
       degrees.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``10``-array of coefficients.
    """
    evaluated = [
        eval_intersection_polynomial(nodes1, nodes2, t_val)
        for t_val in _CHEB10
    ]
    return polynomial.polyfit(_CHEB10, evaluated, 9)


def to_power_basis(nodes1, nodes2):
    """Compute the coefficients of an **intersection polynomial**.

    .. note::

       This requires that the degree of the curve given by ``nodes1`` is
       less than or equal to the degree of that given by ``nodes2``.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: Array of coefficients.

    Raises:
        NotImplementedError: If the degree pair is not ``1-1``, ``1-2``,
            ``1-3``, ``1-4``, ``2-2``, ``2-3``, ``2-4`` or ``3-3``.
    """
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-return-statements
    _, num_nodes1 = nodes1.shape
    _, num_nodes2 = nodes2.shape
    if num_nodes1 == 2:
        if num_nodes2 == 2:
            return _to_power_basis11(nodes1, nodes2)

        elif num_nodes2 == 3:
            return _to_power_basis12(nodes1, nodes2)

        elif num_nodes2 == 4:
            return _to_power_basis13(nodes1, nodes2)

        elif num_nodes2 == 5:
            return _to_power_basis_degree4(nodes1, nodes2)

    elif num_nodes1 == 3:
        if num_nodes2 == 3:
            return _to_power_basis_degree4(nodes1, nodes2)

        elif num_nodes2 == 4:
            return _to_power_basis23(nodes1, nodes2)

        elif num_nodes2 == 5:
            return _to_power_basis_degree8(nodes1, nodes2)

    elif num_nodes1 == 4:
        if num_nodes2 == 4:
            return _to_power_basis33(nodes1, nodes2)

    raise NotImplementedError(
        "Degree 1",
        num_nodes1 - 1,
        "Degree 2",
        num_nodes2 - 1,
        _POWER_BASIS_ERR,
    )


def polynomial_norm(coeffs):
    r"""Computes :math:`L_2` norm of polynomial on :math:`\left[0, 1\right]`.

    We have

    .. math::

       \left\langle f, f \right\rangle = \sum_{i, j}
           \int_0^1 c_i c_j x^{i + j} \, dx = \sum_{i, j}
           \frac{c_i c_j}{i + j + 1} = \sum_{i} \frac{c_i^2}{2 i + 1}
           + 2 \sum_{j > i} \frac{c_i c_j}{i + j + 1}.

    Args:
        coeffs (numpy.ndarray): ``d + 1``-array of coefficients in monomial /
            power basis.

    Returns:
        float: The :math:`L_2` norm of the polynomial.
    """
    (num_coeffs,) = coeffs.shape
    result = 0.0
    for i in range(num_coeffs):
        coeff_i = coeffs[i]
        result += coeff_i * coeff_i / (2.0 * i + 1.0)
        for j in range(i + 1, num_coeffs):
            coeff_j = coeffs[j]
            result += 2.0 * coeff_i * coeff_j / (i + j + 1.0)
    return np.sqrt(result)


def normalize_polynomial(coeffs, threshold=_L2_THRESHOLD):
    r"""Normalizes a polynomial in the :math:`L_2` sense.

    Does so on the interval :math:`\left[0, 1\right]` via
    :func:`polynomial_norm`.

    Args:
        coeffs (numpy.ndarray): ``d + 1``-array of coefficients in monomial /
            power basis.
        threshold (Optional[float]): The point :math:`\tau` below which a
            polynomial will be considered to be numerically equal to zero,
            applies to all :math:`f` with :math:`\| f \|_{L_2} < \tau`.

    Returns:
        numpy.ndarray: The normalized polynomial.
    """
    l2_norm = polynomial_norm(coeffs)
    if l2_norm < threshold:
        return np.zeros(coeffs.shape, order="F")

    else:
        coeffs /= l2_norm
        return coeffs


def _get_sigma_coeffs(coeffs):
    r"""Compute "transformed" form of polynomial in Bernstein form.

    Makes sure the "transformed" polynomial is monic (i.e. by dividing by the
    lead coefficient) and just returns the coefficients of the non-lead terms.

    For example, the polynomial

    .. math::

       f(s) = 4 (1 - s)^2 + 3 \cdot 2 s(1 - s) + 7 s^2

    can be transformed into :math:`f(s) = 7 (1 - s)^2 g\left(\sigma\right)`
    for :math:`\sigma = \frac{s}{1 - s}` and :math:`g(\sigma)` a
    polynomial in the power basis:

    .. math::

       g(\sigma) = \frac{4}{7} + \frac{6}{7} \sigma + \sigma^2.

    In cases where terms with "low exponents" of :math:`(1 - s)` have
    coefficient zero, the degree of :math:`g(\sigma)` may not be the
    same as the degree of :math:`f(s)`:

    .. math::

       \begin{align*}
       f(s) &= 5 (1 - s)^4 - 3 \cdot 4 s(1 - s)^3 + 11 \cdot 6 s^2 (1 - s)^2 \\
       \Longrightarrow g(\sigma) &= \frac{5}{66} - \frac{2}{11} \sigma +
           \sigma^2.
       \end{align*}

    Args:
        coeffs (numpy.ndarray): A 1D array of coefficients in
            the Bernstein basis.

    Returns:
        Tuple[Optional[numpy.ndarray], int, int]: A triple of

        * 1D array of the transformed coefficients (will be unset if
          ``effective_degree == 0``)
        * the "full degree" based on the size of ``coeffs``
        * the "effective degree" determined by the number of leading zeros.
    """
    (num_nodes,) = coeffs.shape
    degree = num_nodes - 1
    effective_degree = None
    for index in range(degree, -1, -1):
        if coeffs[index] != 0.0:
            effective_degree = index
            break

    if effective_degree is None:
        # NOTE: This means ``np.all(coeffs == 0.0)``.
        return None, 0, 0

    if effective_degree == 0:
        return None, degree, 0

    sigma_coeffs = coeffs[:effective_degree] / coeffs[effective_degree]
    # Now we need to add the binomial coefficients, but we avoid
    # computing actual binomial coefficients. Starting from the largest
    # exponent, the first ratio of binomial coefficients is
    #       (d C (e - 1)) / (d C e)
    #     = e / (d - e + 1)
    binom_numerator = effective_degree
    binom_denominator = degree - effective_degree + 1
    for exponent in range(effective_degree - 1, -1, -1):
        sigma_coeffs[exponent] *= binom_numerator
        sigma_coeffs[exponent] /= binom_denominator
        # We swap (d C j) with (d C (j - 1)), so `p / q` becomes
        #       (p / q) (d C (j - 1)) / (d C j)
        #     = (p j) / (q (d - j + 1))
        binom_numerator *= exponent
        binom_denominator *= degree - exponent + 1
    return sigma_coeffs, degree, effective_degree


def bernstein_companion(coeffs):
    r"""Compute a companion matrix for a polynomial in Bernstein basis.

    .. note::

       This assumes the caller passes in a 1D array but does not check.

    This takes the polynomial

    .. math::

       f(s) = \sum_{j = 0}^n b_{j, n} \cdot C_j.

    and uses the variable :math:`\sigma = \frac{s}{1 - s}` to rewrite as

    .. math::

       f(s) = (1 - s)^n \sum_{j = 0}^n \binom{n}{j} C_j \sigma^j.

    This converts the Bernstein coefficients :math:`C_j` into "generalized
    Bernstein" coefficients :math:`\widetilde{C_j} = \binom{n}{j} C_j`.

    Args:
        coeffs (numpy.ndarray): A 1D array of coefficients in
            the Bernstein basis.

    Returns:
        Tuple[numpy.ndarray, int, int]: A triple of

        * 2D NumPy array with the companion matrix.
        * the "full degree" based on the size of ``coeffs``
        * the "effective degree" determined by the number of leading zeros.
    """
    sigma_coeffs, degree, effective_degree = _get_sigma_coeffs(coeffs)
    if effective_degree == 0:
        return np.empty((0, 0), order="F"), degree, 0

    companion = np.zeros((effective_degree, effective_degree), order="F")
    # pylint: disable=unsupported-assignment-operation
    companion.flat[effective_degree :: effective_degree + 1] = (  # noqa: E203
        1.0
    )
    # pylint: enable=unsupported-assignment-operation
    companion[0, :] = -sigma_coeffs[::-1]
    return companion, degree, effective_degree


def bezier_roots(coeffs):
    r"""Compute polynomial roots from a polynomial in the Bernstein basis.

    .. note::

       This assumes the caller passes in a 1D array but does not check.

    This takes the polynomial

    .. math::

       f(s) = \sum_{j = 0}^n b_{j, n} \cdot C_j.

    and uses the variable :math:`\sigma = \frac{s}{1 - s}` to rewrite as

    .. math::

       \begin{align*}
       f(s) &= (1 - s)^n \sum_{j = 0}^n \binom{n}{j} C_j \sigma^j \\
            &= (1 - s)^n \sum_{j = 0}^n \widetilde{C_j} \sigma^j.
       \end{align*}

    Then it uses an eigenvalue solver to find the roots of

    .. math::

       g(\sigma) = \sum_{j = 0}^n \widetilde{C_j} \sigma^j

    and convert them back into roots of :math:`f(s)` via
    :math:`s = \frac{\sigma}{1 + \sigma}`.

    For example, consider

    .. math::

       \begin{align*}
       f_0(s) &= 2 (2 - s)(3 + s) \\
              &= 12(1 - s)^2 + 11 \cdot 2s(1 - s) + 8 s^2
       \end{align*}

    First, we compute the companion matrix for

    .. math::

       g_0(\sigma) = 12 + 22 \sigma + 8 \sigma^2

    .. testsetup:: bezier-roots0, bezier-roots1, bezier-roots2

       import numpy as np
       import numpy.linalg
       from bezier.hazmat.algebraic_intersection import bernstein_companion
       from bezier.hazmat.algebraic_intersection import bezier_roots

       machine_eps = np.finfo(np.float64).eps

    .. doctest:: bezier-roots0

       >>> import numpy as np
       >>> coeffs0 = np.asfortranarray([12.0, 11.0, 8.0])
       >>> companion0, _, _ = bernstein_companion(coeffs0)
       >>> companion0
       array([[-2.75, -1.5 ],
              [ 1.  ,  0.  ]])

    then take the eigenvalues of the companion matrix:

    .. doctest:: bezier-roots0

       >>> sigma_values0 = np.linalg.eigvals(companion0)
       >>> sigma_values0
       array([-2.  , -0.75])

    after transforming them, we have the roots of :math:`f(s)`:

    .. doctest:: bezier-roots0

       >>> sigma_values0 / (1.0 + sigma_values0)
       array([ 2., -3.])
       >>> bezier_roots(coeffs0)
       array([ 2., -3.])

    In cases where :math:`s = 1` is a root, the lead coefficient of
    :math:`g` would be :math:`0`, so there is a reduction in the
    companion matrix.

    .. math::

       \begin{align*}
       f_1(s) &= 6 (s - 1)^2 (s - 3) (s - 5) \\
              &= 90 (1 - s)^4 + 33 \cdot 4s(1 - s)^3 + 8 \cdot 6s^2(1 - s)^2
       \end{align*}

    .. doctest:: bezier-roots1
       :options: +NORMALIZE_WHITESPACE

       >>> coeffs1 = np.asfortranarray([90.0, 33.0, 8.0, 0.0, 0.0])
       >>> companion1, degree1, effective_degree1 = bernstein_companion(
       ...     coeffs1)
       >>> companion1
       array([[-2.75 , -1.875],
              [ 1.   ,  0.   ]])
       >>> degree1
       4
       >>> effective_degree1
       2

    so the roots are a combination of the roots determined from
    :math:`s = \frac{\sigma}{1 + \sigma}` and the number of factors
    of :math:`(1 - s)` (i.e. the difference between the degree and
    the effective degree):

    .. doctest:: bezier-roots1

       >>> bezier_roots(coeffs1)
       array([3., 5., 1., 1.])

    In some cases, a polynomial is represented with an "elevated" degree:

    .. math::

       \begin{align*}
       f_2(s) &= 3 (s^2 + 1) \\
              &= 3 (1 - s)^3 + 3 \cdot 3s(1 - s)^2 +
                 4 \cdot 3s^2(1 - s) + 6 s^3
       \end{align*}

    This results in a "point at infinity"
    :math:`\sigma = -1 \Longleftrightarrow s = \infty`:

    .. doctest:: bezier-roots2

       >>> coeffs2 = np.asfortranarray([3.0, 3.0, 4.0, 6.0])
       >>> companion2, _, _ = bernstein_companion(coeffs2)
       >>> companion2
       array([[-2. , -1.5, -0.5],
              [ 1. ,  0. ,  0. ],
              [ 0. ,  1. ,  0. ]])
       >>> sigma_values2 = np.linalg.eigvals(companion2)
       >>> sigma_values2
       array([-1. +0.j , -0.5+0.5j, -0.5-0.5j])


    so we drop any values :math:`\sigma` that are sufficiently close to
    :math:`-1`:

    .. doctest:: bezier-roots2

       >>> expected2 = np.asfortranarray([1.0j, -1.0j])
       >>> roots2 = bezier_roots(coeffs2)
       >>> np.allclose(expected2, roots2, rtol=8 * machine_eps, atol=0.0)
       True

    Args:
        coeffs (numpy.ndarray): A 1D array of coefficients in
            the Bernstein basis.

    Returns:
        numpy.ndarray: A 1D array containing the roots.
    """
    companion, degree, effective_degree = bernstein_companion(coeffs)
    if effective_degree:
        sigma_roots = np.linalg.eigvals(companion)
        # Filter out `sigma = -1`, i.e. "points at infinity".
        # We want the error ||(sigma - (-1))|| ~= 2^{-52}
        to_keep = np.abs(sigma_roots + 1.0) > _SIGMA_THRESHOLD
        sigma_roots = sigma_roots[to_keep]
        s_vals = sigma_roots / (1.0 + sigma_roots)
    else:
        s_vals = np.empty((0,), order="F")
    if effective_degree != degree:
        delta = degree - effective_degree
        s_vals = np.hstack([s_vals, [1] * delta])
    return s_vals


def lu_companion(top_row, value):
    r"""Compute an LU-factored :math:`C - t I` and its 1-norm.

    .. _dgecon:
        http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga188b8d30443d14b1a3f7f8331d87ae60.html#ga188b8d30443d14b1a3f7f8331d87ae60
    .. _dgetrf:
        http://www.netlib.org/lapack/explore-html/dd/d9a/group__double_g_ecomputational_ga0019443faea08275ca60a734d0593e60.html#ga0019443faea08275ca60a734d0593e60

    .. note::

       The output of this function is intended to be used with `dgecon`_
       from LAPACK. ``dgecon`` expects both the 1-norm of the matrix
       and expects the matrix to be passed in an already LU-factored form
       (via `dgetrf`_).

    The companion matrix :math:`C` is given by the ``top_row``, for
    example, the polynomial :math:`t^3 + 3 t^2 - t + 2` has
    a top row of ``-3, 1, -2`` and the corresponding companion matrix
    is:

    .. math::

       \left[\begin{array}{c c c}
         -3 & 1 & -2 \\
          1 & 0 &  0 \\
          0 & 1 &  0
       \end{array}\right]

    After doing a full cycle of the rows (shifting the first to the last and
    moving all other rows up), row reduction of :math:`C - t I` yields

    .. math::

       \left[\begin{array}{c c c}
          1     & -t &  0 \\
          0     &  1 & -t \\
         -3 - t &  1 & -2
       \end{array}\right] =
       \left[\begin{array}{c c c}
          1     & 0             & 0 \\
          0     & 1             & 0 \\
         -3 - t & 1 + t(-3 - t) & 1
       \end{array}\right]
       \left[\begin{array}{c c c}
          1 & -t &  0                    \\
          0 &  1 & -t                    \\
          0 &  0 & -2 + t(1 + t(-3 - t))
       \end{array}\right]

    and in general, the terms in the bottom row correspond to the intermediate
    values involved in evaluating the polynomial via `Horner's method`_.

    .. _Horner's method: https://en.wikipedia.org/wiki/Horner%27s_method

    .. testsetup:: lu-companion

       import numpy as np
       import numpy.linalg
       from bezier.hazmat.algebraic_intersection import lu_companion

    .. doctest:: lu-companion

       >>> top_row = np.asfortranarray([-3.0, 1.0, -2.0])
       >>> t_val = 0.5
       >>> lu_mat, one_norm = lu_companion(top_row, t_val)
       >>> lu_mat
       array([[ 1.   , -0.5  ,  0.   ],
              [ 0.   ,  1.   , -0.5  ],
              [-3.5  , -0.75 , -2.375]])
       >>> float(one_norm)
       4.5
       >>> l_mat = np.tril(lu_mat, k=-1) + np.eye(3)
       >>> u_mat = np.triu(lu_mat)
       >>> a_mat = l_mat.dot(u_mat)
       >>> a_mat
       array([[ 1. , -0.5,  0. ],
              [ 0. ,  1. , -0.5],
              [-3.5,  1. , -2. ]])
       >>> float(np.linalg.norm(a_mat, ord=1))
       4.5

    Args:
        top_row (numpy.ndarray): 1D array, top row of companion matrix.
        value (float): The :math:`t` value used to form :math:`C - t I`.

    Returns:
        Tuple[numpy.ndarray, float]: Pair of

        * 2D array of LU-factored form of :math:`C - t I`, with the
          non-diagonal part of :math:`L` stored in the strictly lower triangle
          and :math:`U` stored in the upper triangle (we skip the permutation
          matrix, as it won't impact the 1-norm)
        * the 1-norm the matrix :math:`C - t I`

        As mentioned above, these two values are meant to be used with
        `dgecon`_.
    """
    (degree,) = top_row.shape
    lu_mat = np.zeros((degree, degree), order="F")
    if degree == 1:
        lu_mat[0, 0] = top_row[0] - value
        return lu_mat, abs(lu_mat[0, 0])

    # Column 0: Special case since it doesn't have ``-t`` above the diagonal.
    horner_curr = top_row[0] - value
    one_norm = 1.0 + abs(horner_curr)
    lu_mat[0, 0] = 1.0
    lu_mat[degree - 1, 0] = horner_curr
    # Columns 1-(end - 1): Three values in LU and C - t I.
    abs_one_plus = 1.0 + abs(value)
    last_row = degree - 1
    for col in range(1, degree - 1):
        curr_coeff = top_row[col]
        horner_curr = value * horner_curr + curr_coeff
        one_norm = max(one_norm, abs_one_plus + abs(curr_coeff))
        lu_mat[col - 1, col] = -value
        lu_mat[col, col] = 1.0
        lu_mat[last_row, col] = horner_curr
    # Last Column: Special case since it doesn't have ``-1`` on the diagonal.
    curr_coeff = top_row[last_row]
    horner_curr = value * horner_curr + curr_coeff
    one_norm = max(one_norm, abs(value) + abs(curr_coeff))
    lu_mat[last_row - 1, last_row] = -value
    lu_mat[last_row, last_row] = horner_curr
    return lu_mat, one_norm


def _reciprocal_condition_number(lu_mat, one_norm):
    r"""Compute reciprocal condition number of a matrix.

    Args:
        lu_mat (numpy.ndarray): A 2D array of a matrix :math:`A` that has been
            LU-factored, with the non-diagonal part of :math:`L` stored in the
            strictly lower triangle and :math:`U` stored in the upper triangle.
        one_norm (float): The 1-norm of the original matrix :math:`A`.

    Returns:
        float: The reciprocal condition number of :math:`A`.

    Raises:
        RuntimeError: If the reciprocal 1-norm condition number could not
            be computed.
    """
    # NOTE: We import SciPy at runtime to avoid the import-time cost for users
    #       that don't need algebraic intersection helpers (e.g. if only using
    #       the geometric intersection strategy). The ``scipy`` import is a
    #       tad expensive.
    # pylint: disable=import-outside-toplevel,no-name-in-module
    import scipy.linalg.lapack

    # pylint: enable=import-outside-toplevel,no-name-in-module

    _dgecon = scipy.linalg.lapack.dgecon  # pylint: disable=no-member
    rcond, info = _dgecon(lu_mat, one_norm)
    if info != 0:
        raise RuntimeError(
            "The reciprocal 1-norm condition number could not be computed."
        )

    return rcond


def bezier_value_check(coeffs, s_val, rhs_val=0.0):
    r"""Check if a polynomial in the Bernstein basis evaluates to a value.

    This is intended to be used for root checking, i.e. for a polynomial
    :math:`f(s)` and a particular value :math:`s_{\ast}`:

       Is it true that :math:`f\left(s_{\ast}\right) = 0`?

    Does so by re-stating as a matrix rank problem. As in
    :func:`~bezier.hazmat.algebraic_intersection.bezier_roots`, we can rewrite

    .. math::

       f(s) = (1 - s)^n g\left(\sigma\right)

    for :math:`\sigma = \frac{s}{1 - s}` and :math:`g\left(\sigma\right)`
    written in the power / monomial basis. Now, checking if
    :math:`g\left(\sigma\right) = 0` is a matter of checking that
    :math:`\det\left(C_g - \sigma I\right) = 0` (where :math:`C_g` is the
    companion matrix of :math:`g`).

    Due to issues of numerical stability, we'd rather ask if
    :math:`C_g - \sigma I` is singular to numerical precision. A typical
    approach to this using the singular values (assuming :math:`C_g` is
    :math:`m \times m`) is that the matrix is singular if

    .. math::

       \sigma_m < m \varepsilon \sigma_1 \Longleftrightarrow
           \frac{1}{\kappa_2} < m \varepsilon

    (where :math:`\kappa_2` is the 2-norm condition number of the matrix).
    Since we also know that :math:`\kappa_2 < \kappa_1`, a stronger
    requirement would be

    .. math::

       \frac{1}{\kappa_1} < \frac{1}{\kappa_2} < m \varepsilon.

    This is more useful since it is **much** easier to compute the 1-norm
    condition number / reciprocal condition number (and methods in LAPACK are
    provided for doing so).

    Args:
        coeffs (numpy.ndarray): A 1D array of coefficients in
            the Bernstein basis representing a polynomial.
        s_val (float): The value to check on the polynomial:
            :math:`f(s) = r`.
        rhs_val (Optional[float]): The value to check that the polynomial
            evaluates to. Defaults to ``0.0``.

    Returns:
        bool: Indicates if :math:`f\left(s_{\ast}\right) = r` (where
        :math:`s_{\ast}` is ``s_val`` and :math:`r` is ``rhs_val``).
    """
    if s_val == 1.0:
        return coeffs[-1] == rhs_val

    shifted_coeffs = coeffs - rhs_val
    sigma_coeffs, _, effective_degree = _get_sigma_coeffs(shifted_coeffs)
    if effective_degree == 0:
        # This means that all coefficients except the ``(1 - s)^n``
        # term are zero, so we have ``f(s) = C (1 - s)^n``. Since we know
        # ``s != 1``, this can only be zero if ``C == 0``.
        return shifted_coeffs[0] == 0.0

    sigma_val = s_val / (1.0 - s_val)
    lu_mat, one_norm = lu_companion(-sigma_coeffs[::-1], sigma_val)
    rcond = _reciprocal_condition_number(lu_mat, one_norm)
    # "Is a root?" IFF Singular IF ``1/kappa_1 < m epsilon``
    return rcond < effective_degree * _SINGULAR_EPS


def roots_in_unit_interval(coeffs):
    r"""Compute roots of a polynomial in the unit interval.

    Args:
        coeffs (numpy.ndarray): A 1D array (size ``d + 1``) of coefficients in
            monomial / power basis.

    Returns:
        numpy.ndarray: ``N``-array of real values in :math:`\left[0, 1\right]`.
    """
    all_roots = polynomial.polyroots(coeffs)
    # Only keep roots inside or very near to the unit interval.
    all_roots = all_roots[
        (_UNIT_INTERVAL_WIGGLE_START < all_roots.real)
        & (all_roots.real < _UNIT_INTERVAL_WIGGLE_END)
    ]
    # Only keep roots with very small imaginary part. (Really only
    # keep the real parts.)
    real_inds = np.abs(all_roots.imag) < _IMAGINARY_WIGGLE
    return all_roots[real_inds].real


def _strip_leading_zeros(coeffs, threshold=_COEFFICIENT_THRESHOLD):
    r"""Strip leading zero coefficients from a polynomial.

    .. note::

       This assumes the polynomial :math:`f` defined by ``coeffs``
       has been normalized (via :func:`.normalize_polynomial`).

    Args:
        coeffs (numpy.ndarray): ``d + 1``-array of coefficients in monomial /
            power basis.
        threshold (Optional[float]): The point :math:`\tau` below which a
            a coefficient will be considered to be numerically zero.

    Returns:
        numpy.ndarray: The same coefficients without any unnecessary zero
        terms.
    """
    while np.abs(coeffs[-1]) < threshold:
        coeffs = coeffs[:-1]
    return coeffs


def _check_non_simple(coeffs):
    r"""Checks that a polynomial has no non-simple roots.

    Does so by computing the companion matrix :math:`A` of :math:`f'`
    and then evaluating the rank of :math:`B = f(A)`. If :math:`B` is not
    full rank, then :math:`f` and :math:`f'` have a shared factor.

    See: https://dx.doi.org/10.1016/0024-3795(70)90023-6

    .. note::

       This assumes that :math:`f \neq 0`.

    Args:
        coeffs (numpy.ndarray): ``d + 1``-array of coefficients in monomial /
            power basis.

    Raises:
        NotImplementedError: If the polynomial has non-simple roots.
    """
    coeffs = _strip_leading_zeros(coeffs)
    (num_coeffs,) = coeffs.shape
    if num_coeffs < 3:
        return

    deriv_poly = polynomial.polyder(coeffs)
    companion = polynomial.polycompanion(deriv_poly)
    # NOTE: ``polycompanion()`` returns a C-contiguous array.
    companion = companion.T
    # Use Horner's method to evaluate f(companion)
    num_companion, _ = companion.shape
    id_mat = np.eye(num_companion, order="F")
    evaluated = coeffs[-1] * id_mat
    for index in range(num_coeffs - 2, -1, -1):
        coeff = coeffs[index]
        evaluated = (
            _py_helpers.matrix_product(evaluated, companion) + coeff * id_mat
        )
    if num_companion == 1:
        # NOTE: This relies on the fact that coeffs is normalized.
        if np.abs(evaluated[0, 0]) > _NON_SIMPLE_THRESHOLD:
            rank = 1
        else:
            rank = 0
    else:
        rank = np.linalg.matrix_rank(evaluated)
    if rank < num_companion:
        raise NotImplementedError(_NON_SIMPLE_ERR, coeffs)


def _resolve_and_add(nodes1, s_val, final_s, nodes2, t_val, final_t):
    """Resolve a computed intersection and add to lists.

    We perform one Newton step to deal with any residual issues of
    high-degree polynomial solves (one of which depends on the already
    approximate ``x_val, y_val``).

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        s_val (float): The approximate intersection parameter
            along ``nodes1``.
        final_s (List[float]): The list of accepted intersection
            parameters ``s``.
        nodes2 (numpy.ndarray): The nodes in the second curve.
        t_val (float): The approximate intersection parameter
            along ``nodes2``.
        final_t (List[float]): The list of accepted intersection
            parameters ``t``.
    """
    s_val, t_val = _intersection_helpers.newton_refine(
        s_val, nodes1, t_val, nodes2
    )
    s_val, success_s = _helpers.wiggle_interval(s_val)
    t_val, success_t = _helpers.wiggle_interval(t_val)
    if not (success_s and success_t):
        return

    final_s.append(s_val)
    final_t.append(t_val)


def intersect_curves(nodes1, nodes2):
    r"""Intersect two parametric B |eacute| zier curves.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``2 x N`` array of intersection parameters.
        Each row contains a pair of values :math:`s` and :math:`t`
        (each in :math:`\left[0, 1\right]`) such that the curves
        intersect: :math:`B_1(s) = B_2(t)`.

    Raises:
        NotImplementedError: If the "intersection polynomial" is
            all zeros -- which indicates coincident curves.
    """
    nodes1 = _curve_helpers.full_reduce(nodes1)
    nodes2 = _curve_helpers.full_reduce(nodes2)
    _, num_nodes1 = nodes1.shape
    _, num_nodes2 = nodes2.shape
    swapped = False
    if num_nodes1 > num_nodes2:
        nodes1, nodes2 = nodes2, nodes1
        swapped = True
    coeffs = normalize_polynomial(to_power_basis(nodes1, nodes2))
    if np.all(coeffs == 0.0):
        raise NotImplementedError(_COINCIDENT_ERR)

    _check_non_simple(coeffs)
    t_vals = roots_in_unit_interval(coeffs)
    final_s = []
    final_t = []
    for t_val in t_vals:
        (x_val,), (y_val,) = _curve_helpers.evaluate_multi(
            nodes2, np.asfortranarray([t_val])
        )
        s_val = locate_point(nodes1, x_val, y_val)
        if s_val is not None:
            _resolve_and_add(nodes1, s_val, final_s, nodes2, t_val, final_t)
    result = np.zeros((2, len(final_s)), order="F")
    if swapped:
        final_s, final_t = final_t, final_s
    result[0, :] = final_s
    result[1, :] = final_t
    return result


def poly_to_power_basis(bezier_coeffs):
    """Convert a B |eacute| zier curve to polynomial in power basis.

    .. note::

       This assumes, but does not verify, that the "B |eacute| zier
       degree" matches the true degree of the curve. Callers can
       guarantee this by calling :func:`.full_reduce`.

    Args:
        bezier_coeffs (numpy.ndarray): A 1D array of coefficients in
            the Bernstein basis.

    Returns:
        numpy.ndarray: 1D array of coefficients in monomial basis.

    Raises:
        UnsupportedDegree: If the degree of the curve is not among
            0, 1, 2 or 3.
    """
    (num_coeffs,) = bezier_coeffs.shape
    if num_coeffs == 1:
        return bezier_coeffs

    elif num_coeffs == 2:
        # C0 (1 - s) + C1 s = C0 + (C1 - C0) s
        coeff0, coeff1 = bezier_coeffs
        return np.asfortranarray([coeff0, coeff1 - coeff0])

    elif num_coeffs == 3:
        #   C0 (1 - s)^2 + C1 2 (1 - s) s + C2 s^2
        # = C0 + 2(C1 - C0) s + (C2 - 2 C1 + C0) s^2
        coeff0, coeff1, coeff2 = bezier_coeffs
        return np.asfortranarray(
            [coeff0, 2.0 * (coeff1 - coeff0), coeff2 - 2.0 * coeff1 + coeff0]
        )

    elif num_coeffs == 4:
        #   C0 (1 - s)^3 + C1 3 (1 - s)^2 + C2 3 (1 - s) s^2 + C3 s^3
        # = C0 + 3(C1 - C0) s + 3(C2 - 2 C1 + C0) s^2 +
        #   (C3 - 3 C2 + 3 C1 - C0) s^3
        coeff0, coeff1, coeff2, coeff3 = bezier_coeffs
        return np.asfortranarray(
            [
                coeff0,
                3.0 * (coeff1 - coeff0),
                3.0 * (coeff2 - 2.0 * coeff1 + coeff0),
                coeff3 - 3.0 * coeff2 + 3.0 * coeff1 - coeff0,
            ]
        )

    else:
        raise _py_helpers.UnsupportedDegree(
            num_coeffs - 1, supported=(0, 1, 2, 3)
        )


def locate_point(nodes, x_val, y_val):
    r"""Find the parameter corresponding to a point on a curve.

    .. note::

       This assumes that the curve :math:`B(s, t)` defined by ``nodes``
       lives in :math:`\mathbf{R}^2`.

    Args:
        nodes (numpy.ndarray): The nodes defining a B |eacute| zier curve.
        x_val (float): The :math:`x`-coordinate of the point.
        y_val (float): The :math:`y`-coordinate of the point.

    Returns:
        Optional[float]: The parameter on the curve (if it exists).
    """
    # First, reduce to the true degree of x(s) and y(s).
    zero1 = _curve_helpers.full_reduce(nodes[[0], :]) - x_val
    zero2 = _curve_helpers.full_reduce(nodes[[1], :]) - y_val
    # Make sure we have the lowest degree in front, to make the polynomial
    # solve have the fewest number of roots.
    if zero1.shape[1] > zero2.shape[1]:
        zero1, zero2 = zero2, zero1
    # If the "smallest" is a constant, we can't find any roots from it.
    if zero1.shape[1] == 1:
        # NOTE: We assume that callers won't pass ``nodes`` that are
        #       degree 0, so if ``zero1`` is a constant, ``zero2`` won't be.
        zero1, zero2 = zero2, zero1
    power_basis1 = poly_to_power_basis(zero1[0, :])
    all_roots = roots_in_unit_interval(power_basis1)
    if all_roots.size == 0:
        return None

    # NOTE: We normalize ``power_basis2`` because we want to check for
    #       "zero" values, i.e. f2(s) == 0.
    power_basis2 = normalize_polynomial(poly_to_power_basis(zero2[0, :]))
    near_zero = np.abs(polynomial.polyval(all_roots, power_basis2))
    index = np.argmin(near_zero)
    if near_zero[index] < _ZERO_THRESHOLD:
        return all_roots[index]

    return None


def all_intersections(nodes_first, nodes_second):
    r"""Find the points of intersection among a pair of curves.

    .. note::

       This assumes both curves are in :math:`\mathbf{R}^2`, but does not
       **explicitly** check this. However, functions used here will fail if
       that assumption fails.

    Args:
        nodes_first (numpy.ndarray): Control points of a curve to be
            intersected with ``nodes_second``.
        nodes_second (numpy.ndarray): Control points of a curve to be
            intersected with ``nodes_first``.

    Returns:
        Tuple[numpy.ndarray, bool]: An array and a flag:

        * A ``2 x N`` array of intersection parameters.
          Each row contains a pair of values :math:`s` and :math:`t`
          (each in :math:`\left[0, 1\right]`) such that the curves
          intersect: :math:`B_1(s) = B_2(t)`.
        * Flag indicating if the curves are coincident. (For now, this
          will always be :data:`False` since :func:`.intersect_curves`
          fails if coincident curves are detected.)
    """
    # Only attempt this if the bounding boxes intersect.
    bbox_int = _geometric_intersection.bbox_intersect(
        nodes_first, nodes_second
    )
    if bbox_int == _DISJOINT:
        return np.empty((2, 0), order="F"), False

    return intersect_curves(nodes_first, nodes_second), False
