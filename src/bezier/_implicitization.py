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
.. _theorem: https://en.wikipedia.org/wiki/B%C3%A9zout's_theorem

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
from numpy.polynomial import chebyshev
from numpy.polynomial import polynomial
import six

from bezier import _curve_helpers
from bezier import _intersection_helpers


_CHEB7, _ = chebyshev.chebgauss(7)
_CHEB7 = 0.5 * (_CHEB7 + 1.0)
_CHEB9, _ = chebyshev.chebgauss(9)
_CHEB9 = 0.5 * (_CHEB9 + 1.0)
_CHEB10, _ = chebyshev.chebgauss(10)
_CHEB10 = 0.5 * (_CHEB10 + 1.0)
# Allow a buffer of sqrt(sqrt(machine precision)) for polynomial roots.
_IMAGINARY_WIGGLE = 0.5**13
_UNIT_INTERVAL_WIGGLE_START = -0.5**13
_UNIT_INTERVAL_WIGGLE_END = 1.0 + 0.5**13
# Detect almost zero polynomials.
_L2_THRESHOLD = 0.5**40  # 4096 (machine precision)
_ZERO_THRESHOLD = 0.5**38  # 16384 (machine precision)
_COEFFICIENT_THRESHOLD = 0.5**26  # sqrt(machine precision)
_PARAM_THRESHOLD = 0.5**51  # 2 (machine precision)
_NON_SIMPLE_THRESHOLD = 0.5**48  # 16 (machine precision)
_COINCIDENT_ERR = 'Coincident curves not currently supported'
_NON_SIMPLE_ERR = 'Polynomial has non-simple roots'
_POWER_BASIS_ERR = (
    'Currently only supporting degree pairs '
    '1-1, 1-2, 1-3, 1-4, 2-2, 2-3, 2-4 and 3-3.')


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
    val0 = eval_intersection_polynomial(nodes1, nodes2, 0.0)
    val1 = eval_intersection_polynomial(nodes1, nodes2, 1.0)
    # [c0] = [ 1 0][n0]
    # [c1] = [-1 1][n1]
    return np.asfortranarray([val0, -val0 + val1])


def _to_power_basis12(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that the first curve is
    degree one and the second is degree two. In this case, B |eacute|
    zout's `theorem`_ tells us that the **intersection polynomial** is
    degree :math:`1 \cdot 2` hence we return three coefficients.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``3``-array of coefficients.
    """
    # We manually invert the Vandermonde matrix:
    # [1 0.0 0.0 ][c0] = [n0]
    # [1 0.5 0.25][c1]   [n1]
    # [1 1.0 1.0 ][c2]   [n2]
    val0 = eval_intersection_polynomial(nodes1, nodes2, 0.0)
    val1 = eval_intersection_polynomial(nodes1, nodes2, 0.5)
    val2 = eval_intersection_polynomial(nodes1, nodes2, 1.0)
    # [c0] = [ 1  0  0][n0]
    # [c1] = [-3  4 -1][n1]
    # [c2] = [ 2 -4  2][n2]
    return np.asfortranarray([
        val0,
        -3.0 * val0 + 4.0 * val1 - val2,
        2.0 * val0 - 4.0 * val1 + 2.0 * val2,
    ])


def _to_power_basis13(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that the first curve is
    degree one and the second is degree three. In this case, B |eacute|
    zout's `theorem`_ tells us that the **intersection polynomial** is
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
    return np.asfortranarray([
        3.0 * val0,
        -19.0 * val0 + 24.0 * val1 - 8.0 * val2 + 3.0 * val3,
        32.0 * val0 - 56.0 * val1 + 40.0 * val2 - 16.0 * val3,
        -16.0 * val0 + 32.0 * val1 - 32.0 * val2 + 16.0 * val3,
    ])


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
    return np.asfortranarray([
        3.0 * val0,
        -25.0 * val0 + 48.0 * val1 - 36.0 * val2 + 16.0 * val3 - 3.0 * val4,
        70.0 * val0 - 208.0 * val1 + 228.0 * val2 - 112.0 * val3 + 22.0 * val4,
        (-80.0 * val0 + 288.0 * val1 - 384.0 * val2 +
         224.0 * val3 - 48.0 * val4),
        32.0 * val0 - 128.0 * val1 + 192.0 * val2 - 128.0 * val3 + 32.0 * val4,
    ])


def _to_power_basis23(nodes1, nodes2):
    r"""Compute the coefficients of an **intersection polynomial**.

    Helper for :func:`to_power_basis` in the case that the first curve is
    degree two and the second is degree three. In this case, B |eacute|
    zout's `theorem`_ tells us that the **intersection polynomial** is
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
    evaluated = [eval_intersection_polynomial(nodes1, nodes2, t_val)
                 for t_val in _CHEB7]
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
    evaluated = [eval_intersection_polynomial(nodes1, nodes2, t_val)
                 for t_val in _CHEB9]
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
    evaluated = [eval_intersection_polynomial(nodes1, nodes2, t_val)
                 for t_val in _CHEB10]
    return polynomial.polyfit(_CHEB10, evaluated, 9)


def to_power_basis(nodes1, nodes2):
    """Compute the coefficients of an **intersection polynomial**.

    .. note::

       This assumes that the degree of the curve given by ``nodes1`` is
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
    # pylint: disable=too-many-return-statements
    num_nodes1, _ = nodes1.shape
    num_nodes2, _ = nodes2.shape
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
        'Degree 1', num_nodes1 - 1, 'Degree2', num_nodes2 - 1,
        _POWER_BASIS_ERR)
    # pylint: enable=too-many-return-statements


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
    num_coeffs, = coeffs.shape
    result = 0.0
    for i in six.moves.xrange(num_coeffs):
        coeff_i = coeffs[i]
        result += coeff_i * coeff_i / (2.0 * i + 1.0)
        for j in six.moves.xrange(i + 1, num_coeffs):
            coeff_j = coeffs[j]
            result += 2.0 * coeff_i * coeff_j / (i + j + 1.0)

    return np.sqrt(result)


def normalize_polynomial(coeffs, threshold=_L2_THRESHOLD):
    r"""Normalizes a polynomial in the :math:`L_2` sense.

    Does so on the interval :math:\left[0, 1\right]` via
    :func:`polynomial_norm`.

    Args:
        coeffs (numpy.ndarray): ``d + 1``-array of coefficients in monomial /
            power basis.
        threshold (Optional[float]): The point :math:`\tau` below which a
            polynomial will be considered to be numerically equal to zero,
            applies to all :math:`f` with :math`\| f \|_{L_2} < \tau`.

    Returns:
        numpy.ndarray: The normalized polynomial.
    """
    l2_norm = polynomial_norm(coeffs)
    if l2_norm < threshold:
        return np.zeros(coeffs.shape, order='F')
    else:
        coeffs /= l2_norm
        return coeffs


def roots_in_unit_interval(coeffs):
    r"""Compute roots of a polynomial in the unit interval.

    Args:
        coeffs (numpy.ndarray): ``d + 1``-array of coefficients in monomial /
            power basis.

    Returns:
        numpy.ndarray: ``N``-array of real values in :math:`\left[0, 1\right]`.
    """
    all_roots = polynomial.polyroots(coeffs)
    # Only keep roots inside or very near to the unit interval.
    all_roots = all_roots[
        (_UNIT_INTERVAL_WIGGLE_START < all_roots) &
        (all_roots < _UNIT_INTERVAL_WIGGLE_END)]
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

    See: http://dx.doi.org/10.1016/0024-3795(70)90023-6

    .. note::

       This assumes that :math:`f \neq 0`.

    Args:
        coeffs (numpy.ndarray): ``d + 1``-array of coefficients in monomial /
            power basis.

    Raises:
        NotImplementedError: If the polynomial has non-simple roots.
    """
    coeffs = _strip_leading_zeros(coeffs)
    num_coeffs, = coeffs.shape
    if num_coeffs < 3:
        return

    deriv_poly = polynomial.polyder(coeffs)

    companion = polynomial.polycompanion(deriv_poly)
    # Use Horner's method to evaluate f(companion)
    num_companion, _ = companion.shape
    id_mat = np.eye(num_companion)
    evaluated = coeffs[-1] * id_mat
    for index in six.moves.xrange(num_coeffs - 2, -1, -1):
        coeff = coeffs[index]
        evaluated = evaluated.dot(companion) + coeff * id_mat

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


def _near_zero(value, threshold=_PARAM_THRESHOLD):
    """Rounds a number to zero if it is within a specified threshold.

    Args:
        value (float): Value to be rounded.
        threshold (Optional[float]): The threshold for rounding.

    Returns:
        float: The original value or ``0.0``.
    """
    if np.abs(value) < threshold:
        return 0.0
    else:
        return value


def _resolve_and_add(nodes1, s_val, final_s, nodes2, t_val, final_t):
    """Resolve a computed intersection and add to lists.

    We perform one Newton step to deal with any residual issues of
    high-degree polynomial solves (one of which depends the already
    approximate ``x_val, y_val``).

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        s_val (float): The approximate intersection parameter
            along ``nodes1``.
        final_s (List[float, ...]): The list of accepted intersection
            parameters ``s``.
        nodes2 (numpy.ndarray): The nodes in the second curve.
        t_val (float): The approximate intersection parameter
            along ``nodes2``.
        final_t (List[float, ...]): The list of accepted intersection
            parameters ``t``.
    """
    s_val, t_val = _intersection_helpers.newton_refine(
        s_val, nodes1, t_val, nodes2)

    s_val = _near_zero(s_val)
    t_val = _near_zero(t_val)
    if s_val < 0.0 or t_val < 0.0:
        return

    final_s.append(s_val)
    final_t.append(t_val)


def intersect_curves(nodes1, nodes2):
    r"""Intersect two parametric B |eacute| zier curves.

    Args:
        nodes1 (numpy.ndarray): The nodes in the first curve.
        nodes2 (numpy.ndarray): The nodes in the second curve.

    Returns:
        numpy.ndarray: ``Nx2`` array of intersection parameters.
        Each row contains a pair of values :math:`s` and :math:`t`
        (each in :math:`\left[0, 1\right]`) such that the curves
        intersect: :math:`B_1(s) = B_2(t)`.

    Raises:
        NotImplementedError: If the "intersection polynomial" is
        all zeros -- which indicates coincident curves.
    """
    nodes1 = _curve_helpers.full_reduce(nodes1)
    nodes2 = _curve_helpers.full_reduce(nodes2)

    num_nodes1, _ = nodes1.shape
    num_nodes2, _ = nodes2.shape
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
        (x_val, y_val), = _curve_helpers.evaluate_multi(
            nodes2, np.asfortranarray([t_val]))
        s_val = locate_point(nodes1, x_val, y_val)
        if s_val is not None:
            _resolve_and_add(
                nodes1, s_val, final_s, nodes2, t_val, final_t)

    result = np.zeros((len(final_s), 2), order='F')
    if swapped:
        final_s, final_t = final_t, final_s

    result[:, 0] = final_s
    result[:, 1] = final_t

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
        NotImplementedError: If the degree of the curve is not among
            0, 1, 2 or 3.
    """
    num_coeffs, = bezier_coeffs.shape
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
        return np.asfortranarray([
            coeff0, 2.0 * (coeff1 - coeff0), coeff2 - 2.0 * coeff1 + coeff0])
    elif num_coeffs == 4:
        #   C0 (1 - s)^3 + C1 3 (1 - s)^2 + C2 3 (1 - s) s^2 + C3 s^3
        # = C0 + 3(C1 - C0) s + 3(C2 - 2 C1 + C0) s^2 +
        #   (C3 - 3 C2 + 3 C1 - C0) s^3
        coeff0, coeff1, coeff2, coeff3 = bezier_coeffs
        return np.asfortranarray([
            coeff0,
            3.0 * (coeff1 - coeff0),
            3.0 * (coeff2 - 2.0 * coeff1 + coeff0),
            coeff3 - 3.0 * coeff2 + 3.0 * coeff1 - coeff0,
        ])
    else:
        raise NotImplementedError(
            'Currently only supports degrees 0, 1, 2 and 3')


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
    zero1 = _curve_helpers.full_reduce(nodes[:, [0]]) - x_val
    zero2 = _curve_helpers.full_reduce(nodes[:, [1]]) - y_val

    # Make sure we have the lowest degree in front, to make the polynomial
    # solve have the fewest number of roots.
    if zero1.shape[0] > zero2.shape[0]:
        zero1, zero2 = zero2, zero1

    # If the "smallest" is a constant, we can't find any roots from it.
    if zero1.shape[0] == 1:
        # NOTE: We assume that callers won't pass ``nodes`` that are
        #       degree 0, so if ``zero1`` is a constant, ``zero2`` won't be.
        zero1, zero2 = zero2, zero1

    power_basis1 = poly_to_power_basis(zero1[:, 0])
    all_roots = roots_in_unit_interval(power_basis1)
    if all_roots.size == 0:
        return

    # NOTE: We normalize ``power_basis2`` because we want to check for
    #       "zero" values, i.e. f2(s) == 0.
    power_basis2 = normalize_polynomial(poly_to_power_basis(zero2[:, 0]))

    near_zero = np.abs(polynomial.polyval(all_roots, power_basis2))
    index = np.argmin(near_zero)

    if near_zero[index] < _ZERO_THRESHOLD:
        return all_roots[index]
