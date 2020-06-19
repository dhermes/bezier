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

"""This implements a test for closed loops.

.. _clipping: https://dx.doi.org/10.1016/0010-4485(90)90039-F
.. _note: https://doi.org/10.1016/0010-4485(89)90058-4

The collinear normal theorem from the B |eacute| zier `clipping`_ paper
(from Sederberg and Nishita) refers to a short `note`_ from
Sederberg, Christiansen, Katz.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import numpy as np

from bezier import _curve_helpers
from bezier import _helpers
from bezier.hazmat import geometric_intersection


DEFAULT_S_MIN = 1.0
DEFAULT_S_MAX = 0.0
BINOMIAL2 = np.asfortranarray([1.0, 2.0, 1.0])
BINOMIAL3 = np.asfortranarray([1.0, 3.0, 3.0, 1.0])
BINOMIAL4 = np.asfortranarray([1.0, 4.0, 6.0, 4.0, 1.0])
BINOMIAL5 = np.asfortranarray([1.0, 5.0, 10.0, 10.0, 5.0, 1.0])
BINOMIAL6 = np.asfortranarray([1.0, 6.0, 15.0, 20.0, 15.0, 6.0, 1.0])
BINOMIAL7 = np.asfortranarray([1.0, 7.0, 21.0, 35.0, 35.0, 21.0, 7.0, 1.0])
ROT90 = np.asfortranarray([[0.0, -1.0], [1.0, 0.0]])


def binomial(n, k):
    if n < 2:
        return 1.0

    elif n == 2:
        return BINOMIAL2[k]

    elif n == 3:
        return BINOMIAL3[k]

    elif n == 4:
        return BINOMIAL4[k]

    elif n == 5:
        return BINOMIAL5[k]

    elif n == 6:
        return BINOMIAL6[k]

    elif n == 7:
        return BINOMIAL7[k]

    else:
        raise NotImplementedError(n)


def compute_focus(nodes):
    r"""Compute the focus curve for a given B |eacute| zier curve.

    For a degree :math:`n` curve :math:`B(t)`, the focus is of the form

    .. math::

       F(t) = B(t) + c(t) N(t)

    where :math:`N(t)` is the hodograph rotated (and scaled) so that it
    is normal to :math:`B(t)`.

    .. math::

       R_{\pi/2} = \left[\begin{array}{c c}
           0 & -1 \\ 1 & 0 \end{array}\right] \qquad
       N(t) = \frac{1}{n} R_{\pi/2} B'(t)

    and :math:`c(t) = c_0 (1 - t) + c_1 t` is a line chosen so that
    :math:`F(0) = F(1)`. We write the control points of the normal curve
    :math:`N_j = R_{\pi/2} \left(P_{j + 1} - P_j\right) =
    R_{\pi/2} \Delta P_j` and then the conditions on :math:`c_0, c_1` become

    .. math::

       P_0 + c_0 N_0 = F(0) = F(1) = P_n + c_1 N_{n - 1} \\
       \Longleftrightarrow c_0 N_0 - c_1 N_{n - 1} = P_n - P_0 \\
       \Longleftrightarrow R_{\pi/2} \left[\begin{array}{c | c}
          \Delta P_0 & \Delta P_{n - 1} \end{array}\right]
       \left[\begin{array}{c c}
           1 & 0 \\ 0 & -1 \end{array}\right]
       \left[\begin{array}{c}
           c_0 \\ c_1 \end{array}\right] = P_n - P_0.

    Once :math:`c_0` and :math:`c_1` have been computed, the control
    points of :math:`F(t)` can be shown to be

    .. math::

       \begin{align*}
       F_0 &= P_0 + c_0 N_0 \\
       F_j &= P_j + \frac{n - j}{n} c_0 N_j + \frac{j}{n} c_1 N_{j - 1} \\
       F_n &= P_n + c_1 N_{n - 1}
       \end{align*}
    """
    # NOTE: This assumes, but does not check, that ``num_nodes >= 2``.
    _, num_nodes = nodes.shape
    degree = num_nodes - 1

    lhs_mat = np.empty((2, 2), order="F")
    N_curr = ROT90.dot(nodes[:, 1] - nodes[:, 0])
    lhs_mat[:, 0] = N_curr
    lhs_mat[:, 1] = ROT90.dot(nodes[:, degree] - nodes[:, degree - 1])
    rhs_mat = nodes[:, degree] - nodes[:, 0]

    # Solve the system.
    c0, c1 = np.linalg.solve(lhs_mat, rhs_mat)
    c1 = -c1
    # TODO: Bound c0 and c1.

    focus_nodes = nodes.copy(order="F")
    # Special handling for first column (N{-1} isn't defined).
    focus_nodes[:, 0] += c0 * N_curr

    denominator = float(degree)
    for index in range(1, degree):
        N_prev = N_curr
        N_curr = ROT90.dot(nodes[:, index + 1] - nodes[:, index])
        focus_nodes[:, index] += (
            c0 * (degree - index) * N_curr + c1 * index * N_prev
        ) / denominator

    # Special handling for last column (Nn isn't defined).
    focus_nodes[:, degree] += c1 * N_curr

    return focus_nodes


def distance_patch(nodes1, nodes2):
    r"""Compute the control points of a distance function.

    .. note::

       This assumes, but does not check, that each curve is degree at least 2.

    Given two B |eacute| zier curves

    .. math::

       B_1(s) = \sum_{i = 0}^n \binom{n}{i} (1 - s)^{n - i} s^i P_i \\
       B_2(t) = \sum_{j = 0}^m \binom{m}{j} (1 - t)^{m - j} t^j Q_j

    we seek to find a collinear normal for the two curves. To do this,
    we'll try to clip the range for :math:`s` by using a focus curve
    :math:`F_2(t)` and computing the distance function:

    .. math::

       \begin{align*}
       D(s, t) &= \left(B_1(s) - F_2(t)\right) \cdot \frac{B_1'(s)}{n} \\
               &= \sum_{i = 0}^{2n - 1} \binom{2n - 1}{i}
               (1 - s)^{2n - 1 - i} s^i \left[\sum_{j = 0}^m \binom{m}{j}
               (1 - t)^{m - j} t^j d_{ij}\right]
       \end{align*}

    in Bernstein form. Once this has been done, a B |eacute| zier triangle
    patch with control points

    .. math::

       D_{ij} = \left[\begin{array}{c}
       \frac{i}{2n - 1} \\ \frac{j}{m} \\ d_{ij} \end{array}\right]

    will be used to clip the :math:`s` range.

    Specifically, these control points will be projected onto :math:`t = 0`
    and we will take the convex hull (in :math:`\mathbf{R}^2`) of

    .. math::

       \left[\begin{array}{c}
           \frac{i}{2n - 1} \\ d_{ij} \end{array}\right]

    and then intersect that convex hull with the line :math:`d = 0`. This
    will provide a segment of "accepted" :math:`s` values (i.e. those that
    lie inside the convex hull).

    If the focus is given with control points :math:`F_j`:
    :math:`F_2(t) = \sum_{j = 0}^m \binom{m}{j} (1 - t)^{m - j} t^j F_j` then
    we can write

    .. math::

       \binom{2n - 1}{i} d_{ij} = \sum_{k + \ell = i} \binom{n}{k}
           \binom{n - 1}{\ell} \left(P_k - F_j\right) \cdot \left(
           P_{\ell + 1} - P_{\ell}\right).
    """
    _, num_nodes1 = nodes1.shape
    degree1 = num_nodes1 - 1
    _, num_nodes2 = nodes2.shape
    degree2 = num_nodes2 - 1

    focus_nodes = compute_focus(nodes2)
    # Row ``j`` will contain ``d_{ij}``, so the number of
    # rows corresponds to the number of indices in ``F_2(t)``.
    # The number of columns corresponds to the number of indices
    # allowed for a degree ``2n - 1`` Bernstein basis.
    control_mesh = np.zeros((num_nodes2, 2 * degree1), order="F")

    for i in range(2 * degree1):
        # control_mesh[j, i] == d_{ij}
        min_k = max(0, 1 + i - degree1)
        max_k = min(i, degree1)
        for k in range(min_k, max_k + 1):
            ell = i - k
            multiplier = binomial(degree1, k) * binomial(degree1 - 1, ell)
            # Use ``.reshape`` to avoid a copy.
            pk = nodes1[:, k].reshape((2, 1), order="F")
            # Use broadcasting for subtraction.
            delta_f = pk - focus_nodes
            delta_p = nodes1[:, ell + 1] - nodes1[:, ell]
            # Intentionally reshape as the transpose (so that we
            # may take a dot product with it).
            delta_p = delta_p.reshape((1, 2), order="F")
            control_mesh[:, i] += (multiplier * delta_p.dot(delta_f)).ravel(
                order="A"
            )

        # After accumulating all terms in the summand, divide by
        # ((2n - 1) C i).
        control_mesh[:, i] /= binomial(2 * degree1 - 1, i)

    return control_mesh


def minmax(vector):
    """Simultaneously compute the minimum and maximum value in an array.

    .. note::

       This assumes, but does not check, that ``vector`` is 1D and has
       at least one value (i.e. isn't empty).

    Args:
        vector (numpy.ndarray): A 1D NumPy array with at least
            one element.

    Returns:
        Tuple[float, float]: Pair of the minimum and maximum values.
    """
    max_val = vector[0]
    min_val = vector[0]
    for value in vector[1:]:
        if value < min_val:
            min_val = value
        elif value > max_val:
            max_val = value

    return min_val, max_val


def clip_collinear(nodes1, nodes2):
    _, num_nodes1 = nodes1.shape
    degree1 = num_nodes1 - 1
    control_mesh = distance_patch(nodes1, nodes2)

    # Before computing the convex hull, limit to the maximum
    # and minimum ``d_{ij}`` for each fixed ``j``.
    _, size_mesh = control_mesh.shape
    hull_candidates = np.empty((2, 2 * size_mesh), order="F")

    # Rather than considering ``(j/n, d)`` we consider ``(j, nd)``.
    multiplier = float(size_mesh - 1)
    for j in range(size_mesh):
        d_min, d_max = minmax(control_mesh[:, j])
        hull_candidates[0, 2 * j] = j
        hull_candidates[1, 2 * j] = multiplier * d_min
        hull_candidates[0, 2 * j + 1] = j
        hull_candidates[1, 2 * j + 1] = multiplier * d_max

    zero_start = np.asfortranarray([0.0, 0.0])
    zero_end = np.asfortranarray([multiplier, 0.0])

    polygon = _helpers.simple_convex_hull(hull_candidates)
    _, size_polygon = polygon.shape
    s_min = DEFAULT_S_MIN
    s_max = DEFAULT_S_MAX
    for start_index in range(-1, size_polygon - 1):
        end_index = start_index + 1
        s, t, success = geometric_intersection.segment_intersection(
            zero_start,
            zero_end,
            polygon[:, start_index],
            polygon[:, end_index],
        )
        if success:
            if _helpers.in_interval(t, 0.0, 1.0):
                if _helpers.in_interval(s, 0.0, s_min):
                    s_min = s
                if _helpers.in_interval(s, s_max, 1.0):
                    s_max = s
        else:
            raise RuntimeError

    if (s_min, s_max) == (DEFAULT_S_MIN, DEFAULT_S_MAX):
        return 0.0, 1.0

    else:
        return s_min, s_max


def foo(value):
    # return value
    return np.spacing(value)


def main():
    # NP = np.asfortranarray([
    #     [536805376.0, 536854528.0, 536903680.0],
    #     [        4.0,        -2.0,         1.0],
    # ])
    # NQ = np.asfortranarray([
    #     [536739840.0, 536870912.0, 537001984.0],
    #     [      -16.0,        16.0,       -16.0],
    # ])
    # nodesP = 0.5**29 * NP
    # nodesQ = 0.5**29 * NQ

    nodesP = np.asfortranarray(
        [[0.5, 1.25, 2.0], [0.125, -0.25, 0.5]]
    )  # nodes38
    nodesQ = np.asfortranarray(
        [[0.5, 1.0, 1.5], [-0.125, 0.125, -0.125]]
    )  # nodes39

    nodes1a = nodesP
    nodes2a = nodesQ
    s1, e1 = 0.0, 1.0
    s2, e2 = 0.0, 1.0

    third = 1.0 / 3.0
    half = 1.0 / 2.0

    # Clip nodes1a with nodes2a.
    start, end = clip_collinear(nodes1a, nodes2a)
    nodes1b = _curve_helpers.specialize_curve(nodes1a, start, end)
    s1, e1 = s1 + start * (e1 - s1), s1 + end * (e1 - s1)
    print(
        (
            start,
            end,
            s1,
            e1,
            (s1 - third) / foo(third),
            (e1 - third) / foo(third),
        )
    )

    # Clip nodes2a with nodes1b.
    start, end = clip_collinear(nodes2a, nodes1b)
    nodes2b = _curve_helpers.specialize_curve(nodes2a, start, end)
    s2, e2 = s2 + start * (e2 - s2), s2 + end * (e2 - s2)
    print(
        (start, end, s2, e2, (s2 - half) / foo(half), (e2 - half) / foo(half))
    )

    # Clip nodes1b with nodes2b.
    start, end = clip_collinear(nodes1b, nodes2b)
    nodes1c = _curve_helpers.specialize_curve(nodes1b, start, end)
    s1, e1 = s1 + start * (e1 - s1), s1 + end * (e1 - s1)
    print(
        (
            start,
            end,
            s1,
            e1,
            (s1 - third) / foo(third),
            (e1 - third) / foo(third),
        )
    )

    # Clip nodes2b with nodes1c.
    start, end = clip_collinear(nodes2b, nodes1c)
    nodes2c = _curve_helpers.specialize_curve(nodes2b, start, end)
    s2, e2 = s2 + start * (e2 - s2), s2 + end * (e2 - s2)
    print(
        (start, end, s2, e2, (s2 - half) / foo(half), (e2 - half) / foo(half))
    )

    # Clip nodes1c with nodes2c.
    start, end = clip_collinear(nodes1c, nodes2c)
    nodes1d = _curve_helpers.specialize_curve(nodes1c, start, end)
    s1, e1 = s1 + start * (e1 - s1), s1 + end * (e1 - s1)
    print(
        (
            start,
            end,
            s1,
            e1,
            (s1 - third) / foo(third),
            (e1 - third) / foo(third),
        )
    )

    # Clip nodes2c with nodes1d.
    start, end = clip_collinear(nodes2c, nodes1d)
    nodes2d = _curve_helpers.specialize_curve(nodes2c, start, end)
    s2, e2 = s2 + start * (e2 - s2), s2 + end * (e2 - s2)
    print(
        (start, end, s2, e2, (s2 - half) / foo(half), (e2 - half) / foo(half))
    )

    # Clip nodes1d with nodes2d.
    start, end = clip_collinear(nodes1d, nodes2d)
    nodes1e = _curve_helpers.specialize_curve(nodes1d, start, end)
    s1, e1 = s1 + start * (e1 - s1), s1 + end * (e1 - s1)
    print(
        (
            start,
            end,
            s1,
            e1,
            (s1 - third) / foo(third),
            (e1 - third) / foo(third),
        )
    )

    # Clip nodes2d with nodes1e.
    start, end = clip_collinear(nodes2d, nodes1e)
    nodes2e = _curve_helpers.specialize_curve(nodes2d, start, end)
    s2, e2 = s2 + start * (e2 - s2), s2 + end * (e2 - s2)
    print(
        (start, end, s2, e2, (s2 - half) / foo(half), (e2 - half) / foo(half))
    )

    # See that we cannot clip nodes1e with nodes2e.
    try:
        clip_collinear(nodes1e, nodes2e)
        raise RuntimeError("Should not reach")

    except np.linalg.LinAlgError:
        pass

    return nodesP, nodesQ


if __name__ == "__main__":
    main()
