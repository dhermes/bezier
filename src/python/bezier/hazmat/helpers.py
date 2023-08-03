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

"""Pure Python generic geometry and floating point helpers."""

import bisect

import numpy as np


_EPS = 0.5**40


def vector_close(vec1, vec2, eps=_EPS):
    r"""Checks that two vectors are equal to some threshold.

    Does so by computing :math:`s_1 = \|v_1\|_2` and
    :math:`s_2 = \|v_2\|_2` and then checking if

    .. math::

       \|v_1 - v_2\|_2 \leq \varepsilon \min(s_1, s_2)

    where :math:`\varepsilon = 2^{-40} \approx 10^{-12}` is a fixed
    threshold. In the rare case that one of ``vec1`` or ``vec2`` is
    the zero vector (i.e. when :math:`\min(s_1, s_2) = 0`) instead
    checks that the other vector is close enough to zero:

    .. math::

       \|v_1\|_2 = 0 \Longrightarrow \|v_2\|_2 \leq \varepsilon

    .. note::

       This function assumes that both vectors have finite values,
       i.e. that no NaN or infinite numbers occur. NumPy provides
       :func:`numpy.allclose` for coverage of **all** cases.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        vec1 (numpy.ndarray): First vector (1D) for comparison.
        vec2 (numpy.ndarray): Second vector (1D) for comparison.
        eps (float): Error threshold. Defaults to :math:`2^{-40}`.

    Returns:
        bool: Flag indicating if they are close to precision.
    """
    # NOTE: This relies on ``vec1`` and ``vec2`` being one-dimensional
    #       vectors so NumPy doesn't try to use a matrix norm.
    size1 = np.linalg.norm(vec1, ord=2)
    size2 = np.linalg.norm(vec2, ord=2)
    if size1 == 0:
        return size2 <= eps

    elif size2 == 0:
        return size1 <= eps

    else:
        upper_bound = eps * min(size1, size2)
        return np.linalg.norm(vec1 - vec2, ord=2) <= upper_bound


def in_interval(value, start, end):
    """Checks if a ``value`` is an interval (inclusive).

    .. note::

       The current implementation does the most basic check,
       however, in the future, a more generic check may be desired
       that allows wiggle room around the endpoints to account
       for round-off.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        value (float): The value to check.
        start (float): The (inclusive) start of the interval.
        end (float): The (inclusive) end of the interval.

    Returns:
        bool: Indicating if the value is in the interval.
    """
    return start <= value <= end


def bbox(nodes):
    """Get the bounding box for set of points.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
       nodes (numpy.ndarray): A set of points.

    Returns:
        Tuple[float, float, float, float]: The left, right,
        bottom and top bounds for the box.
    """
    left, bottom = np.min(nodes, axis=1)
    right, top = np.max(nodes, axis=1)
    return left, right, bottom, top


def contains_nd(nodes, point):
    r"""Predicate indicating if a point is within a bounding box.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
       nodes (numpy.ndarray): A set of points.
       point (numpy.ndarray): A 1D NumPy array representing a point
           in the same dimension as ``nodes``.

    Returns:
        bool: Indicating containment.
    """
    min_vals = np.min(nodes, axis=1)
    if not np.all(min_vals <= point):
        return False

    max_vals = np.max(nodes, axis=1)
    if not np.all(point <= max_vals):
        return False

    return True


def cross_product(vec0, vec1):
    r"""Compute the cross product of vectors in :math:`\mathbf{R}^2`.

    Utilizes the fact that

    .. math::

       \left[\begin{array}{c} A \\ B \\ 0 \end{array}\right] \times
           \left[\begin{array}{c} C \\ D \\ 0 \end{array}\right] =
           \left[\begin{array}{c} 0 \\ 0 \\ AD - BC \end{array}\right]

    and just returns the :math:`z` component.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        vec0 (numpy.ndarray): A vector as a 1D NumPy array with two values.
        vec1 (numpy.ndarray): A vector as a 1D NumPy array with two values.

    Returns:
        float: The cross product (or rather, its :math:`z` component).
    """
    return vec0[0] * vec1[1] - vec0[1] * vec1[0]


def matrix_product(mat1, mat2):
    """Compute the product of two Fortran contiguous matrices.

    This is to avoid the overhead of NumPy converting to C-contiguous
    before computing a matrix product.

    Does so via ``A B = (B^T A^T)^T`` since ``B^T`` and ``A^T`` will be
    C-contiguous without a copy, then the product ``P = B^T A^T`` will
    be C-contiguous and we can return the view ``P^T`` without a copy.

    Args:
        mat1 (numpy.ndarray): The left-hand side matrix.
        mat2 (numpy.ndarray): The right-hand side matrix.

    Returns:
        numpy.ndarray: The product of the two matrices.
    """
    return np.dot(mat2.T, mat1.T).T  # pylint: disable=no-member


def wiggle_interval(value, wiggle=0.5**44):
    r"""Check if ``value`` is in :math:`\left[0, 1\right]`.

    Allows a little bit of wiggle room outside the interval. A value
    within ``wiggle`` of ``0.0`` will be converted to ``0.0`` and similar
    for ``1.0``.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        value (float): Value to check in interval.
        wiggle (Optional[float]): The amount of wiggle room around the
            the endpoints ``0.0`` and ``1.0``. Defaults to :math:`2^{-44}`.

    Returns:
        Tuple[float, bool]: Pair of

        * The ``value`` if it's in the interval, or ``0.0`` or ``1.0``
          if the value lies slightly outside. If the ``value`` is
          too far outside the unit interval, will be NaN.
        * Boolean indicating if the ``value`` is inside the unit interval.
    """
    if -wiggle < value < wiggle:
        return 0.0, True

    elif wiggle <= value <= 1.0 - wiggle:
        return value, True

    elif 1.0 - wiggle < value < 1.0 + wiggle:
        return 1.0, True

    else:
        return np.nan, False


def cross_product_compare(start, candidate1, candidate2):
    """Compare two relative changes by their cross-product.

    This is meant to be a way to determine which vector is more "inside"
    relative to ``start``.

    .. note::

       This is a helper for :func:`simple_convex_hull`.

    Args:
        start (numpy.ndarray): The start vector (as 1D NumPy array with
            2 elements).
        candidate1 (numpy.ndarray): The first candidate vector (as 1D
            NumPy array with 2 elements).
        candidate2 (numpy.ndarray): The second candidate vector (as 1D
            NumPy array with 2 elements).

    Returns:
        float: The cross product of the two differences.
    """
    delta1 = candidate1 - start
    delta2 = candidate2 - start
    return cross_product(delta1, delta2)


def in_sorted(values, value):
    """Checks if a value is in a sorted list.

    Uses the :mod:`bisect` builtin to find the insertion point for
    ``value``.

    Args:
        values (List[int]): Integers sorted in ascending order.
        value (int): Value to check if contained in ``values``.

    Returns:
        bool: Indicating if the value is contained.
    """
    index = bisect.bisect_left(values, value)
    if index >= len(values):
        return False

    return values[index] == value


def simple_convex_hull(points):
    r"""Compute the convex hull for a set of points.

    .. _wikibooks: https://en.wikibooks.org/wiki/Algorithm_Implementation/\
                   Geometry/Convex_hull/Monotone_chain

    This uses Andrew's monotone chain convex hull algorithm and this code
    used a `wikibooks`_ implementation as motivation. The code there
    is licensed CC BY-SA 3.0.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built. Note that ``scipy.spatial.ConvexHull``
       can do this as well (via Qhull), but that would require a hard
       dependency on ``scipy`` and that helper computes much more than we need.

    .. note::

       This computes the convex hull in a "naive" way. It's expected that
       internal callers of this function will have a small number of points
       so ``n log n`` vs. ``n^2`` vs. ``n`` aren't that relevant.

    Args:
        points (numpy.ndarray): A ``2 x N`` array (``float64``) of points.

    Returns:
        numpy.ndarray: The ``2 x N`` array (``float64``) of ordered points in
        the polygonal convex hull.
    """
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-branches
    if points.size == 0:
        return points

    # First, drop duplicates.
    unique_points = np.unique(points, axis=1)
    _, num_points = unique_points.shape
    if num_points < 2:
        return unique_points

    # Then sort the data in left-to-right order (and break ties by y-value).
    points = np.empty((2, num_points), order="F")
    for index, xy_val in enumerate(
        sorted(tuple(column) for column in unique_points.T)
    ):
        points[:, index] = xy_val
    # After sorting, if there are only 2 points, return.
    if num_points < 3:
        return points

    # Build lower hull
    lower = [0, 1]
    for index in range(2, num_points):
        point2 = points[:, index]
        while len(lower) >= 2:
            point0 = points[:, lower[-2]]
            point1 = points[:, lower[-1]]
            if cross_product_compare(point0, point1, point2) > 0:
                break

            lower.pop()

        lower.append(index)
    # Build upper hull
    upper = [num_points - 1]
    for index in range(num_points - 2, -1, -1):
        # Don't consider indices from the lower hull (other than the ends).
        if index > 0 and in_sorted(lower, index):
            continue

        point2 = points[:, index]
        while len(upper) >= 2:
            point0 = points[:, upper[-2]]
            point1 = points[:, upper[-1]]
            if cross_product_compare(point0, point1, point2) > 0:
                break

            upper.pop()

        upper.append(index)
    # **Both** corners are double counted.
    size_polygon = len(lower) + len(upper) - 2
    polygon = np.empty((2, size_polygon), order="F")
    for index, column in enumerate(lower[:-1]):
        polygon[:, index] = points[:, column]
    index_start = len(lower) - 1
    for index, column in enumerate(upper[:-1]):
        polygon[:, index + index_start] = points[:, column]
    return polygon


def is_separating(direction, polygon1, polygon2):
    """Checks if a given ``direction`` is a separating line for two polygons.

    .. note::

       This is a helper for :func:`polygon_collide`.

    Args:
        direction (numpy.ndarray): A 1D ``2``-array (``float64``) of a
            potential separating line for the two polygons.
        polygon1 (numpy.ndarray): A ``2 x N`` array (``float64``) of ordered
            points in a polygon.
        polygon2 (numpy.ndarray): A ``2 x N`` array (``float64``) of ordered
            points in a polygon.

    Returns:
        bool: Flag indicating if ``direction`` is a separating line.
    """
    # NOTE: We assume throughout that ``norm_squared != 0``. If it **were**
    #       zero that would mean the ``direction`` corresponds to an
    #       invalid edge.
    norm_squared = direction[0] * direction[0] + direction[1] * direction[1]
    params = []
    vertex = np.empty((2,), order="F")
    for polygon in (polygon1, polygon2):
        _, polygon_size = polygon.shape
        min_param = np.inf
        max_param = -np.inf
        for index in range(polygon_size):
            vertex[:] = polygon[:, index]
            param = cross_product(direction, vertex) / norm_squared
            min_param = min(min_param, param)
            max_param = max(max_param, param)
        params.append((min_param, max_param))
    # NOTE: The indexing is based on:
    #       params[0] = (min_param1, max_param1)
    #       params[1] = (min_param2, max_param2)
    return params[0][0] > params[1][1] or params[0][1] < params[1][0]


def polygon_collide(polygon1, polygon2):
    """Determines if two **convex** polygons collide.

    .. _SAT: https://en.wikipedia.org/wiki/Hyperplane_separation_theorem
    .. _see also: https://hackmd.io/s/ryFmIZrsl

    This code uses the Separating axis theorem (`SAT`_) to quickly
    determine if the polygons intersect. `See also`_.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        polygon1 (numpy.ndarray): A ``2 x N`` array (``float64``) of ordered
            points in a polygon.
        polygon2 (numpy.ndarray): A ``2 x N`` array (``float64``) of ordered
            points in a polygon.

    Returns:
        bool: Flag indicating if the two polygons collide.
    """
    direction = np.empty((2,), order="F")
    for polygon in (polygon1, polygon2):
        _, polygon_size = polygon.shape
        for index in range(polygon_size):
            # NOTE: When ``index == 0`` this will "wrap around" and refer
            #       to index ``-1``.
            direction[:] = polygon[:, index] - polygon[:, index - 1]
            if is_separating(direction, polygon1, polygon2):
                return False

    return True


def solve2x2(lhs, rhs):
    """Solve a square 2 x 2 system via LU factorization.

    This is meant to be a stand-in for LAPACK's ``dgesv``, which just wraps
    two calls to ``dgetrf`` and ``dgetrs``. We wrap for two reasons:

    * We seek to avoid exceptions as part of the control flow (which is
      what :func:`numpy.linalg.solve` does).
    * We seek to avoid excessive type- and size-checking, since this
      special case is already known.

    Args:
        lhs (numpy.ndarray): A ``2 x 2`` array of real numbers.
        rhs (numpy.ndarray): A 1D array of 2 real numbers.

    Returns:
        Tuple[bool, float, float]: A triple of

        * A flag indicating if ``lhs`` is a singular matrix.
        * The first component of the solution.
        * The second component of the solution.
    """
    # A <--> lhs[0, 0]
    # B <--> lhs[0, 1]
    # C <--> lhs[1, 0]
    # D <--> lhs[1, 1]
    # E <--> rhs[0]
    # F <--> rhs[1]
    if np.abs(lhs[1, 0]) > np.abs(lhs[0, 0]):
        # NOTE: We know there is no division by zero here since ``C``
        #       is **strictly** bigger than **some** value (in magnitude).
        # [A | B][x] = [E]
        # [C | D][y]   [F]
        ratio = lhs[0, 0] / lhs[1, 0]
        # r = A / C
        # [A - rC | B - rD][x]   [E - rF]
        # [C      | D     ][y] = [F     ]
        # ==> 0x + (B - rD) y = E - rF
        denominator = lhs[0, 1] - ratio * lhs[1, 1]
        if denominator == 0.0:
            return True, None, None

        y_val = (rhs[0] - ratio * rhs[1]) / denominator
        # Cx + Dy = F ==> x = (F - Dy) / C
        x_val = (rhs[1] - lhs[1, 1] * y_val) / lhs[1, 0]
        return False, x_val, y_val

    else:
        if lhs[0, 0] == 0.0:
            return True, None, None

        # [A | B][x] = [E]
        # [C | D][y]   [F]
        ratio = lhs[1, 0] / lhs[0, 0]
        # r = C / A
        # [A      | B     ][x] = [E     ]
        # [C - rA | D - rB][y]   [F - rE]
        # ==> 0x + (D - rB) y = F - rE
        denominator = lhs[1, 1] - ratio * lhs[0, 1]
        if denominator == 0.0:
            return True, None, None

        y_val = (rhs[1] - ratio * rhs[0]) / denominator
        # Ax + By = E ==> x = (E - B y) / A
        x_val = (rhs[0] - lhs[0, 1] * y_val) / lhs[0, 0]
        return False, x_val, y_val


class UnsupportedDegree(NotImplementedError):
    """Custom exception to indicate the given degree is unsupported.

    This is intentionally a subclass of a :exc:`NotImplementedError`
    since it's intended to indicate a lack of an implementation. For
    example, :meth:`.Curve.reduce_` uses hard-coded matrices for
    a small subset of possible degrees, so the implementation is
    **degree-specific**:

    .. doctest:: unsupported-degree
       :options: +NORMALIZE_WHITESPACE

       >>> import bezier
       >>> import numpy as np
       >>> degree = 5
       >>> nodes = np.empty((2, degree + 1), order="F")
       >>> curve = bezier.Curve(nodes, degree=degree)
       >>> curve.reduce_()
       Traceback (most recent call last):
         ...
       bezier.hazmat.helpers.UnsupportedDegree: The only degrees supported at
                                                this time are 1, 2, 3 and
                                                4 (degree=5)

    Args:
        degree (int): The degree that is not possible to support.
        supported (Tuple[int, ...]): The degrees that are
            actually supported by the failing method.
    """

    def __init__(self, degree, supported=()):
        super().__init__()
        self.degree = degree
        """int: The degree that the caller attempted to use."""
        self.supported = supported
        """Tuple[int, ...]: The degrees supported by the failing method."""

    def __str__(self):
        num_supported = len(self.supported)
        if num_supported == 0:
            return f"degree={self.degree}"

        degrees_str = [str(degree) for degree in self.supported]
        if num_supported == 1:
            msg = "The only degree supported at this time is " + degrees_str[0]
        else:
            msg = (
                "The only degrees supported at this time are "
                + ", ".join(degrees_str[:-1])
                + " and "
                + degrees_str[-1]
            )
        return f"{msg} (degree={self.degree})"
