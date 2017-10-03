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

"""Generic geometry and floating point helpers.

As a convention, the functions defined here with a leading underscore
(e.g. :func:`_bbox`) have a special meaning.

Each of these functions have a Cython speedup with the exact same
interface which calls out to a Fortran implementation. The speedup
will be used if the extension can be built. The name **without** the
leading underscore will be surfaced as the actual interface (e.g.
``bbox``) whether that is the pure Python implementation or the speedup.
"""


import numpy as np

try:
    from bezier import _helpers_speedup
except ImportError:  # pragma: NO COVER
    _helpers_speedup = None


_EPS = 0.5**40


def _vector_close(vec1, vec2, eps=_EPS):
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
       :func:`np.allclose` for coverage of **all** cases.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        vec1 (numpy.ndarray): First vector for comparison.
        vec2 (numpy.ndarray): Second vector for comparison.
        eps (float): Error threshold. Defaults to :math:`2^{-40}`.

    Returns:
        bool: Flag indicating if they are close to precision.
    """
    # NOTE: We assume the caller sends a 1xD vector. We turn it into
    #       a one-dimensional vector so NumPy doesn't use a matrix norm.
    size1 = np.linalg.norm(vec1[0, :], ord=2)
    size2 = np.linalg.norm(vec2[0, :], ord=2)
    if size1 == 0:
        return size2 <= eps
    elif size2 == 0:
        return size1 <= eps
    else:
        upper_bound = eps * min(size1, size2)
        return np.linalg.norm(vec1[0, :] - vec2[0, :], ord=2) <= upper_bound


def _in_interval(value, start, end):
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


def _bbox(nodes):
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
    left, bottom = np.min(nodes, axis=0)
    right, top = np.max(nodes, axis=0)
    return left, right, bottom, top


def _contains_nd(nodes, point):
    r"""Predicate indicating if a point is within a bounding box.

    Like :func:`contains` but supports points in arbitrary dimension.
    Unlike :func:`contains`, this function directly uses ``<=`` and
    ``>=`` for comparison (:func:`contains` uses :func:`in_interval`).

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
    min_vals = np.min(nodes, axis=0)
    if not np.all(min_vals <= point):
        return False

    max_vals = np.max(nodes, axis=0)
    if not np.all(point <= max_vals):
        return False

    return True


def _cross_product(vec0, vec1):
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
        vec0 (numpy.ndarray): A vector as a 1x2 NumPy array.
        vec1 (numpy.ndarray): A vector as a 1x2 NumPy array.

    Returns:
        float: The cross product (or rather, its :math:`z` component).
    """
    return vec0[0, 0] * vec1[0, 1] - vec0[0, 1] * vec1[0, 0]


def _ulps_away(value1, value2, num_bits=1):
    r"""Determines if ``value1`` is within ``n`` ULPs of ``value2``.

    Uses ``np.spacing`` to determine the unit of least precision (ULP)
    for ``value1`` and then checks that the different between the values
    does not exceed ``n`` ULPs.

    When ``value1 == 0`` or ``value2 == 0``, we instead check that the other
    is less than :math:`2^{-40}` (``_EPS``) in magnitude.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        value1 (float): The first value that being compared.
        value2 (float): The second value that being compared.
        num_bits (Optional[int]): The number of bits allowed to differ.
            Defaults to ``1``.

    Returns:
        bool: Predicate indicating if the values agree to ``n`` bits.
    """
    if value1 == 0.0:
        return abs(value2) < _EPS
    elif value2 == 0.0:
        return abs(value1) < _EPS
    else:
        local_epsilon = np.spacing(value1)  # pylint: disable=no-member
        return abs(value1 - value2) <= num_bits * abs(local_epsilon)


def eye(num_elts):
    """Get ``n x n`` identity matrix.

    Provided here because ``np.eye`` doesn't allow creation in Fortran order.

    Args:
        num_elts (int): The number of elements (``n``) on each side.

    Returns:
        numpy.ndarray: The identity matrix.
    """
    id_mat = np.zeros((num_elts, num_elts), order='F')
    id_mat.flat[::num_elts + 1] = 1.0
    return id_mat


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
    return np.dot(mat2.T, mat1.T).T


def _wiggle_interval(value, wiggle=0.5**45):
    r"""Check if ``value`` is in :math:`\left[0, 1\right]`.

    Allows a little bit of wiggle room outside the interval. Any value
    within ``wiggle`` of ``0.0` will be converted to ``0.0` and similar
    for ``1.0``.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Args:
        value (float): Value to check in interval.
        wiggle (Optional[float]): The amount of wiggle room around the
            the endpoints ``0.0`` and ``1.0``.

    Returns:
        Tuple[float, bool]: Pair of

        * The ``value`` if it's in the interval, or ``0`` or ``1``
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


# pylint: disable=invalid-name
if _helpers_speedup is None:  # pragma: NO COVER
    vector_close = _vector_close
    in_interval = _in_interval
    bbox = _bbox
    contains_nd = _contains_nd
    cross_product = _cross_product
    ulps_away = _ulps_away
    wiggle_interval = _wiggle_interval
else:
    vector_close = _helpers_speedup.vector_close
    in_interval = _helpers_speedup.in_interval
    bbox = _helpers_speedup.bbox
    contains_nd = _helpers_speedup.contains_nd
    cross_product = _helpers_speedup.cross_product
    ulps_away = _helpers_speedup.ulps_away
    wiggle_interval = _helpers_speedup.wiggle_interval
# pylint: enable=invalid-name
