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

"""Generic geometry and floating point helpers."""


import numpy as np


_EPS = 2.0**(-40)


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
       :func:`np.allclose` for coverage of **all** cases.

    Args:
        vec1 (numpy.ndarray): First vector for comparison.
        vec2 (numpy.ndarray): Second vector for comparison.
        eps (float): Error threshold. Defaults to :math:`2^{-40}`.

    Returns:
        bool: Flag indicating if they are close to precision.
    """
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

    Args:
       nodes (numpy.ndarray): A set of points.

    Returns:
        Tuple[float, float, float, float]: The left, right,
        bottom and top bounds for the box.
    """
    left, bottom = np.min(nodes, axis=0)
    right, top = np.max(nodes, axis=0)
    return left, right, bottom, top


def contains(nodes, x_val, y_val):
    r"""Predicate indicating if a point is within a bounding box.

    Args:
       nodes (numpy.ndarray): A set of points.
       x_val (float): The :math:`x`-coordinate of the point.
       y_val (float): The :math:`y`-coordinate of the point.

    Returns:
        bool: Indicating containment.
    """
    left, right, bottom, top = bbox(nodes)
    if not in_interval(x_val, left, right):
        return False
    if not in_interval(y_val, bottom, top):
        return False
    return True
