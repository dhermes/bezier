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

"""Common features of B |eacute| zier shapes.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import numpy as np


class Base:
    """Base shape object.

    Args:
        nodes (Sequence[Sequence[numbers.Number]]): The control points for the
            shape. Must be convertible to a 2D NumPy array of floating point
            values, where the columns are the nodes and the rows correspond to
            each dimension the shape occurs in.
        copy (bool): Flag indicating if the nodes should be copied before
            being stored. Defaults to :data:`True` since callers may
            freely mutate ``nodes`` after passing in.

    Raises:
        ValueError: If the ``nodes`` are not 2D.
    """

    __slots__ = ("_dimension", "_nodes")
    _degree = -1

    def __init__(self, nodes, copy=True):
        nodes_np = sequence_to_array(nodes)
        dimension, _ = nodes_np.shape
        self._dimension = dimension
        if copy:
            self._nodes = nodes_np.copy(order="F")
        else:
            self._nodes = nodes_np

    @property
    def degree(self):
        """int: The degree of the current shape."""
        return self._degree

    @property
    def dimension(self):
        r"""int: The dimension that the shape lives in.

        For example, if the shape lives in :math:`\mathbf{R}^3`, then
        the dimension is ``3``.
        """
        return self._dimension

    @property
    def nodes(self):
        """numpy.ndarray: The nodes that define the current shape."""
        return self._nodes.copy(order="F")

    def __repr__(self):
        """Representation of current object.

        Returns:
            str: Object representation.
        """
        return (
            f"<{self.__class__.__name__} "
            f"(degree={self._degree:d}, dimension={self._dimension:d})>"
        )


def _lossless_to_float(array):
    """Convert a NumPy array to ``np.float64`` data type.

    Args:
        array (numpy.ndarray): The NumPy array to convert to float.

    Returns:
        numpy.ndarray: The converted array (or the original if already a
        float).

    Raises:
        ValueError: If ``array`` can't be directly converted without rounding.
    """
    if array.dtype == np.float64:
        return array

    converted = array.astype(np.float64)
    if not np.all(array == converted):
        raise ValueError("Array cannot be converted to floating point")

    return converted


def sequence_to_array(nodes):
    """Convert a sequence to a Fortran-ordered ``np.float64`` NumPy array.

    Args:
        nodes (Sequence[Sequence[numbers.Number]]): The control points for a
            shape. Must be convertible to a 2D NumPy array of floating point
            values, where the columns are the nodes and the rows correspond to
            each dimension the shape occurs in.

    Returns:
        numpy.ndarray: The converted array (or the original if already a
        float array).

    Raises:
        ValueError: If the ``nodes`` are not 2D.
    """
    nodes_np = np.asarray(nodes, order="F")
    if nodes_np.ndim != 2:
        raise ValueError("Nodes must be 2-dimensional, not", nodes_np.ndim)

    return _lossless_to_float(nodes_np)
