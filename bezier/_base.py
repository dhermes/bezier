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

"""Common features of B |eacute| zier shapes.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


class Base(object):
    """Base shape object.

    Args:
        nodes (numpy.ndarray): The control points for the shape.
            Must be a 2D array, where the rows are the nodes and the
            columns correspond to each dimension the shape occurs in.
        _copy (bool): Flag indicating if the nodes should be copied before
            being stored. Defaults to :data:`True` since callers may
            freely mutate ``nodes`` after passing in.

    Raises:
        ValueError: If the ``nodes`` are not 2D.
        ValueError: If the ``degree`` is less than ``1``.
    """

    __slots__ = ('_degree', '_dimension', '_nodes')

    def __init__(self, nodes, _copy=True):
        if nodes.ndim != 2:
            raise ValueError('Nodes must be 2-dimensional, not',
                             nodes.ndim)

        num_nodes, dimension = nodes.shape
        degree = self._get_degree(num_nodes)
        if degree < 1:
            raise ValueError('Shape must be at least degree 1')

        self._degree = degree
        self._dimension = dimension
        if _copy:
            self._nodes = nodes.copy()
        else:
            self._nodes = nodes

    @staticmethod
    def _get_degree(num_nodes):
        """Get the degree of the current shape.

        .. note::

           Subclasses are expected to implement this method.

        Args:
            num_nodes (int): The number of nodes provided.

        Returns:
            int: The degree of the current shape.
        """
        return num_nodes

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
        return self._nodes.copy()

    def __repr__(self):
        """Representation of current object.

        Returns:
            str: Object representation.
        """
        return '<{} (degree={:d}, dimension={:d})>'.format(
            self.__class__.__name__, self._degree, self._dimension)
