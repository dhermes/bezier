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

r"""Curved polygon and associated helpers.

A curved polygon (in :math:`\mathbf{R}^2`) is defined by the
collection of B |eacute| zier curves that determine the
boundary.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


class CurvedPolygon(object):
    """Represents an object defined by its curved boundary.

    The boundary is a piecewise defined collection of
    B |eacute| zier curves.

    .. note::

       The direction of the nodes in each :class:`.Curve`
       on the boundary is important: we check that one curve
       begins where the last one ended.

    .. doctest:: curved-polygon-ctor

       >>> import bezier
       >>> edge0 = bezier.Curve(np.array([
       ...     [0.0,  0.0],
       ...     [0.5, -1.0],
       ...     [1.0,  0.0],
       ... ]))
       >>> edge1 = bezier.Curve(np.array([
       ...     [1.0, 0.0],
       ...     [1.0, 2.0],
       ... ]))
       >>> edge2 = bezier.Curve(np.array([
       ...     [1.0, 0.0],
       ...     [0.5, 3.0],
       ...     [0.0, 0.0],
       ... ]))
       >>> edge3 = bezier.Curve(np.array([
       ...     [0.0, 2.0],
       ...     [0.0, 0.0],
       ... ]))
       >>> curved_poly = bezier.CurvedPolygon(
       ...     edge0, edge1, edge2, edge3)
       >>> curved_poly
       <CurvedPolygon (num_sides=4)>

    Args:
        edges (Tuple[.Curve, ...]): The boundary edges of the
            curved polygon.
    """

    def __init__(self, *edges):
        self._edges = edges
        self._num_sides = len(edges)
        self._verify()

    def _verify(self):
        pass

    @property
    def num_sides(self):
        """int: The number of sides in the current polygon."""
        return self._num_sides

    def __repr__(self):
        """Representation of current object.

        Returns:
            str: Object representation.
        """
        return '<{} (num_sides={:d})>'.format(
            self.__class__.__name__, self.num_sides)
