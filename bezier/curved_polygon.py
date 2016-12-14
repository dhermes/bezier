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

.. testsetup:: *

   import numpy as np
   import bezier
"""


from matplotlib import patches
from matplotlib import path as _path_mod
import matplotlib.pyplot as plt
import numpy as np
import six

from bezier import _helpers


class CurvedPolygon(object):
    """Represents an object defined by its curved boundary.

    The boundary is a piecewise defined collection of
    B |eacute| zier curves.

    .. note::

       The direction of the nodes in each :class:`.Curve`
       on the boundary is important: we check that one curve
       begins where the last one ended.

    .. image:: ../images/curved_polygon_constructor1.png
       :align: center

    .. doctest:: curved-polygon-constructor

       >>> import bezier
       >>> edge0 = bezier.Curve(np.array([
       ...     [0.0,  0.0],
       ...     [1.0, -1.0],
       ...     [2.0,  0.0],
       ... ]))
       >>> edge1 = bezier.Curve(np.array([
       ...     [2.0, 0.0],
       ...     [2.0, 1.0],
       ... ]))
       >>> edge2 = bezier.Curve(np.array([
       ...     [2.0, 1.0],
       ...     [1.0, 2.0],
       ...     [0.0, 1.0],
       ... ]))
       >>> edge3 = bezier.Curve(np.array([
       ...     [0.0, 1.0],
       ...     [0.0, 0.0],
       ... ]))
       >>> curved_poly = bezier.CurvedPolygon(
       ...     edge0, edge1, edge2, edge3)
       >>> curved_poly
       <CurvedPolygon (num_sides=4)>

    .. testcleanup:: curved-polygon-constructor

       import make_images
       make_images.curved_polygon_constructor1(curved_poly)

    Though the endpoints of edge pair of edges are verified to match,
    the curved polygon as a whole is not verified, so creating
    a curved polygon with self-intersections is possible:

    .. image:: ../images/curved_polygon_constructor2.png
       :align: center

    .. doctest:: curved-polygon-constructor-invalid

       >>> edge0 = bezier.Curve(np.array([
       ...     [0.0, 0.0],
       ...     [1.0, 0.0],
       ... ]))
       >>> edge1 = bezier.Curve(np.array([
       ...     [1.0 , 0.0],
       ...     [1.25, 0.5],
       ...     [1.0 , 1.0],
       ... ]))
       >>> edge2 = bezier.Curve(np.array([
       ...     [1.0, 1.0],
       ...     [2.0, 1.0],
       ... ]))
       >>> edge3 = bezier.Curve(np.array([
       ...     [2.0, 1.0 ],
       ...     [1.0, 0.75],
       ...     [0.0, 0.0 ],
       ... ]))
       >>> curved_poly = bezier.CurvedPolygon(
       ...     edge0, edge1, edge2, edge3)
       >>> curved_poly
       <CurvedPolygon (num_sides=4)>

    .. testcleanup:: curved-polygon-constructor-invalid

       import make_images
       make_images.curved_polygon_constructor2(curved_poly)

    Args:
        edges (Tuple[~bezier.curve.Curve, ...]): The boundary edges
            of the curved polygon.
    """

    def __init__(self, *edges):
        self._edges = edges
        self._num_sides = len(edges)
        self._verify()

    @staticmethod
    def _verify_pair(prev, curr):
        """Verify a pair of sides share an endpoint.

        .. note::

           This currently checks that edge endpoints match **exactly**
           but allowing some roundoff may be desired.

        Args:
            prev (.Curve): "Previous" curve at piecewise junction.
            curr (.Curve): "Next" curve at piecewise junction.

        Raises:
            ValueError: If the previous side is not in 2D.
            ValueError: If consecutive sides don't share an endpoint.
        """
        if prev.dimension != 2:
            raise ValueError('Curve not in R^2', prev)

        # pylint: disable=protected-access
        prev_nodes = prev._nodes
        curr_nodes = curr._nodes
        # pylint: enable=protected-access
        if not _helpers.vector_close(prev_nodes[-1, :], curr_nodes[0, :]):
            raise ValueError(
                'Not sufficiently close',
                'Consecutive sides do not have common endpoint',
                prev, curr)

    def _verify(self):
        """Verify that the edges define a curved polygon.

        This may not be entirely comprehensive, e.g. won't check
        self-intersection of the defined polygon.

        .. note::

           This currently checks that edge endpoints match **exactly**
           but allowing some roundoff may be desired.

        Raises:
            ValueError: If there are fewer than two sides.
            ValueError: If one of the sides is not in 2D.
            ValueError: If consecutive sides don't share an endpoint.
        """
        if self.num_sides < 2:
            raise ValueError('At least two sides required.')

        for prev, curr in six.moves.zip(self._edges, self._edges[1:]):
            self._verify_pair(prev, curr)

        # Now we check that the final edge wraps around.
        prev = self._edges[-1]
        curr = self._edges[0]
        self._verify_pair(prev, curr)

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

    def plot(self, pts_per_edge, color=None, ax=None, show=False):
        """Plot the current curved polygon.

        Args:
            pts_per_edge (int): Number of points to plot per curved edge.
            color (Optional[Tuple[float, float, float]]): Color as RGB profile.
            ax (Optional[matplotlib.artist.Artist]): matplotlib axis object
                to add plot to.
            show (Optional[bool]): Flag indicating if the plot should be
                shown.

        Returns:
            matplotlib.artist.Artist: The axis containing the plot. This
            may be a newly created axis.
        """
        if ax is None:
            fig = plt.figure()
            ax = fig.gca()

        s_vals = np.linspace(0.0, 1.0, pts_per_edge)

        # Evaluate points on each edge.
        all_points = []
        for edge in self._edges:
            points = edge.evaluate_multi(s_vals)
            line, = ax.plot(points[:, 0], points[:, 1], color=color)
            color = line.get_color()
            # Since the edges overlap, we leave out the first point in each.
            all_points.append(points[1:, :])

        polygon = np.vstack(all_points)
        path = _path_mod.Path(polygon)
        patch = patches.PathPatch(
            path, facecolor=color, alpha=0.6)
        ax.add_patch(patch)

        if show:
            plt.show()

        return ax
