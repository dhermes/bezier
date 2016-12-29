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

"""Helper for B |eacute| zier Surfaces / Triangles.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

.. testsetup:: *

   import numpy as np
   import bezier
"""


import itertools

from matplotlib import patches
from matplotlib import path as _path_mod
import matplotlib.pyplot as plt
import numpy as np
import six

from bezier import _base
from bezier import _intersection_helpers
from bezier import _surface_helpers
from bezier import curve as _curve_mod


_REPR_TEMPLATE = (
    '<{} (degree={:d}, dimension={:d}, base=({:g}, {:g}), width={:g})>')


class Surface(_base.Base):
    r"""Represents a B |eacute| zier `surface`_.

    .. _surface: https://en.wikipedia.org/wiki/B%C3%A9zier_triangle
    .. _unit simplex:
        https://en.wikipedia.org/wiki/Simplex#The_standard_simplex
    .. _barycentric coordinates:
        https://en.wikipedia.org/wiki/Barycentric_coordinate_system

    We define a B |eacute| zier triangle as a mapping from the
    `unit simplex`_ in 2D (i.e. the unit triangle) onto a surface in an
    arbitrary dimension. We use `barycentric coordinates`_

    .. math::

       \lambda_1 = 1 - s - t, \lambda_2 = s, \lambda_3 = t

    for points in

    .. math::

       \left\{(s, t) \mid 0 \leq s, t, s + t \leq 1\right\}.

    As with curves, using these weights we get convex combinations
    of points :math:`v_{i, j, k}` in some vector space:

    .. math::

       B\left(\lambda_1, \lambda_2, \lambda_3\right) =
           \sum_{i + j + k = d} \binom{d}{i \, j \, k}
           \lambda_1^i \lambda_2^j \lambda_3^k \cdot v_{i, j, k}

    .. note::

       We assume the nodes are ordered from left-to-right and from
       bottom-to-top. So for example, the linear triangle:

       .. code-block:: rest

          (0,0,1)

          (1,0,0)  (0,1,0)

       is ordered as

       .. math::

          \left[\begin{array}{c c c}
              v_{1,0,0} & v_{0,1,0} & v_{0,0,1} \end{array}\right]^T

       the quadratic triangle:

       .. code-block:: rest

          (0,0,2)

          (1,0,1)  (0,1,1)

          (2,0,0)  (1,1,0)  (0,2,0)

       is ordered as

       .. math::

          \left[\begin{array}{c c c c c c}
              v_{2,0,0} & v_{1,1,0} &
              v_{0,2,0} & v_{1,0,1} &
              v_{0,1,1} & v_{0,0,2} \end{array}\right]^T

       the cubic triangle:

       .. code-block:: rest

          (0,0,3)

          (1,0,2)  (0,1,2)

          (2,0,1)  (1,1,1)  (0,2,1)

          (3,0,0)  (2,1,0)  (1,2,0)  (0,3,0)

       is ordered as

       .. math::

          \left[\begin{array}{c c c c c c c c c c}
              v_{3,0,0} & v_{2,1,0} &
              v_{1,2,0} & v_{0,3,0} &
              v_{2,0,1} & v_{1,1,1} &
              v_{0,2,1} & v_{1,0,2} &
              v_{0,1,2} & v_{0,0,3} \end{array}\right]^T

       and so on.

    .. image:: ../images/surface_constructor.png
       :align: center

    .. doctest:: surface-constructor

       >>> import bezier
       >>> nodes = np.array([
       ...     [0.0  , 0.0  ],
       ...     [0.5  , 0.0  ],
       ...     [1.0  , 0.25 ],
       ...     [0.125, 0.5  ],
       ...     [0.375, 0.375],
       ...     [0.25 , 1.0  ],
       ... ])
       >>> surface = bezier.Surface(nodes)
       >>> surface
       <Surface (degree=2, dimension=2)>

    .. testcleanup:: surface-constructor

       import make_images
       make_images.surface_constructor(surface)

    Args:
        nodes (numpy.ndarray): The nodes in the surface. The rows
            represent each node while the columns are the dimension
            of the ambient space.
        base_x (Optional[float]): The :math:`x`-coordinate of the base
           vertex of the sub-triangle that this surface represents.
        base_y (Optional[float]): The :math:`y`-coordinate of the base
           vertex of the sub-triangle that this surface represents.
        width (Optional[float]): The width of the sub-triangle that
           this surface represents.
        _copy (bool): Flag indicating if the nodes should be copied before
            being stored. Defaults to :data:`True` since callers may
            freely mutate ``nodes`` after passing in.
    """

    _area = None
    _edges = None
    _is_valid = None

    def __init__(self, nodes, base_x=0.0, base_y=0.0, width=1.0, _copy=True):
        super(Surface, self).__init__(nodes, _copy=_copy)
        self._base_x = base_x
        self._base_y = base_y
        self._width = width

    def __repr__(self):
        """Representation of current object.

        Returns:
            str: Object representation.
        """
        if self._base_x == 0.0 and self._base_y == 0.0 and self._width == 1.0:
            return super(Surface, self).__repr__()
        else:
            return _REPR_TEMPLATE.format(
                self.__class__.__name__, self.degree, self.dimension,
                self._base_x, self._base_y, self._width)

    @staticmethod
    def _get_degree(num_nodes):
        """Get the degree of the current surface.

        Args:
            num_nodes (int): The number of control points for a
                B |eacute| zier surface.

        Returns:
            int: The degree :math:`d` such that :math:`(d + 1)(d + 2)/2`
            equals ``num_nodes``.

        Raises:
            ValueError: If ``num_nodes`` isn't a triangular number.
        """
        # 8 * num_nodes = 4(d + 1)(d + 2)
        #               = 4d^2 + 12d + 8
        #               = (2d + 3)^2 - 1
        d_float = 0.5 * (np.sqrt(8.0 * num_nodes + 1.0) - 3.0)
        d_int = int(np.round(d_float))
        if (d_int + 1) * (d_int + 2) == 2 * num_nodes:
            return d_int
        else:
            raise ValueError(num_nodes, 'not a triangular number')

    @property
    def area(self):  # pylint: disable=missing-returns-doc
        """float: The area of the current surface.

        Raises:
            NotImplementedError: If the area isn't already cached.
        """
        if self._area is None:
            raise NotImplementedError(
                'Area computation not yet implemented.')
        return self._area

    @property
    def width(self):
        """float: The "width" of the parameterized triangle.

        When re-parameterizing (e.g. via :meth:`subdivide`) we
        specialize the surface from the unit triangle to some
        sub-triangle. After doing this, we re-parameterize so that
        that sub-triangle is treated like the unit triangle.

        To track which sub-triangle we are in during the subdivision
        process, we use the coordinates of the base vertex as well
        as the "width" of each leg.

        .. image:: ../images/surface_width1.png
           :align: center

        .. testsetup:: surface-width1, surface-width2

           import numpy as np
           import bezier

           surface = bezier.Surface(np.array([
               [0.0, 0.0],
               [1.0, 0.0],
               [0.0, 1.0],
           ]))

        .. doctest:: surface-width1

           >>> surface.base_x, surface.base_y
           (0.0, 0.0)
           >>> surface.width
           1.0

        .. testcleanup:: surface-width1

           import make_images
           make_images.surface_width1(surface)

        Upon subdivision, the width halves (and potentially changes sign) and
        the vertex moves to one of four points:

        .. image:: ../images/surface_width2.png
           :align: center

        .. doctest:: surface-width2

           >>> _, sub_surface_b, sub_surface_c, _ = surface.subdivide()
           >>> sub_surface_b.base_x, sub_surface_b.base_y
           (0.5, 0.5)
           >>> sub_surface_b.width
           -0.5
           >>> sub_surface_c.base_x, sub_surface_c.base_y
           (0.5, 0.0)
           >>> sub_surface_c.width
           0.5

        .. testcleanup:: surface-width2

           import make_images
           make_images.surface_width2(sub_surface_b, sub_surface_c)
        """
        return self._width

    @property
    def base_x(self):
        """float: The ``x``-coordinate of the base vertex.

        See :meth:`width` for more detail.
        """
        return self._base_x

    @property
    def base_y(self):
        """float: The ``y``-coordinate of the base vertex.

        See :meth:`width` for more detail.
        """
        return self._base_y

    def _compute_edges(self):
        """Compute the edges of the current surface.

        Returns:
            Tuple[~curve.Curve, ~curve.Curve, ~curve.Curve]: The edges of
            the surface.
        """
        indices1 = slice(0, self.degree + 1)
        indices2 = np.empty(self.degree + 1, dtype=int)
        indices3 = np.empty(self.degree + 1, dtype=int)

        curr2 = self.degree
        curr3 = -1
        for i in six.moves.xrange(self.degree + 1):
            indices2[i] = curr2
            indices3[i] = curr3
            curr2 += self.degree - i
            curr3 -= i + 2

        edge1 = _curve_mod.Curve(self._nodes[indices1, :], _copy=False)
        edge2 = _curve_mod.Curve(self._nodes[indices2, :], _copy=False)
        edge3 = _curve_mod.Curve(self._nodes[indices3, :], _copy=False)
        _surface_helpers.edge_cycle(edge1, edge2, edge3)
        return edge1, edge2, edge3

    def _get_edges(self):
        """Get the edges for the current surface.

        If they haven't been computed yet, first compute and store them.

        This is provided as a means for internal calls to get the edges
        without copying (since :attr:`.edges` copies before giving to
        a user to keep the stored data immutable).

        Returns:
            Tuple[~bezier.curve.Curve, ~bezier.curve.Curve, \
                  ~bezier.curve.Curve]: The edges of
            the surface.
        """
        if self._edges is None:
            self._edges = self._compute_edges()
        return self._edges

    @property
    def edges(self):
        """tuple: The edges of the surface.

        .. doctest:: surface-edges
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.array([
           ...     [0.0   ,  0.0   ],
           ...     [0.5   , -0.1875],
           ...     [1.0   ,  0.0   ],
           ...     [0.1875,  0.5   ],
           ...     [0.625 ,  0.625 ],
           ...     [0.0   ,  1.0   ],
           ... ])
           >>> surface = bezier.Surface(nodes)
           >>> edge1, _, _ = surface.edges
           >>> edge1
           <Curve (degree=2, dimension=2)>
           >>> edge1.nodes
           array([[ 0.  ,  0.    ],
                  [ 0.5 , -0.1875],
                  [ 1.  ,  0.    ]])

        Returns:
            Tuple[~bezier.curve.Curve, ~bezier.curve.Curve, \
                  ~bezier.curve.Curve]: The edges of
            the surface.
        """
        edge1, edge2, edge3 = self._get_edges()
        # NOTE: It is crucial that we return copies here. Since the edges
        #       are cached, if they were mutable, callers could
        #       inadvertently mutate the cached value.
        edge1 = edge1.copy()
        edge2 = edge2.copy()
        edge3 = edge3.copy()
        _surface_helpers.edge_cycle(edge1, edge2, edge3)
        return edge1, edge2, edge3

    def evaluate_barycentric(self, lambda1, lambda2, lambda3):
        r"""Compute a point on the surface.

        Evaluates :math:`B\left(\lambda_1, \lambda_2, \lambda_3\right)`.

        .. image:: ../images/surface_evaluate_barycentric.png
           :align: center

        .. testsetup:: surface-barycentric, surface-barycentric-fail

           import numpy as np
           import bezier
           nodes = np.array([
               [0.0  , 0.0  ],
               [0.5  , 0.0  ],
               [1.0  , 0.25 ],
               [0.125, 0.5  ],
               [0.375, 0.375],
               [0.25 , 1.0  ],
           ])
           surface = bezier.Surface(nodes)

        .. doctest:: surface-barycentric
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.array([
           ...     [0.0  , 0.0  ],
           ...     [0.5  , 0.0  ],
           ...     [1.0  , 0.25 ],
           ...     [0.125, 0.5  ],
           ...     [0.375, 0.375],
           ...     [0.25 , 1.0  ],
           ... ])
           >>> surface = bezier.Surface(nodes)
           >>> point = surface.evaluate_barycentric(0.125, 0.125, 0.75)
           >>> point
           array([ 0.265625 , 0.73046875])

        .. testcleanup:: surface-barycentric

           import make_images
           make_images.surface_evaluate_barycentric(surface, point)

        However, this can't be used for points **outside** the
        reference triangle:

        .. doctest:: surface-barycentric-fail

           >>> surface.evaluate_barycentric(-0.25, 0.75, 0.5)
           Traceback (most recent call last):
             ...
           ValueError: ('Parameters must be positive', -0.25, 0.75, 0.5)

        or for non-Barycentric coordinates;

        .. doctest:: surface-barycentric-fail

           >>> surface.evaluate_barycentric(0.25, 0.25, 0.25)
           Traceback (most recent call last):
             ...
           ValueError: ('Values do not sum to 1', 0.25, 0.25, 0.25)

        Args:
            lambda1 (float): Parameter along the reference triangle.
            lambda2 (float): Parameter along the reference triangle.
            lambda3 (float): Parameter along the reference triangle.

        Returns:
            numpy.ndarray: The point on the surface (as a one dimensional
            NumPy array).

        Raises:
            ValueError: If the weights are not valid barycentric
                coordinates, i.e. they don't sum to ``1``.
            ValueError: If some weights are negative.
        """
        weights_total = lambda1 + lambda2 + lambda3
        if not np.allclose(weights_total, 1.0):
            raise ValueError('Weights do not sum to 1',
                             lambda1, lambda2, lambda3)
        if lambda1 < 0.0 or lambda2 < 0.0 or lambda3 < 0.0:
            raise ValueError('Weights must be positive',
                             lambda1, lambda2, lambda3)

        if self.degree == 1:
            weights = np.array([
                [lambda1, lambda2, lambda3],
            ])
        elif self.degree == 2:
            weights = np.array([
                [
                    lambda1 * lambda1,
                    2.0 * lambda1 * lambda2,
                    lambda2 * lambda2,
                    2.0 * lambda1 * lambda3,
                    2.0 * lambda2 * lambda3,
                    lambda3 * lambda3,
                ]
            ])
        elif self.degree == 3:
            weights = np.array([
                [
                    lambda1 * lambda1 * lambda1,
                    3.0 * lambda1 * lambda1 * lambda2,
                    3.0 * lambda1 * lambda2 * lambda2,
                    lambda2 * lambda2 * lambda2,
                    3.0 * lambda1 * lambda1 * lambda3,
                    6.0 * lambda1 * lambda2 * lambda3,
                    3.0 * lambda2 * lambda2 * lambda3,
                    3.0 * lambda1 * lambda3 * lambda3,
                    3.0 * lambda2 * lambda3 * lambda3,
                    lambda3 * lambda3 * lambda3,
                ]
            ])
        else:
            result = self._nodes
            for reduced_deg in six.moves.xrange(self.degree, 0, -1):
                result = _surface_helpers.de_casteljau_one_round(
                    result, reduced_deg, lambda1, lambda2, lambda3)
            return result.flatten()

        return weights.dot(self._nodes).flatten()  # pylint: disable=no-member

    def evaluate_cartesian(self, s, t):
        r"""Compute a point on the surface.

        Evaluates :math:`B\left(1 - s - t, s, t\right)` by calling
        :meth:`evaluate_barycentric`.

        Args:
            s (float): Parameter along the reference triangle.
            t (float): Parameter along the reference triangle.

        Returns:
            numpy.ndarray: The point on the surface (as a one dimensional
            NumPy array).
        """
        return self.evaluate_barycentric(1.0 - s - t, s, t)

    def evaluate_multi(self, param_vals):
        r"""Compute multiple points on the surface.

        If ``param_vals`` has two columns, this method treats
        them as Cartesian:

        .. image:: ../images/surface_evaluate_multi1.png
           :align: center

        .. doctest:: surface-eval-multi1
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.array([
           ...     [ 0.0, 0.0],
           ...     [ 2.0, 1.0],
           ...     [-3.0, 2.0],
           ... ])
           >>> surface = bezier.Surface(nodes)
           >>> surface
           <Surface (degree=1, dimension=2)>
           >>> param_vals = np.array([
           ...     [0.0  , 0.0  ],
           ...     [0.125, 0.625],
           ...     [0.5  , 0.5  ],
           ... ])
           >>> points = surface.evaluate_multi(param_vals)
           >>> points
           array([[ 0.   , 0.   ],
                  [-1.625, 1.375],
                  [-0.5  , 1.5  ]])

        .. testcleanup:: surface-eval-multi1

           import make_images
           make_images.surface_evaluate_multi1(surface, points)

        and if ``param_vals`` has three columns, treats them as Barycentric:

        .. image:: ../images/surface_evaluate_multi2.png
           :align: center

        .. doctest:: surface-eval-multi2
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.array([
           ...     [ 0. , 0.  ],
           ...     [ 1. , 0.75],
           ...     [ 2. , 1.  ],
           ...     [-1.5, 1.  ],
           ...     [-0.5, 1.5 ],
           ...     [-3. , 2.  ],
           ... ])
           >>> surface = bezier.Surface(nodes)
           >>> surface
           <Surface (degree=2, dimension=2)>
           >>> param_vals = np.array([
           ...     [0.   , 0.25, 0.75 ],
           ...     [1.   , 0.  , 0.   ],
           ...     [0.25 , 0.5 , 0.25 ],
           ...     [0.375, 0.25, 0.375],
           ... ])
           >>> points = surface.evaluate_multi(param_vals)
           >>> points
           array([[-1.75  , 1.75    ],
                  [ 0.    , 0.      ],
                  [ 0.25  , 1.0625  ],
                  [-0.625 , 1.046875]])

        .. testcleanup:: surface-eval-multi2

           import make_images
           make_images.surface_evaluate_multi2(surface, points)

        .. note::

           This currently just uses :meth:`evaluate_cartesian` and
           :meth:`evaluate_barycentric` so is less
           performant than it could be.

        Args:
            param_vals (numpy.ndarray): Array of parameter values (as a
                2D array).

        Returns:
            numpy.ndarray: The point on the surface.

        Raises:
            ValueError: If ``param_vals`` is not a 2D array.
            ValueError: If ``param_vals`` doesn't have 2 or 3 columns.
        """
        if param_vals.ndim != 2:
            raise ValueError('Parameter values must be 2D array')
        num_vals, num_cols = param_vals.shape

        if num_cols == 2:
            transform = self.evaluate_cartesian
        elif num_cols == 3:
            transform = self.evaluate_barycentric
        else:
            raise ValueError(
                'Parameter values must either be Barycentric or Cartesian')

        result = np.empty((num_vals, self.dimension))
        for index in six.moves.xrange(num_vals):
            result[index, :] = transform(*param_vals[index, :])
        return result

    @staticmethod
    def _add_patch(ax, color, edge1, edge2, edge3):
        """Add a polygonal surface patch to a plot.

        Args:
            ax (matplotlib.artist.Artist): A matplotlib axis.
            color (Tuple[float, float, float]): Color as RGB profile.
            edge1 (numpy.ndarray): 2D array of points along a curved edge.
            edge2 (numpy.ndarray): 2D array of points along a curved edge.
            edge3 (numpy.ndarray): 2D array of points along a curved edge.
        """
        # Since the edges overlap, we leave out the first point in each.
        polygon = np.vstack([
            edge1[1:, :],
            edge2[1:, :],
            edge3[1:, :],
        ])
        path = _path_mod.Path(polygon)
        patch = patches.PathPatch(
            path, facecolor=color, alpha=0.6)
        ax.add_patch(patch)

    def plot(self, pts_per_edge, color=None, ax=None,
             with_nodes=False, show=False):
        """Plot the current surface.

        Args:
            pts_per_edge (int): Number of points to plot per edge.
            color (Optional[Tuple[float, float, float]]): Color as RGB profile.
            ax (Optional[matplotlib.artist.Artist]): matplotlib axis object
                to add plot to.
            with_nodes (Optional[bool]): Determines if the control points
                should be added to the plot. Off by default.
            show (Optional[bool]): Flag indicating if the plot should be
                shown.

        Returns:
            matplotlib.artist.Artist: The axis containing the plot. This
            may be a newly created axis.

        Raises:
            NotImplementedError: If the surface's dimension is not ``2``.
        """
        if self.dimension != 2:
            raise NotImplementedError('2D is the only supported dimension',
                                      'Current dimension', self.dimension)

        edge1, edge2, edge3 = self.edges
        s_vals = np.linspace(0.0, 1.0, pts_per_edge)

        points1 = edge1.evaluate_multi(s_vals)
        points2 = edge2.evaluate_multi(s_vals)
        points3 = edge3.evaluate_multi(s_vals)

        if ax is None:
            fig = plt.figure()
            ax = fig.gca()

        line, = ax.plot(points1[:, 0], points1[:, 1], color=color)
        color = line.get_color()
        ax.plot(points2[:, 0], points2[:, 1], color=color)
        ax.plot(points3[:, 0], points3[:, 1], color=color)

        self._add_patch(ax, color, points1, points2, points3)

        if with_nodes:
            ax.plot(self._nodes[:, 0], self._nodes[:, 1],
                    color='black', marker='o', linestyle='None')

        if show:
            plt.show()

        return ax

    def subdivide(self):
        r"""Split the surface into four sub-surfaces.

        Does so by taking the unit triangle (i.e. the domain
        of the surface) and splitting it into four sub-triangles

        .. image:: ../images/surface_subdivide1.png
           :align: center

        Then the surface is re-parameterized via the map to/from the
        given sub-triangles and the unit triangle.

        For example, when a degree two surface is subdivided:

        .. image:: ../images/surface_subdivide2.png
           :align: center

        .. doctest:: surface-subdivide
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.array([
           ...     [-1.0 , 0.0 ],
           ...     [ 0.5 , 0.5 ],
           ...     [ 2.0 , 0.0 ],
           ...     [ 0.25, 1.75],
           ...     [ 2.0 , 3.0 ],
           ...     [ 0.0 , 4.0 ],
           ... ])
           >>> surface = bezier.Surface(nodes)
           >>> _, sub_surface_b, _, _ = surface.subdivide()
           >>> sub_surface_b
           <Surface (degree=2, dimension=2, base=(0.5, 0.5), width=-0.5)>
           >>> sub_surface_b.nodes
           array([[ 1.5   , 2.5   ],
                  [ 0.6875, 2.3125],
                  [-0.125 , 1.875 ],
                  [ 1.1875, 1.3125],
                  [ 0.4375, 1.3125],
                  [ 0.5   , 0.25  ]])

        .. testcleanup:: surface-subdivide

           import make_images
           make_images.surface_subdivide1()
           make_images.surface_subdivide2(surface, sub_surface_b)

        Returns:
            Tuple[Surface, Surface, Surface, Surface]: The lower left, central,
            lower right and upper left sub-surfaces (in that order).
        """
        if self.degree == 1:
            # pylint: disable=no-member
            new_nodes = _surface_helpers.LINEAR_SUBDIVIDE.dot(self._nodes)
            # pylint: enable=no-member
            nodes_a = new_nodes[(0, 1, 3), :]
            nodes_b = new_nodes[(4, 3, 1), :]
            nodes_c = new_nodes[(1, 2, 4), :]
            nodes_d = new_nodes[3:, :]
        elif self.degree == 2:
            # pylint: disable=no-member
            new_nodes = _surface_helpers.QUADRATIC_SUBDIVIDE.dot(self._nodes)
            # pylint: enable=no-member
            nodes_a = new_nodes[(0, 1, 2, 5, 6, 9), :]
            nodes_b = new_nodes[(11, 10, 9, 7, 6, 2), :]
            nodes_c = new_nodes[(2, 3, 4, 7, 8, 11), :]
            nodes_d = new_nodes[9:, :]
        elif self.degree == 3:
            # pylint: disable=no-member
            new_nodes = _surface_helpers.CUBIC_SUBDIVIDE.dot(self._nodes)
            # pylint: enable=no-member
            nodes_a = new_nodes[(0, 1, 2, 3, 7, 8, 9, 13, 14, 18), :]
            nodes_b = new_nodes[(21, 20, 19, 18, 16, 15, 14, 10, 9, 3), :]
            nodes_c = new_nodes[(3, 4, 5, 6, 10, 11, 12, 16, 17, 21), :]
            nodes_d = new_nodes[18:, :]
        elif self.degree == 4:
            # pylint: disable=no-member
            new_nodes = _surface_helpers.QUARTIC_SUBDIVIDE.dot(self._nodes)
            # pylint: enable=no-member
            nodes_a = new_nodes[
                (0, 1, 2, 3, 4, 9, 10, 11, 12, 17, 18, 19, 24, 25, 30), :]
            nodes_b = new_nodes[
                (34, 33, 32, 31, 30, 28, 27, 26, 25, 21, 20, 19, 13, 12, 4), :]
            nodes_c = new_nodes[
                (4, 5, 6, 7, 8, 13, 14, 15, 16, 21, 22, 23, 28, 29, 34), :]
            nodes_d = new_nodes[30:, :]
        else:
            nodes_a = _surface_helpers.specialize_surface(
                self._nodes, self.degree,
                (1.0, 0.0, 0.0), (0.5, 0.5, 0.0), (0.5, 0.0, 0.5))
            nodes_b = _surface_helpers.specialize_surface(
                self._nodes, self.degree,
                (0.0, 0.5, 0.5), (0.5, 0.0, 0.5), (0.5, 0.5, 0.0))
            nodes_c = _surface_helpers.specialize_surface(
                self._nodes, self.degree,
                (0.5, 0.5, 0.0), (0.0, 1.0, 0.0), (0.0, 0.5, 0.5))
            nodes_d = _surface_helpers.specialize_surface(
                self._nodes, self.degree,
                (0.5, 0.0, 0.5), (0.0, 0.5, 0.5), (0.0, 0.0, 1.0))

        half_width = 0.5 * self._width
        shifted_x = self._base_x + half_width
        shifted_y = self._base_y + half_width
        return (
            Surface(nodes_a, base_x=self._base_x, base_y=self._base_y,
                    width=half_width, _copy=False),
            Surface(nodes_b, base_x=shifted_x, base_y=shifted_y,
                    width=-half_width, _copy=False),
            Surface(nodes_c, base_x=shifted_x, base_y=self._base_y,
                    width=half_width, _copy=False),
            Surface(nodes_d, base_x=self._base_x, base_y=shifted_y,
                    width=half_width, _copy=False),
        )

    def _compute_valid(self):
        r"""Determines if the current surface is "valid".

        Does this by checking if the Jacobian of the map from the
        reference triangle is nonzero.

        Returns:
            bool: Flag indicating if the current surface is valid.

        Raises:
            NotImplementedError: If the degree is greater than 3.
            NotImplementedError: If the surface is quadratic or cubic
                but in dimension other than :math:`\mathbf{R}^2`.
        """
        if self.degree == 1:
            # In the linear case, we are only invalid if the points
            # are collinear.
            # pylint: disable=no-member
            first_deriv = self._nodes[1:, :] - self._nodes[:-1, :]
            # pylint: enable=no-member
            return np.linalg.matrix_rank(first_deriv) == 2
        elif self.degree in (2, 3):
            if self.dimension != 2:
                raise NotImplementedError(
                    'Cubic/quadratic validity check only implemented in R^2')

            if self.degree == 2:
                bernstein = _surface_helpers.quadratic_jacobian_polynomial(
                    self._nodes)
            else:
                bernstein = _surface_helpers.cubic_jacobian_polynomial(
                    self._nodes)

            # Form the polynomial p(s, t) as a Surface.
            jac_poly = Surface(bernstein, _copy=False)
            # Find the sign of the polynomial, where 0 means mixed.
            poly_sign = _surface_helpers.polynomial_sign(jac_poly)
            return poly_sign != 0
        else:
            raise NotImplementedError(
                'Degrees 1, 2 and 3 only supported at this time')

    @property
    def is_valid(self):
        """bool: Flag indicating if the surface is "valid".

        Here, "valid" means there are no self-intersections or
        singularities.

        This checks if the Jacobian of the map from the reference
        triangle is nonzero. For example, a linear "surface"
        with collinear points is invalid:

        .. image:: ../images/surface_is_valid1.png
           :align: center

        .. doctest:: surface-is-valid1

           >>> nodes = np.array([
           ...     [0.0, 0.0],
           ...     [1.0, 1.0],
           ...     [2.0, 2.0],
           ... ])
           >>> surface = bezier.Surface(nodes)
           >>> surface.is_valid
           False

        .. testcleanup:: surface-is-valid1

           import make_images
           make_images.surface_is_valid1(surface)

        while a quadratic surface with one straight side:

        .. image:: ../images/surface_is_valid2.png
           :align: center

        .. doctest:: surface-is-valid2

           >>> nodes = np.array([
           ...     [ 0.0  , 0.0  ],
           ...     [ 0.5  , 0.125],
           ...     [ 1.0  , 0.0  ],
           ...     [-0.125, 0.5  ],
           ...     [ 0.5  , 0.5  ],
           ...     [ 0.0  , 1.0  ],
           ... ])
           >>> surface = bezier.Surface(nodes)
           >>> surface.is_valid
           True

        .. testcleanup:: surface-is-valid2

           import make_images
           make_images.surface_is_valid2(surface)

        though not all higher degree surfaces are valid:

        .. image:: ../images/surface_is_valid3.png
           :align: center

        .. doctest:: surface-is-valid3

           >>> nodes = np.array([
           ...     [1.0, 0.0],
           ...     [0.0, 0.0],
           ...     [1.0, 1.0],
           ...     [0.0, 0.0],
           ...     [0.0, 0.0],
           ...     [0.0, 1.0],
           ... ])
           >>> surface = bezier.Surface(nodes)
           >>> surface.is_valid
           False

        .. testcleanup:: surface-is-valid3

           import make_images
           make_images.surface_is_valid3(surface)
        """
        if self._is_valid is None:
            self._is_valid = self._compute_valid()
        return self._is_valid

    def locate(self, point):
        r"""Find a point on the current surface.

        Solves for :math:`s` and :math:`t` in :math:`B(s, t) = p`.

        .. note::

           A unique solution is only guaranteed if the current surface is
           valid. This code assumes a valid surface, but doesn't check.

        .. image:: ../images/surface_locate.png
           :align: center

        .. doctest:: surface-locate

           >>> surface = bezier.Surface(np.array([
           ...     [0.0 ,  0.0 ],
           ...     [0.5 , -0.25],
           ...     [1.0 ,  0.0 ],
           ...     [0.25,  0.5 ],
           ...     [0.75,  0.75],
           ...     [0.0 ,  1.0 ],
           ... ]))
           >>> point = np.array([[0.59375, 0.25]])
           >>> s, t = surface.locate(point)
           >>> s
           0.5
           >>> t
           0.25

        .. testcleanup:: surface-locate

           import make_images
           make_images.surface_locate(surface, point)

        Args:
            point (numpy.ndarray): A (``1xD``) point on the surface,
                where :math:`D` is the dimension of the surface.

        Returns:
            Optional[Tuple[float, float]]: The :math:`s` and :math:`t`
            values corresponding to ``x_val`` and ``y_val`` or
            :data:`None` if the point is not on the ``surface``.

        Raises:
            NotImplementedError: If the surface isn't in :math:`\mathbf{R}^2`.
            ValueError: If the dimension of the ``point`` doesn't match the
                dimension of the current surface.
        """
        if self.dimension != 2:
            raise NotImplementedError('Only 2D surfaces supported.')

        if point.shape != (1, self.dimension):
            raise ValueError('Point is not in same dimension as surface',
                             point, 'Shape expected:', (1, self.dimension))

        return _surface_helpers.locate_point(self, point[0, 0], point[0, 1])

    def intersect(self, other):
        """Find the common intersection with another surface.

        Args:
            other (Surface): Other surface to intersect with.

        Returns:
            list: List of intersection objects (possibly empty).

        Raises:
            TypeError: If ``other`` is not a surface.
            NotImplementedError: If both surfaces aren't two-dimensional.
        """
        if not isinstance(other, Surface):
            raise TypeError('Can only intersect with another surface',
                            'Received', other)
        if self.dimension != 2 or other.dimension != 2:
            raise NotImplementedError(
                'Intersection only implemented in 2D')

        edges1 = self._get_edges()
        lin1 = six.moves.map(
            _intersection_helpers.Linearization.from_shape,
            edges1)
        edges2 = other._get_edges()  # pylint: disable=protected-access
        lin2 = six.moves.map(
            _intersection_helpers.Linearization.from_shape,
            edges2)

        edge_pairs = itertools.product(lin1, lin2)
        intersections = _intersection_helpers.all_intersections(
            edge_pairs)

        # Classify each intersection.
        duplicates = []
        uniques = []
        for intersection in intersections:
            changed = _surface_helpers.handle_corners(intersection)
            if changed:
                duplicates.append(intersection)
            else:
                interior = _surface_helpers.classify_intersection(
                    intersection)
                intersection.interior_curve = interior
                uniques.append(intersection)

        # NOTE: This is an extra check to make sure we found duplicate
        #       intersections at corners. It is not strictly needed,
        #       it is just provided as a sanity check.
        _surface_helpers.verify_duplicates(duplicates, uniques)
        return _surface_helpers.combine_intersections(
            uniques, self, other)

    def elevate(self):
        r"""Return a degree-elevated version of the current surface.

        Does this by converting the current nodes
        :math:`\left\{v_{i, j, k}\right\}_{i + j + k = d}` to new nodes
        :math:`\left\{w_{i, j, k}\right\}_{i + j + k = d + 1}`. Does so
        by re-writing

        .. math::

           E\left(\lambda_1, \lambda_2, \lambda_3\right) =
               \left(\lambda_1 + \lambda_2 + \lambda_3\right)
               B\left(\lambda_1, \lambda_2, \lambda_3\right) =
               \sum_{i + j + k = d + 1} \binom{d + 1}{i \, j \, k}
               \lambda_1^i \lambda_2^j \lambda_3^k \cdot w_{i, j, k}

        In this form, we must have

        .. math::

           \begin{align*}
           \binom{d + 1}{i \, j \, k} \cdot w_{i, j, k} &=
               \binom{d}{i - 1 \, j \, k} \cdot v_{i - 1, j, k} +
               \binom{d}{i \, j - 1 \, k} \cdot v_{i, j - 1, k} +
               \binom{d}{i \, j \, k - 1} \cdot v_{i, j, k - 1} \\
           \Longleftrightarrow (d + 1) \cdot w_{i, j, k} &=
               i \cdot v_{i - 1, j, k} + j \cdot v_{i, j - 1, k} +
               k \cdot v_{i, j, k - 1}
           \end{align*}

        where we assume that, for example, :math:`v_{i, j, k - 1}` is
        :math:`0` (or any other unused value) if :math:`k = 0`.

        .. doctest:: surface-elevate
           :options: +NORMALIZE_WHITESPACE

           >>> surface = bezier.Surface(np.array([
           ...     [0.0, 0.0],
           ...     [1.0, 0.0],
           ...     [0.0, 1.0],
           ... ]))
           >>> surface
           <Surface (degree=1, dimension=2)>
           >>> new_surface = surface.elevate()
           >>> new_surface
           <Surface (degree=2, dimension=2)>
           >>> new_surface.nodes
           array([[ 0. , 0. ],
                  [ 0.5, 0. ],
                  [ 1. , 0. ],
                  [ 0. , 0.5],
                  [ 0.5, 0.5],
                  [ 0. , 1. ]])

        Returns:
            Surface: The degree-elevated surface.
        """
        num_nodes, _ = self._nodes.shape
        # (d + 1)(d + 2)/2 --> (d + 2)(d + 3)/2
        num_new = num_nodes + self.degree + 2
        new_nodes = np.zeros((num_new, self.dimension))

        # NOTE: We start from the index triples (i, j, k) for the current
        #       nodes and map them onto (i + 1, j, k), etc. This index
        #       tracking is also done in :func:`.de_casteljau_one_round`.
        index = 0
        # parent_i1 = index + k
        # parent_i2 = index + k + 1
        # parent_i3 = index + degree + 2
        parent_i1 = 0
        parent_i2 = 1
        parent_i3 = self.degree + 2
        for k in six.moves.xrange(self.degree + 1):
            for j in six.moves.xrange(self.degree + 1 - k):
                i = self.degree - j - k
                new_nodes[parent_i1, :] += (i + 1) * self._nodes[index, :]
                new_nodes[parent_i2, :] += (j + 1) * self._nodes[index, :]
                new_nodes[parent_i3, :] += (k + 1) * self._nodes[index, :]
                # Update all the indices.
                parent_i1 += 1
                parent_i2 += 1
                parent_i3 += 1
                index += 1

            # Update the indices that depend on k.
            parent_i1 += 1
            parent_i2 += 1

        # Hold off on division until the end, to (attempt to) avoid round-off.
        denominator = self.degree + 1.0
        new_nodes /= denominator

        return Surface(new_nodes, _copy=False)
