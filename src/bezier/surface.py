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

"""Helper for B |eacute| zier Surfaces / Triangles.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

.. testsetup:: *

   import numpy as np
   import bezier
"""


import numpy as np
import six

from bezier import _algebraic_intersection
from bezier import _base
from bezier import _geometric_intersection
from bezier import _intersection_helpers
from bezier import _plot_helpers
from bezier import _surface_helpers
from bezier import _surface_intersection
from bezier import curve as _curve_mod


_REPR_TEMPLATE = (
    '<{} (degree={:d}, dimension={:d}, base=({:g}, {:g}), width={:g})>')
_LOCATE_ERROR_TEMPLATE = (
    'Dimension mismatch: This surface is {:d}-dimensional, so the point '
    'should be a 1x{:d} NumPy array. Instead the point {} has dimensions {}.')
_STRATEGY = _intersection_helpers.IntersectionStrategy
_INTERSECTION_T = _geometric_intersection.BoxIntersectionType.INTERSECTION


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

    for points in the unit triangle
    :math:`\left\{(s, t) \mid 0 \leq s, t, s + t \leq 1\right\}`:

    .. image:: ../images/unit_triangle.png
       :align: center

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

       The index formula

       .. math::

          j + \frac{k}{2} \left(2 (i + j) + k + 3\right)

       can be used to map a triple :math:`(i, j, k)` onto the
       corresponding linear index, but it is not particularly insightful
       or useful.

    .. image:: ../images/surface_constructor.png
       :align: center

    .. doctest:: surface-constructor

       >>> import bezier
       >>> nodes = np.asfortranarray([
       ...     [0.0  , 0.0  ],
       ...     [0.5  , 0.0  ],
       ...     [1.0  , 0.25 ],
       ...     [0.125, 0.5  ],
       ...     [0.375, 0.375],
       ...     [0.25 , 1.0  ],
       ... ])
       >>> surface = bezier.Surface(nodes, degree=2)
       >>> surface
       <Surface (degree=2, dimension=2)>

    .. testcleanup:: surface-constructor

       import make_images
       make_images.unit_triangle()
       make_images.surface_constructor(surface)

    Args:
        nodes (numpy.ndarray): The nodes in the surface. The rows
            represent each node while the columns are the dimension
            of the ambient space.
        degree (int): The degree of the surface. This is assumed to
            correctly correspond to the number of ``nodes``. Use
            :meth:`from_nodes` if the degree has not yet been computed.
        base_x (Optional[float]): The :math:`x`-coordinate of the base
           vertex of the sub-triangle that this surface represents.
           See :attr:`width` for more info.
        base_y (Optional[float]): The :math:`y`-coordinate of the base
           vertex of the sub-triangle that this surface represents.
           See :attr:`width` for more info.
        width (Optional[float]): The width of the sub-triangle that
           this surface represents. See :attr:`width` for more info.
        _copy (bool): Flag indicating if the nodes should be copied before
            being stored. Defaults to :data:`True` since callers may
            freely mutate ``nodes`` after passing in.
    """

    __slots__ = (
        '_dimension', '_nodes',  # From base class
        '_degree', '_base_x', '_base_y', '_width',  # From constructor
        '_area', '_edges', '_is_valid',  # Empty defaults
    )

    def __init__(self, nodes, degree, base_x=0.0, base_y=0.0,
                 width=1.0, _copy=True):
        super(Surface, self).__init__(nodes, _copy=_copy)
        self._degree = degree
        self._base_x = base_x
        self._base_y = base_y
        self._width = width
        self._area = None
        self._edges = None
        self._is_valid = None

    @classmethod
    def from_nodes(cls, nodes, base_x=0.0, base_y=0.0, width=1.0, _copy=True):
        """Create a :class:`.Surface` from nodes.

        Computes the ``degree`` based on the shape of ``nodes``.

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

        Returns:
            Surface: The constructed surface.
        """
        num_nodes, _ = nodes.shape
        degree = cls._get_degree(num_nodes)
        return cls(nodes, degree, base_x=base_x, base_y=base_y,
                   width=width, _copy=True)

    def __repr__(self):
        """Representation of current object.

        Returns:
            str: Object representation.
        """
        if self._base_x == 0.0 and self._base_y == 0.0 and self._width == 1.0:
            return super(Surface, self).__repr__()
        else:
            return _REPR_TEMPLATE.format(
                self.__class__.__name__, self._degree, self._dimension,
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

    # pylint: disable=missing-return-doc,missing-return-type-doc
    @property
    def area(self):
        """float: The area of the current surface.

        Raises:
            NotImplementedError: If the area isn't already cached.
        """
        if self._area is None:
            raise NotImplementedError(
                'Area computation not yet implemented.')
        return self._area
    # pylint: enable=missing-return-doc,missing-return-type-doc

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

           nodes = np.asfortranarray([
               [0.0, 0.0],
               [1.0, 0.0],
               [0.0, 1.0],
           ])
           surface = bezier.Surface(nodes, degree=1)

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

        See :attr:`width` for more detail.
        """
        return self._base_x

    @property
    def base_y(self):
        """float: The ``y``-coordinate of the base vertex.

        See :attr:`width` for more detail.
        """
        return self._base_y

    def _compute_edges(self):
        """Compute the edges of the current surface.

        Returns:
            Tuple[~curve.Curve, ~curve.Curve, ~curve.Curve]: The edges of
            the surface.
        """
        nodes1, nodes2, nodes3 = _surface_helpers.compute_edge_nodes(
            self._nodes, self._degree)
        edge1 = _curve_mod.Curve(nodes1, self._degree, _copy=False)
        edge2 = _curve_mod.Curve(nodes2, self._degree, _copy=False)
        edge3 = _curve_mod.Curve(nodes3, self._degree, _copy=False)
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
        """The edges of the surface.

        .. doctest:: surface-edges
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0   ,  0.0   ],
           ...     [0.5   , -0.1875],
           ...     [1.0   ,  0.0   ],
           ...     [0.1875,  0.5   ],
           ...     [0.625 ,  0.625 ],
           ...     [0.0   ,  1.0   ],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
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
        edge1 = edge1._copy()  # pylint: disable=protected-access
        edge2 = edge2._copy()  # pylint: disable=protected-access
        edge3 = edge3._copy()  # pylint: disable=protected-access
        return edge1, edge2, edge3

    @staticmethod
    def _verify_barycentric(lambda1, lambda2, lambda3):
        """Verifies that weights are barycentric and on the reference triangle.

        I.e., checks that they sum to one and are all non-negative.

        Args:
            lambda1 (float): Parameter along the reference triangle.
            lambda2 (float): Parameter along the reference triangle.
            lambda3 (float): Parameter along the reference triangle.

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

    def evaluate_barycentric(self, lambda1, lambda2, lambda3, _verify=True):
        r"""Compute a point on the surface.

        Evaluates :math:`B\left(\lambda_1, \lambda_2, \lambda_3\right)`.

        .. image:: ../images/surface_evaluate_barycentric.png
           :align: center

        .. testsetup:: surface-barycentric, surface-barycentric-fail,
                       surface-barycentric-no-verify

           import numpy as np
           import bezier
           nodes = np.asfortranarray([
               [0.0  , 0.0  ],
               [0.5  , 0.0  ],
               [1.0  , 0.25 ],
               [0.125, 0.5  ],
               [0.375, 0.375],
               [0.25 , 1.0  ],
           ])
           surface = bezier.Surface(nodes, degree=2)

        .. doctest:: surface-barycentric
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0  , 0.0  ],
           ...     [0.5  , 0.0  ],
           ...     [1.0  , 0.25 ],
           ...     [0.125, 0.5  ],
           ...     [0.375, 0.375],
           ...     [0.25 , 1.0  ],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
           >>> point = surface.evaluate_barycentric(0.125, 0.125, 0.75)
           >>> point
           array([[ 0.265625 , 0.73046875]])

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

        or for non-barycentric coordinates;

        .. doctest:: surface-barycentric-fail

           >>> surface.evaluate_barycentric(0.25, 0.25, 0.25)
           Traceback (most recent call last):
             ...
           ValueError: ('Values do not sum to 1', 0.25, 0.25, 0.25)

        However, these "invalid" inputs can be used if ``_verify`` is
        :data:`False`.

        .. doctest:: surface-barycentric-no-verify
           :options: +NORMALIZE_WHITESPACE

           >>> surface.evaluate_barycentric(-0.25, 0.75, 0.5, _verify=False)
           array([[ 0.6875 , 0.546875]])
           >>> surface.evaluate_barycentric(0.25, 0.25, 0.25, _verify=False)
           array([[ 0.203125, 0.1875 ]])

        Args:
            lambda1 (float): Parameter along the reference triangle.
            lambda2 (float): Parameter along the reference triangle.
            lambda3 (float): Parameter along the reference triangle.
            _verify (Optional[bool]): Indicates if the barycentric coordinates
                should be verified as summing to one and all non-negative (i.e.
                verified as barycentric). Can either be used to evaluate at
                points outside the domain, or to save time when the caller
                already knows the input is verified. Defaults to :data:`True`.

        Returns:
            numpy.ndarray: The point on the surface (as a two dimensional
            NumPy array with a single row).

        Raises:
            ValueError: If the weights are not valid barycentric
                coordinates, i.e. they don't sum to ``1``. (Won't raise if
                ``_verify=False``.)
            ValueError: If some weights are negative. (Won't raise if
                ``_verify=False``.)
        """
        if _verify:
            self._verify_barycentric(lambda1, lambda2, lambda3)

        return _surface_helpers.evaluate_barycentric(
            self._nodes, self._degree, lambda1, lambda2, lambda3)

    def evaluate_barycentric_multi(self, param_vals, _verify=True):
        r"""Compute multiple points on the surface.

        Assumes ``param_vals`` has three columns of barycentric coordinates.
        See :meth:`evaluate_barycentric` for more details on how each row of
        parameter values is evaluated.

        .. image:: ../images/surface_evaluate_barycentric_multi.png
           :align: center

        .. doctest:: surface-eval-multi2
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [ 0. , 0.  ],
           ...     [ 1. , 0.75],
           ...     [ 2. , 1.  ],
           ...     [-1.5, 1.  ],
           ...     [-0.5, 1.5 ],
           ...     [-3. , 2.  ],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
           >>> surface
           <Surface (degree=2, dimension=2)>
           >>> param_vals = np.asfortranarray([
           ...     [0.   , 0.25, 0.75 ],
           ...     [1.   , 0.  , 0.   ],
           ...     [0.25 , 0.5 , 0.25 ],
           ...     [0.375, 0.25, 0.375],
           ... ])
           >>> points = surface.evaluate_barycentric_multi(param_vals)
           >>> points
           array([[-1.75  , 1.75    ],
                  [ 0.    , 0.      ],
                  [ 0.25  , 1.0625  ],
                  [-0.625 , 1.046875]])

        .. testcleanup:: surface-eval-multi2

           import make_images
           make_images.surface_evaluate_barycentric_multi(surface, points)

        Args:
            param_vals (numpy.ndarray): Array of parameter values (as a
                ``Nx3`` array).
            _verify (Optional[bool]): Indicates if the coordinates should be
                verified. See :meth:`evaluate_barycentric`. Defaults to
                :data:`True`. Will also double check that ``param_vals``
                is the right shape.

        Returns:
            numpy.ndarray: The points on the surface.

        Raises:
            ValueError: If ``param_vals`` is not a 2D array and
                ``_verify=True``.
        """
        if _verify:
            if param_vals.ndim != 2:
                raise ValueError('Parameter values must be 2D array')
            for lambda1, lambda2, lambda3 in param_vals:
                self._verify_barycentric(lambda1, lambda2, lambda3)

        return _surface_helpers.evaluate_barycentric_multi(
            self._nodes, self._degree, param_vals, self._dimension)

    @staticmethod
    def _verify_cartesian(s, t):
        """Verifies that a point is in the reference triangle.

        I.e., checks that they sum to <= one and are each non-negative.

        Args:
            s (float): Parameter along the reference triangle.
            t (float): Parameter along the reference triangle.

        Raises:
            ValueError: If the point lies outside the reference triangle.
        """
        if s < 0.0 or t < 0.0 or s + t > 1.0:
            raise ValueError('Point lies outside reference triangle', s, t)

    def evaluate_cartesian(self, s, t, _verify=True):
        r"""Compute a point on the surface.

        Evaluates :math:`B\left(1 - s - t, s, t\right)` by calling
        :meth:`evaluate_barycentric`:

        This method acts as a (partial) inverse to :meth:`locate`.

        .. testsetup:: surface-cartesian

           import numpy as np
           import bezier

        .. doctest:: surface-cartesian
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0 , 0.0  ],
           ...     [0.5 , 0.5  ],
           ...     [1.0 , 0.625],
           ...     [0.0 , 0.5  ],
           ...     [0.5 , 0.5  ],
           ...     [0.25, 1.0  ],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
           >>> point = surface.evaluate_cartesian(0.125, 0.375)
           >>> point
           array([[ 0.16015625, 0.44726562]])
           >>> surface.evaluate_barycentric(0.5, 0.125, 0.375)
           array([[ 0.16015625, 0.44726562]])

        Args:
            s (float): Parameter along the reference triangle.
            t (float): Parameter along the reference triangle.
            _verify (Optional[bool]): Indicates if the coordinates should be
                verified inside of the reference triangle. Defaults to
                :data:`True`.

        Returns:
            numpy.ndarray: The point on the surface (as a two dimensional
            NumPy array).
        """
        if _verify:
            self._verify_cartesian(s, t)

        return _surface_helpers.evaluate_barycentric(
            self._nodes, self._degree, 1.0 - s - t, s, t)

    def evaluate_cartesian_multi(self, param_vals, _verify=True):
        r"""Compute multiple points on the surface.

        Assumes ``param_vals`` has two columns of Cartesian coordinates.
        See :meth:`evaluate_cartesian` for more details on how each row of
        parameter values is evaluated.

        .. image:: ../images/surface_evaluate_cartesian_multi.png
           :align: center

        .. doctest:: surface-eval-multi1
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [ 0.0, 0.0],
           ...     [ 2.0, 1.0],
           ...     [-3.0, 2.0],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=1)
           >>> surface
           <Surface (degree=1, dimension=2)>
           >>> param_vals = np.asfortranarray([
           ...     [0.0  , 0.0  ],
           ...     [0.125, 0.625],
           ...     [0.5  , 0.5  ],
           ... ])
           >>> points = surface.evaluate_cartesian_multi(param_vals)
           >>> points
           array([[ 0.   , 0.   ],
                  [-1.625, 1.375],
                  [-0.5  , 1.5  ]])

        .. testcleanup:: surface-eval-multi1

           import make_images
           make_images.surface_evaluate_cartesian_multi(surface, points)

        Args:
            param_vals (numpy.ndarray): Array of parameter values (as a
                ``Nx2`` array).
            _verify (Optional[bool]): Indicates if the coordinates should be
                verified. See :meth:`evaluate_cartesian`. Defaults to
                :data:`True`. Will also double check that ``param_vals``
                is the right shape.

        Returns:
            numpy.ndarray: The points on the surface.

        Raises:
            ValueError: If ``param_vals`` is not a 2D array and
                ``_verify=True``.
        """
        if _verify:
            if param_vals.ndim != 2:
                raise ValueError('Parameter values must be 2D array')
            for s, t in param_vals:
                self._verify_cartesian(s, t)

        return _surface_helpers.evaluate_cartesian_multi(
            self._nodes, self._degree, param_vals, self._dimension)

    def plot(self, pts_per_edge, color=None, ax=None, with_nodes=False):
        """Plot the current surface.

        Args:
            pts_per_edge (int): Number of points to plot per edge.
            color (Optional[Tuple[float, float, float]]): Color as RGB profile.
            ax (Optional[matplotlib.artist.Artist]): matplotlib axis object
                to add plot to.
            with_nodes (Optional[bool]): Determines if the control points
                should be added to the plot. Off by default.

        Returns:
            matplotlib.artist.Artist: The axis containing the plot. This
            may be a newly created axis.

        Raises:
            NotImplementedError: If the surface's dimension is not ``2``.
        """
        if self._dimension != 2:
            raise NotImplementedError('2D is the only supported dimension',
                                      'Current dimension', self._dimension)

        if ax is None:
            ax = _plot_helpers.new_axis()

        _plot_helpers.add_patch(ax, color, pts_per_edge, *self._get_edges())

        if with_nodes:
            ax.plot(self._nodes[:, 0], self._nodes[:, 1],
                    color='black', marker='o', linestyle='None')

        return ax

    def subdivide(self):
        r"""Split the surface into four sub-surfaces.

        Does so by taking the unit triangle (i.e. the domain
        of the surface) and splitting it into four sub-triangles

        .. image:: ../images/surface_subdivide1.png
           :align: center

        Then the surface is re-parameterized via the map to / from the
        given sub-triangles and the unit triangle.

        For example, when a degree two surface is subdivided:

        .. image:: ../images/surface_subdivide2.png
           :align: center

        .. doctest:: surface-subdivide
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [-1.0 , 0.0 ],
           ...     [ 0.5 , 0.5 ],
           ...     [ 2.0 , 0.0 ],
           ...     [ 0.25, 1.75],
           ...     [ 2.0 , 3.0 ],
           ...     [ 0.0 , 4.0 ],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
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
        nodes_a, nodes_b, nodes_c, nodes_d = _surface_helpers.subdivide_nodes(
            self._nodes, self._degree)

        half_width = 0.5 * self._width
        shifted_x = self._base_x + half_width
        shifted_y = self._base_y + half_width
        return (
            Surface(nodes_a, self._degree, base_x=self._base_x,
                    base_y=self._base_y, width=half_width, _copy=False),
            Surface(nodes_b, self._degree, base_x=shifted_x, base_y=shifted_y,
                    width=-half_width, _copy=False),
            Surface(nodes_c, self._degree, base_x=shifted_x,
                    base_y=self._base_y, width=half_width, _copy=False),
            Surface(nodes_d, self._degree, base_x=self._base_x,
                    base_y=shifted_y, width=half_width, _copy=False),
        )

    def _compute_valid(self):
        r"""Determines if the current surface is "valid".

        Does this by checking if the Jacobian of the map from the
        reference triangle is everywhere positive.

        Returns:
            bool: Flag indicating if the current surface is valid.

        Raises:
            NotImplementedError: If the surface is in a dimension other
                than :math:`\mathbf{R}^2`.
            NotImplementedError: If the degree is not 1, 2 or 3.
        """
        if self._dimension != 2:
            raise NotImplementedError(
                'Validity check only implemented in R^2')

        poly_sign = None
        if self._degree == 1:
            # In the linear case, we are only invalid if the points
            # are collinear.
            first_deriv = self._nodes[1:, :] - self._nodes[:-1, :]
            poly_sign = np.sign(np.linalg.det(first_deriv))
        elif self._degree == 2:
            bernstein = _surface_helpers.quadratic_jacobian_polynomial(
                self._nodes)
            poly_sign = _surface_helpers.polynomial_sign(bernstein, 2)
        elif self._degree == 3:
            bernstein = _surface_helpers.cubic_jacobian_polynomial(
                self._nodes)
            poly_sign = _surface_helpers.polynomial_sign(bernstein, 4)
        else:
            raise NotImplementedError(
                'Degrees 1, 2 and 3 only supported at this time')

        return poly_sign == 1

    @property
    def is_valid(self):
        """bool: Flag indicating if the surface is "valid".

        Here, "valid" means there are no self-intersections or
        singularities and the edges are oriented with the interior
        (i.e. a 90 degree rotation of the tangent vector to the left is
        the interior).

        This checks if the Jacobian of the map from the reference
        triangle is everywhere positive. For example, a linear "surface"
        with collinear points is invalid:

        .. image:: ../images/surface_is_valid1.png
           :align: center

        .. doctest:: surface-is-valid1

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.0],
           ...     [1.0, 1.0],
           ...     [2.0, 2.0],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=1)
           >>> surface.is_valid
           False

        .. testcleanup:: surface-is-valid1

           import make_images
           make_images.surface_is_valid1(surface)

        while a quadratic surface with one straight side:

        .. image:: ../images/surface_is_valid2.png
           :align: center

        .. doctest:: surface-is-valid2

           >>> nodes = np.asfortranarray([
           ...     [ 0.0  , 0.0  ],
           ...     [ 0.5  , 0.125],
           ...     [ 1.0  , 0.0  ],
           ...     [-0.125, 0.5  ],
           ...     [ 0.5  , 0.5  ],
           ...     [ 0.0  , 1.0  ],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
           >>> surface.is_valid
           True

        .. testcleanup:: surface-is-valid2

           import make_images
           make_images.surface_is_valid2(surface)

        though not all higher degree surfaces are valid:

        .. image:: ../images/surface_is_valid3.png
           :align: center

        .. doctest:: surface-is-valid3

           >>> nodes = np.asfortranarray([
           ...     [1.0, 0.0],
           ...     [0.0, 0.0],
           ...     [1.0, 1.0],
           ...     [0.0, 0.0],
           ...     [0.0, 0.0],
           ...     [0.0, 1.0],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
           >>> surface.is_valid
           False

        .. testcleanup:: surface-is-valid3

           import make_images
           make_images.surface_is_valid3(surface)
        """
        if self._is_valid is None:
            self._is_valid = self._compute_valid()
        return self._is_valid

    @property
    def __dict__(self):
        """dict: Dictionary of current surface's property namespace.

        This is just a stand-in property for the usual ``__dict__``. This
        class defines ``__slots__`` so by default would not provide a
        ``__dict__``.

        This also means that the current object can't be modified by the
        returned dictionary.
        """
        return {
            '_dimension': self._dimension,
            '_nodes': self._nodes,
            '_degree': self._degree,
            '_base_x': self._base_x,
            '_base_y': self._base_y,
            '_width': self._width,
            '_area': self._area,
            '_edges': self._edges,
            '_is_valid': self._is_valid,
        }

    def locate(self, point, _verify=True):
        r"""Find a point on the current surface.

        Solves for :math:`s` and :math:`t` in :math:`B(s, t) = p`.

        This method acts as a (partial) inverse to :meth:`evaluate_cartesian`.

        .. warning::

           A unique solution is only guaranteed if the current surface is
           valid. This code assumes a valid surface, but doesn't check.

        .. image:: ../images/surface_locate.png
           :align: center

        .. doctest:: surface-locate

           >>> nodes = np.asfortranarray([
           ...     [0.0 ,  0.0 ],
           ...     [0.5 , -0.25],
           ...     [1.0 ,  0.0 ],
           ...     [0.25,  0.5 ],
           ...     [0.75,  0.75],
           ...     [0.0 ,  1.0 ],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
           >>> point = np.asfortranarray([[0.59375, 0.25]])
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
            _verify (Optional[bool]): Indicates if extra caution should be
                used to verify assumptions about the inputs. Can be
                disabled to speed up execution time. Defaults to :data:`True`.

        Returns:
            Optional[Tuple[float, float]]: The :math:`s` and :math:`t`
            values corresponding to ``point`` or :data:`None` if the point
            is not on the surface.

        Raises:
            NotImplementedError: If the surface isn't in :math:`\mathbf{R}^2`.
            ValueError: If the dimension of the ``point`` doesn't match the
                dimension of the current surface.
        """
        if _verify:
            if self._dimension != 2:
                raise NotImplementedError('Only 2D surfaces supported.')

            if point.shape != (1, self._dimension):
                point_dimensions = ' x '.join(
                    str(dimension) for dimension in point.shape)
                msg = _LOCATE_ERROR_TEMPLATE.format(
                    self._dimension, self._dimension, point, point_dimensions)
                raise ValueError(msg)

        return _surface_intersection.locate_point(
            self._nodes, self._degree, point[0, 0], point[0, 1])

    def intersect(self, other, strategy=_STRATEGY.GEOMETRIC, _verify=True):
        """Find the common intersection with another surface.

        Args:
            other (Surface): Other surface to intersect with.
            strategy (Optional[~bezier.curve.IntersectionStrategy]): The
                intersection algorithm to use. Defaults to geometric.
            _verify (Optional[bool]): Indicates if extra caution should be
                used to verify assumptions about the algorithm as it
                proceeds. Can be disabled to speed up execution time.
                Defaults to :data:`True`.

        Returns:
            List[~bezier.curved_polygon.CurvedPolygon]: List of
            intersections (possibly empty).

        Raises:
            TypeError: If ``other`` is not a surface (and ``_verify=True``).
            NotImplementedError: If at least one of the surfaces
                isn't two-dimensional (and ``_verify=True``).
        """
        if _verify:
            if not isinstance(other, Surface):
                raise TypeError('Can only intersect with another surface',
                                'Received', other)
            if self._dimension != 2 or other._dimension != 2:
                raise NotImplementedError(
                    'Intersection only implemented in 2D')

        bbox_int = _geometric_intersection.bbox_intersect(
            self._nodes, other._nodes)
        if bbox_int != _INTERSECTION_T:
            return []

        # We need **all** pairs of edges.
        edges1 = self._get_edges()
        edges2 = other._get_edges()  # pylint: disable=protected-access
        intersections, duplicates = _surface_intersections(
            edges1, edges2, strategy)

        if _verify:
            _surface_helpers.verify_duplicates(duplicates, intersections)
        return _surface_helpers.combine_intersections(
            intersections, self, edges1, other, edges2)

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

        where we define, for example, :math:`v_{i, j, k - 1} = 0`
        if :math:`k = 0`.

        .. image:: ../images/surface_elevate.png
           :align: center

        .. doctest:: surface-elevate
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0 , 0.0 ],
           ...     [1.5 , 0.0 ],
           ...     [3.0 , 0.0 ],
           ...     [0.75, 1.5 ],
           ...     [2.25, 2.25],
           ...     [0.0 , 3.0 ],
           ... ])
           >>> surface = bezier.Surface(nodes, degree=2)
           >>> elevated = surface.elevate()
           >>> elevated
           <Surface (degree=3, dimension=2)>
           >>> elevated.nodes
           array([[ 0.  , 0.  ],
                  [ 1.  , 0.  ],
                  [ 2.  , 0.  ],
                  [ 3.  , 0.  ],
                  [ 0.5 , 1.  ],
                  [ 1.5 , 1.25],
                  [ 2.5 , 1.5 ],
                  [ 0.5 , 2.  ],
                  [ 1.5 , 2.5 ],
                  [ 0.  , 3.  ]])

        .. testcleanup:: surface-elevate

           import make_images
           make_images.surface_elevate(surface, elevated)

        Returns:
            Surface: The degree-elevated surface.
        """
        num_nodes, _ = self._nodes.shape
        # (d + 1)(d + 2)/2 --> (d + 2)(d + 3)/2
        num_new = num_nodes + self._degree + 2
        new_nodes = np.zeros((num_new, self._dimension), order='F')

        # NOTE: We start from the index triples (i, j, k) for the current
        #       nodes and map them onto (i + 1, j, k), etc. This index
        #       tracking is also done in :func:`.de_casteljau_one_round`.
        index = 0
        # parent_i1 = index + k
        # parent_i2 = index + k + 1
        # parent_i3 = index + degree + 2
        parent_i1 = 0
        parent_i2 = 1
        parent_i3 = self._degree + 2
        for k in six.moves.xrange(self._degree + 1):
            for j in six.moves.xrange(self._degree + 1 - k):
                i = self._degree - j - k
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
        denominator = self._degree + 1.0
        new_nodes /= denominator

        return Surface(new_nodes, self._degree + 1, _copy=False)


def _surface_intersections(edges1, edges2, strategy):
    """Find all intersections among edges of two surfaces.

    This treats intersections which have ``s == 1.0`` or ``t == 1.0``
    as duplicates. The duplicates may be checked by the caller, e.g.
    by :func:`~bezier._surface_helpers.verify_duplicates`.

    Args:
        edges1 (Tuple[.Curve, .Curve, .Curve]): The three edges
            of the first surface being intersected.
        edges2 (Tuple[.Curve, .Curve, .Curve]): The three edges
            of the second surface being intersected.
        strategy (Optional[~bezier.curve.IntersectionStrategy]): The
            intersection algorithm to use. Defaults to geometric.

    Returns:
        Tuple[list, list]: Pair of lists (both of which could be empty)

        * The actual "unique" :class:`Intersection`-s
        * Duplicate :class:`Intersection`-s encountered (these will be
          corner intersections)

    Raises:
        ValueError: If ``strategy`` is not a valid
            :attr:`.IntersectionStrategy`.
    """
    if strategy == _STRATEGY.GEOMETRIC:
        all_intersections = _geometric_intersection.all_intersections
    elif strategy == _STRATEGY.ALGEBRAIC:
        all_intersections = _algebraic_intersection.all_intersections
    else:
        raise ValueError('Unexpected strategy.', strategy)

    intersections = []
    duplicates = []
    for index1, edge1 in enumerate(edges1):
        for index2, edge2 in enumerate(edges2):
            st_vals = all_intersections(edge1._nodes, edge2._nodes)
            for s, t in st_vals:
                _add_intersection(
                    index1, s, index2, t, edges1, edges2,
                    duplicates, intersections)

    return intersections, duplicates


def _add_intersection(  # pylint: disable=too-many-arguments
        index1, s, index2, t, edges1, edges2, duplicates, intersections):
    """Create an :class:`Intersection` and append.

    The intersection will be classified as either a duplicate or a valid
    intersection and appended to one of ``duplicates`` or ``intersections``
    depending on that classification.

    Args:
        index1 (int): The index (among 0, 1, 2) of the first edge in the
            intersection.
        s (float): The parameter along the first curve of the intersection.
        index2 (int): The index (among 0, 1, 2) of the second edge in the
            intersection.
        t (float): The parameter along the second curve of the intersection.
        edges1 (Tuple[.Curve, .Curve, .Curve]): The three edges
            of the first surface being intersected.
        edges2 (Tuple[.Curve, .Curve, .Curve]): The three edges
            of the second surface being intersected.
        duplicates (List[.Intersection]): List of duplicate intersections.
        intersections (List[.Intersection]): List of "accepted" (i.e.
            non-duplicate) intersections.
    """
    edge_end, intersection_args = _surface_helpers.handle_ends(
        index1, s, index2, t)
    if edge_end:
        intersection = _intersection_helpers.Intersection(
            *intersection_args)
        duplicates.append(intersection)
    else:
        intersection = _intersection_helpers.Intersection(
            index1, s, index2, t)
        # Classify the intersection.
        interior = _surface_helpers.classify_intersection(
            intersection, edges1, edges2)
        intersection.interior_curve = interior
        intersections.append(intersection)
