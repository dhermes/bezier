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

"""Helper for B |eacute| zier Triangles.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

.. testsetup:: *

   import numpy as np
   import bezier
"""

import numpy as np

from bezier import _base
from bezier import _curve_helpers
from bezier import _plot_helpers
from bezier import _symbolic
from bezier import _triangle_helpers
from bezier import _triangle_intersection
from bezier import curve as _curve_mod
from bezier import curved_polygon
from bezier.hazmat import helpers as _py_helpers
from bezier.hazmat import intersection_helpers
from bezier.hazmat import triangle_helpers as _py_triangle_helpers
from bezier.hazmat import triangle_intersection as _py_triangle_intersection


_SIGN = np.sign  # pylint: disable=no-member
_LOCATE_ERROR_TEMPLATE = (
    "Dimension mismatch: This triangle is {:d}-dimensional, "
    "so the point should be a {:d} x 1 NumPy array. "
    "Instead the point {} has dimensions {}."
)
_STRATEGY = intersection_helpers.IntersectionStrategy


class Triangle(_base.Base):
    r"""Represents a B |eacute| zier `triangle`_.

    .. _triangle: https://en.wikipedia.org/wiki/B%C3%A9zier_triangle
    .. _unit simplex:
        https://en.wikipedia.org/wiki/Simplex#The_standard_simplex
    .. _barycentric coordinates:
        https://en.wikipedia.org/wiki/Barycentric_coordinate_system

    We define a B |eacute| zier triangle as a mapping from the
    `unit simplex`_ in :math:`\mathbf{R}^2` (i.e. the unit triangle) onto a
    triangle in an arbitrary dimension. We use `barycentric coordinates`_

    .. math::

       \lambda_1 = 1 - s - t, \lambda_2 = s, \lambda_3 = t

    for points in the unit triangle
    :math:`\left\{(s, t) \mid 0 \leq s, t, s + t \leq 1\right\}`:

    .. image:: ../../images/unit_triangle.png
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
              v_{1,0,0} & v_{0,1,0} & v_{0,0,1} \end{array}\right]

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
              v_{0,1,1} & v_{0,0,2} \end{array}\right]

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
              v_{0,1,2} & v_{0,0,3} \end{array}\right]

       and so on.

       The index formula

       .. math::

          j + \frac{k}{2} \left(2 (i + j) + k + 3\right)

       can be used to map a triple :math:`(i, j, k)` onto the
       corresponding linear index, but it is not particularly insightful
       or useful.

    .. image:: ../../images/triangle_constructor.png
       :align: center

    .. doctest:: triangle-constructor

       >>> import bezier
       >>> import numpy as np
       >>> nodes = np.asfortranarray([
       ...     [0.0, 0.5, 1.0 , 0.125, 0.375, 0.25],
       ...     [0.0, 0.0, 0.25, 0.5  , 0.375, 1.0 ],
       ... ])
       >>> triangle = bezier.Triangle(nodes, degree=2)
       >>> triangle
       <Triangle (degree=2, dimension=2)>

    .. testcleanup:: triangle-constructor

       import make_images
       make_images.unit_triangle()
       make_images.triangle_constructor(triangle)

    Args:
        nodes (Sequence[Sequence[numbers.Number]]): The nodes in the triangle.
            Must be convertible to a 2D NumPy array of floating point values,
            where the columns represent each node while the rows are the
            dimension of the ambient space.
        degree (int): The degree of the triangle. This is assumed to
            correctly correspond to the number of ``nodes``. Use
            :meth:`from_nodes` if the degree has not yet been computed.
        copy (bool): Flag indicating if the nodes should be copied before
            being stored. Defaults to :data:`True` since callers may
            freely mutate ``nodes`` after passing in.
        verify (bool): Flag indicating if the degree should be verified against
            the number of nodes. Defaults to :data:`True`.
    """

    __slots__ = (
        "_degree",  # From constructor
        "_edges",  # Empty default
    )

    def __init__(self, nodes, degree, *, copy=True, verify=True):
        super().__init__(nodes, copy=copy)
        self._degree = degree
        self._edges = None
        self._verify_degree(verify)

    @classmethod
    def from_nodes(cls, nodes, copy=True):
        """Create a :class:`.Triangle` from nodes.

        Computes the ``degree`` based on the shape of ``nodes``.

        Args:
            nodes (Sequence[Sequence[numbers.Number]]): The nodes in the
                triangle. Must be convertible to a 2D NumPy array of floating
                point values, where the columns represent each node while the
                rows are the dimension of the ambient space.
            copy (bool): Flag indicating if the nodes should be copied before
                being stored. Defaults to :data:`True` since callers may
                freely mutate ``nodes`` after passing in.

        Returns:
            Triangle: The constructed triangle.
        """
        nodes_np = _base.sequence_to_array(nodes)
        _, num_nodes = nodes_np.shape
        degree = cls._get_degree(num_nodes)
        # NOTE: **Explicitly** verify because ``_get_degree`` does not.
        return cls(nodes_np, degree, copy=copy, verify=True)

    @staticmethod
    def _get_degree(num_nodes):
        """Get the degree of the current triangle.

        .. note::

            If ``num_nodes`` isn't a triangular number, no degree can be
            correct so the return value will be invalid. Callers should use
            ``verify`` in the constructor to ensure correctness.

        Args:
            num_nodes (int): The number of control points for a
                B |eacute| zier triangle.

        Returns:
            int: The degree :math:`d` such that :math:`(d + 1)(d + 2)/2`
            equals ``num_nodes``.
        """
        # 8 * num_nodes = 4(d + 1)(d + 2)
        #               = 4d^2 + 12d + 8
        #               = (2d + 3)^2 - 1
        d_float = 0.5 * (np.sqrt(8.0 * num_nodes + 1.0) - 3.0)
        return int(np.round(d_float))

    def _verify_degree(self, verify):
        """Verify that the number of nodes matches the degree.

        Args:
            verify (bool): Flag indicating if the degree should be verified
                against the number of nodes.

        Raises:
            ValueError: If ``verify`` is :data:`True` and the number of nodes
                does not match the degree.
        """
        if not verify:
            return

        _, num_nodes = self._nodes.shape
        twice_expected_nodes = (self._degree + 1) * (self._degree + 2)
        # Avoid rounding by division by 2.
        if twice_expected_nodes == 2 * num_nodes:
            return

        msg = (
            f"A degree {self._degree} triangle should have "
            f"{0.5 * twice_expected_nodes:g} nodes, not {num_nodes}."
        )
        raise ValueError(msg)

    @property
    def area(self):
        r"""The area of the current triangle.

        For triangles in :math:`\mathbf{R}^2`, this computes the area via
        Green's theorem. Using the vector field :math:`\mathbf{F} =
        \left[-y, x\right]^T`, since :math:`\partial_x(x) - \partial_y(-y) = 2`
        Green's theorem says twice the area is equal to

        .. math::

           \int_{B\left(\mathcal{U}\right)} 2 \, d\mathbf{x} =
           \int_{\partial B\left(\mathcal{U}\right)} -y \, dx + x \, dy.

        This relies on the assumption that the current triangle is valid, which
        implies that the image of the unit triangle under the B |eacute| zier
        map --- :math:`B\left(\mathcal{U}\right)`  --- has the edges of the
        triangle as its boundary.

        Note that for a given edge :math:`C(r)` with control points
        :math:`x_j, y_j`, the integral can be simplified:

        .. math::

            \int_C -y \, dx + x \, dy = \int_0^1 (x y' - y x') \, dr
                = \sum_{i < j} (x_i y_j - y_i x_j) \int_0^1 b_{i, d}
                b'_{j, d} \, dr

        where :math:`b_{i, d}, b_{j, d}` are Bernstein basis polynomials.

        Returns:
            float: The area of the current triangle.

        Raises:
            NotImplementedError: If the current triangle isn't in
                :math:`\mathbf{R}^2`.
        """
        if self._dimension != 2:
            raise NotImplementedError(
                "2D is the only supported dimension",
                "Current dimension",
                self._dimension,
            )

        edge1, edge2, edge3 = self._get_edges()
        return _triangle_helpers.compute_area(
            (edge1._nodes, edge2._nodes, edge3._nodes)
        )

    def _compute_edges(self):
        """Compute the edges of the current triangle.

        Returns:
            Tuple[~curve.Curve, ~curve.Curve, ~curve.Curve]: The edges of
            the triangle.
        """
        nodes1, nodes2, nodes3 = _triangle_helpers.compute_edge_nodes(
            self._nodes, self._degree
        )
        edge1 = _curve_mod.Curve(
            nodes1, self._degree, copy=False, verify=False
        )
        edge2 = _curve_mod.Curve(
            nodes2, self._degree, copy=False, verify=False
        )
        edge3 = _curve_mod.Curve(
            nodes3, self._degree, copy=False, verify=False
        )
        return edge1, edge2, edge3

    def _get_edges(self):
        """Get the edges for the current triangle.

        If they haven't been computed yet, first compute and store them.

        This is provided as a means for internal calls to get the edges
        without copying (since :attr:`.edges` copies before giving to
        a user to keep the stored data immutable).

        Returns:
            Tuple[~bezier.curve.Curve, ~bezier.curve.Curve, \
                  ~bezier.curve.Curve]: The edges of
            the triangle.
        """
        if self._edges is None:
            self._edges = self._compute_edges()
        return self._edges

    @property
    def edges(self):
        """The edges of the triangle.

        .. doctest:: triangle-edges
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0,  0.5   , 1.0, 0.1875, 0.625, 0.0],
           ...     [0.0, -0.1875, 0.0, 0.5   , 0.625, 1.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> edge1, _, _ = triangle.edges
           >>> edge1
           <Curve (degree=2, dimension=2)>
           >>> edge1.nodes
           array([[ 0. ,  0.5   , 1. ],
                  [ 0. , -0.1875, 0. ]])

        Returns:
            Tuple[~bezier.curve.Curve, ~bezier.curve.Curve, \
                  ~bezier.curve.Curve]: The edges of
            the triangle.
        """
        edge1, edge2, edge3 = self._get_edges()
        # NOTE: It is crucial that we return copies here. Since the edges
        #       are cached, if they were mutable, callers could
        #       inadvertently mutate the cached value.
        edge1 = edge1.copy()
        edge2 = edge2.copy()
        edge3 = edge3.copy()
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
        if not np.allclose(weights_total, 1.0, atol=0.0):
            raise ValueError(
                "Weights do not sum to 1", lambda1, lambda2, lambda3
            )

        if lambda1 < 0.0 or lambda2 < 0.0 or lambda3 < 0.0:
            raise ValueError(
                "Weights must be positive", lambda1, lambda2, lambda3
            )

    def evaluate_barycentric(self, lambda1, lambda2, lambda3, verify=True):
        r"""Compute a point on the triangle.

        Evaluates :math:`B\left(\lambda_1, \lambda_2, \lambda_3\right)`.

        .. image:: ../../images/triangle_evaluate_barycentric.png
           :align: center

        .. testsetup:: triangle-barycentric, triangle-barycentric-fail1,
                       triangle-barycentric-fail2,
                       triangle-barycentric-no-verify

           import numpy as np
           import bezier
           nodes = np.asfortranarray([
               [0.0, 0.5, 1.0 , 0.125, 0.375, 0.25],
               [0.0, 0.0, 0.25, 0.5  , 0.375, 1.0 ],
           ])
           triangle = bezier.Triangle(nodes, degree=2)

        .. doctest:: triangle-barycentric
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.5, 1.0 , 0.125, 0.375, 0.25],
           ...     [0.0, 0.0, 0.25, 0.5  , 0.375, 1.0 ],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> point = triangle.evaluate_barycentric(0.125, 0.125, 0.75)
           >>> point
           array([[0.265625  ],
                  [0.73046875]])

        .. testcleanup:: triangle-barycentric

           import make_images
           make_images.triangle_evaluate_barycentric(triangle, point)

        However, this can't be used for points **outside** the
        reference triangle:

        .. doctest:: triangle-barycentric-fail1

           >>> triangle.evaluate_barycentric(-0.25, 0.75, 0.5)
           Traceback (most recent call last):
             ...
           ValueError: ('Weights must be positive', -0.25, 0.75, 0.5)

        or for non-barycentric coordinates;

        .. doctest:: triangle-barycentric-fail2

           >>> triangle.evaluate_barycentric(0.25, 0.25, 0.25)
           Traceback (most recent call last):
             ...
           ValueError: ('Weights do not sum to 1', 0.25, 0.25, 0.25)

        However, these "invalid" inputs can be used if ``verify`` is
        :data:`False`.

        .. doctest:: triangle-barycentric-no-verify
           :options: +NORMALIZE_WHITESPACE

           >>> triangle.evaluate_barycentric(-0.25, 0.75, 0.5, verify=False)
           array([[0.6875  ],
                  [0.546875]])
           >>> triangle.evaluate_barycentric(0.25, 0.25, 0.25, verify=False)
           array([[0.203125],
                  [0.1875  ]])

        Args:
            lambda1 (float): Parameter along the reference triangle.
            lambda2 (float): Parameter along the reference triangle.
            lambda3 (float): Parameter along the reference triangle.
            verify (Optional[bool]): Indicates if the barycentric coordinates
                should be verified as summing to one and all non-negative (i.e.
                verified as barycentric). Can either be used to evaluate at
                points outside the domain, or to save time when the caller
                already knows the input is verified. Defaults to :data:`True`.

        Returns:
            numpy.ndarray: The point on the triangle (as a two dimensional
            NumPy array with a single column).

        Raises:
            ValueError: If the weights are not valid barycentric
                coordinates, i.e. they don't sum to ``1``. (Won't raise if
                ``verify=False``.)
            ValueError: If some weights are negative. (Won't raise if
                ``verify=False``.)
        """
        if verify:
            self._verify_barycentric(lambda1, lambda2, lambda3)
        return _triangle_helpers.evaluate_barycentric(
            self._nodes, self._degree, lambda1, lambda2, lambda3
        )

    def evaluate_barycentric_multi(self, param_vals, verify=True):
        r"""Compute multiple points on the triangle.

        Assumes ``param_vals`` has three columns of barycentric coordinates.
        See :meth:`evaluate_barycentric` for more details on how each row of
        parameter values is evaluated.

        .. image:: ../../images/triangle_evaluate_barycentric_multi.png
           :align: center

        .. doctest:: triangle-eval-multi2
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 1.0 , 2.0, -1.5, -0.5, -3.0],
           ...     [0.0, 0.75, 1.0,  1.0,  1.5,  2.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> triangle
           <Triangle (degree=2, dimension=2)>
           >>> param_vals = np.asfortranarray([
           ...     [0.   , 0.25, 0.75 ],
           ...     [1.   , 0.  , 0.   ],
           ...     [0.25 , 0.5 , 0.25 ],
           ...     [0.375, 0.25, 0.375],
           ... ])
           >>> points = triangle.evaluate_barycentric_multi(param_vals)
           >>> points
           array([[-1.75 , 0. , 0.25   , -0.625   ],
                  [ 1.75 , 0. , 1.0625 ,  1.046875]])

        .. testcleanup:: triangle-eval-multi2

           import make_images
           make_images.triangle_evaluate_barycentric_multi(triangle, points)

        Args:
            param_vals (numpy.ndarray): Array of parameter values (as a
                ``N x 3`` array).
            verify (Optional[bool]): Indicates if the coordinates should be
                verified. See :meth:`evaluate_barycentric`. Defaults to
                :data:`True`. Will also double check that ``param_vals``
                is the right shape.

        Returns:
            numpy.ndarray: The points on the triangle.

        Raises:
            ValueError: If ``param_vals`` is not a 2D array and
                ``verify=True``.
        """
        if verify:
            if param_vals.ndim != 2:
                raise ValueError("Parameter values must be 2D array")

            for lambda1, lambda2, lambda3 in param_vals:
                self._verify_barycentric(lambda1, lambda2, lambda3)
        return _triangle_helpers.evaluate_barycentric_multi(
            self._nodes, self._degree, param_vals, self._dimension
        )

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
            raise ValueError("Point lies outside reference triangle", s, t)

    def evaluate_cartesian(self, s, t, verify=True):
        r"""Compute a point on the triangle.

        Evaluates :math:`B\left(1 - s - t, s, t\right)` by calling
        :meth:`evaluate_barycentric`:

        This method acts as a (partial) inverse to :meth:`locate`.

        .. testsetup:: triangle-cartesian

           import numpy as np
           import bezier

        .. doctest:: triangle-cartesian
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.5, 1.0  , 0.0, 0.5, 0.25],
           ...     [0.0, 0.5, 0.625, 0.5, 0.5, 1.0 ],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> point = triangle.evaluate_cartesian(0.125, 0.375)
           >>> point
           array([[0.16015625],
                  [0.44726562]])
           >>> triangle.evaluate_barycentric(0.5, 0.125, 0.375)
           array([[0.16015625],
                  [0.44726562]])

        Args:
            s (float): Parameter along the reference triangle.
            t (float): Parameter along the reference triangle.
            verify (Optional[bool]): Indicates if the coordinates should be
                verified inside of the reference triangle. Defaults to
                :data:`True`.

        Returns:
            numpy.ndarray: The point on the triangle (as a two dimensional
            NumPy array).
        """
        if verify:
            self._verify_cartesian(s, t)
        return _triangle_helpers.evaluate_barycentric(
            self._nodes, self._degree, 1.0 - s - t, s, t
        )

    def evaluate_cartesian_multi(self, param_vals, verify=True):
        r"""Compute multiple points on the triangle.

        Assumes ``param_vals`` has two columns of Cartesian coordinates.
        See :meth:`evaluate_cartesian` for more details on how each row of
        parameter values is evaluated.

        .. image:: ../../images/triangle_evaluate_cartesian_multi.png
           :align: center

        .. doctest:: triangle-eval-multi1
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 2.0, -3.0],
           ...     [0.0, 1.0,  2.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=1)
           >>> triangle
           <Triangle (degree=1, dimension=2)>
           >>> param_vals = np.asfortranarray([
           ...     [0.0  , 0.0  ],
           ...     [0.125, 0.625],
           ...     [0.5  , 0.5  ],
           ... ])
           >>> points = triangle.evaluate_cartesian_multi(param_vals)
           >>> points
           array([[ 0. , -1.625, -0.5 ],
                  [ 0. ,  1.375,  1.5 ]])

        .. testcleanup:: triangle-eval-multi1

           import make_images
           make_images.triangle_evaluate_cartesian_multi(triangle, points)

        Args:
            param_vals (numpy.ndarray): Array of parameter values (as a
                ``N x 2`` array).
            verify (Optional[bool]): Indicates if the coordinates should be
                verified. See :meth:`evaluate_cartesian`. Defaults to
                :data:`True`. Will also double check that ``param_vals``
                is the right shape.

        Returns:
            numpy.ndarray: The points on the triangle.

        Raises:
            ValueError: If ``param_vals`` is not a 2D array and
                ``verify=True``.
        """
        if verify:
            if param_vals.ndim != 2:
                raise ValueError("Parameter values must be 2D array")

            for s, t in param_vals:
                self._verify_cartesian(s, t)
        return _triangle_helpers.evaluate_cartesian_multi(
            self._nodes, self._degree, param_vals, self._dimension
        )

    def plot(
        self, pts_per_edge, color=None, ax=None, with_nodes=False, alpha=0.625
    ):
        """Plot the current triangle.

        Args:
            pts_per_edge (int): Number of points to plot per edge.
            color (Optional[Tuple[float, float, float]]): Color as RGB profile.
            ax (Optional[matplotlib.artist.Artist]): matplotlib axis object
                to add plot to.
            with_nodes (Optional[bool]): Determines if the control points
                should be added to the plot. Off by default.
            alpha (Optional[float]): Alpha value of patch center, between 0 and
                1 inclusive.

        Returns:
            matplotlib.artist.Artist: The axis containing the plot. This
            may be a newly created axis.

        Raises:
            NotImplementedError: If the triangle's dimension is not ``2``.
        """
        if self._dimension != 2:
            raise NotImplementedError(
                "2D is the only supported dimension",
                "Current dimension",
                self._dimension,
            )

        if ax is None:
            ax = _plot_helpers.new_axis()
        _plot_helpers.add_patch(
            ax, color, pts_per_edge, *self._get_edges(), alpha=alpha
        )
        if with_nodes:
            ax.plot(
                self._nodes[0, :],
                self._nodes[1, :],
                color="black",
                marker="o",
                linestyle="None",
            )
        return ax

    def subdivide(self):
        r"""Split the triangle into four sub-triangles.

        Does so by taking the unit triangle (i.e. the domain
        of the triangle) and splitting it into four sub-triangles

        .. image:: ../../images/triangle_subdivide1.png
           :align: center

        Then the triangle is re-parameterized via the map to / from the
        given sub-triangles and the unit triangle.

        For example, when a degree two triangle is subdivided:

        .. image:: ../../images/triangle_subdivide2.png
           :align: center

        .. doctest:: triangle-subdivide
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [-1.0, 0.5, 2.0, 0.25, 2.0, 0.0],
           ...     [ 0.0, 0.5, 0.0, 1.75, 3.0, 4.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> _, sub_triangle_b, _, _ = triangle.subdivide()
           >>> sub_triangle_b
           <Triangle (degree=2, dimension=2)>
           >>> sub_triangle_b.nodes
           array([[ 1.5 ,  0.6875, -0.125 , 1.1875, 0.4375, 0.5  ],
                  [ 2.5 ,  2.3125,  1.875 , 1.3125, 1.3125, 0.25 ]])

        .. testcleanup:: triangle-subdivide

           import make_images
           make_images.triangle_subdivide1()
           make_images.triangle_subdivide2(triangle, sub_triangle_b)

        Returns:
            Tuple[Triangle, Triangle, Triangle, Triangle]: The lower left,
            central, lower right and upper left sub-triangles (in that order).
        """
        (
            nodes_a,
            nodes_b,
            nodes_c,
            nodes_d,
        ) = _triangle_helpers.subdivide_nodes(self._nodes, self._degree)
        return (
            Triangle(nodes_a, self._degree, copy=False, verify=False),
            Triangle(nodes_b, self._degree, copy=False, verify=False),
            Triangle(nodes_c, self._degree, copy=False, verify=False),
            Triangle(nodes_d, self._degree, copy=False, verify=False),
        )

    def _compute_valid(self):
        r"""Determines if the current triangle is "valid".

        Does this by checking if the Jacobian of the map from the
        reference triangle is everywhere positive.

        Returns:
            bool: Flag indicating if the current triangle is valid.

        Raises:
            NotImplementedError: If the triangle is in a dimension other
                than :math:`\mathbf{R}^2`.
            .UnsupportedDegree: If the degree is not 1, 2 or 3.
        """
        if self._dimension != 2:
            raise NotImplementedError("Validity check only implemented in R^2")

        poly_sign = None
        if self._degree == 1:
            # In the linear case, we are only invalid if the points
            # are collinear.
            first_deriv = self._nodes[:, 1:] - self._nodes[:, :-1]
            # pylint: disable=assignment-from-no-return
            poly_sign = _SIGN(np.linalg.det(first_deriv))
            # pylint: enable=assignment-from-no-return
        elif self._degree == 2:
            bernstein = _py_triangle_helpers.quadratic_jacobian_polynomial(
                self._nodes
            )
            poly_sign = _py_triangle_helpers.polynomial_sign(bernstein, 2)
        elif self._degree == 3:
            bernstein = _py_triangle_helpers.cubic_jacobian_polynomial(
                self._nodes
            )
            poly_sign = _py_triangle_helpers.polynomial_sign(bernstein, 4)
        else:
            raise _py_helpers.UnsupportedDegree(
                self._degree, supported=(1, 2, 3)
            )

        return poly_sign == 1

    @property
    def is_valid(self):
        """bool: Flag indicating if the triangle is "valid".

        Here, "valid" means there are no self-intersections or
        singularities and the edges are oriented with the interior
        (i.e. a 90 degree rotation of the tangent vector to the left is
        the interior).

        This checks if the Jacobian of the map from the reference
        triangle is everywhere positive. For example, a linear "triangle"
        with collinear points is invalid:

        .. image:: ../../images/triangle_is_valid1.png
           :align: center

        .. doctest:: triangle-is-valid1

           >>> nodes = np.asfortranarray([
           ...     [0.0, 1.0, 2.0],
           ...     [0.0, 1.0, 2.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=1)
           >>> bool(triangle.is_valid)
           False

        .. testcleanup:: triangle-is-valid1

           import make_images
           make_images.triangle_is_valid1(triangle)

        while a quadratic triangle with one straight side:

        .. image:: ../../images/triangle_is_valid2.png
           :align: center

        .. doctest:: triangle-is-valid2

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.5  , 1.0, -0.125, 0.5, 0.0],
           ...     [0.0, 0.125, 0.0,  0.5  , 0.5, 1.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> bool(triangle.is_valid)
           True

        .. testcleanup:: triangle-is-valid2

           import make_images
           make_images.triangle_is_valid2(triangle)

        though not all higher degree triangles are valid:

        .. image:: ../../images/triangle_is_valid3.png
           :align: center

        .. doctest:: triangle-is-valid3

           >>> nodes = np.asfortranarray([
           ...     [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
           ...     [0.0, 0.0, 1.0, 0.0, 0.0, 1.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> triangle.is_valid
           False

        .. testcleanup:: triangle-is-valid3

           import make_images
           make_images.triangle_is_valid3(triangle)
        """
        return self._compute_valid()

    @property
    def __dict__(self):
        """dict: Dictionary of current triangle's property namespace.

        This is just a stand-in property for the usual ``__dict__``. This
        class defines ``__slots__`` so by default would not provide a
        ``__dict__``.

        This also means that the current object can't be modified by the
        returned dictionary.
        """
        return {
            "_dimension": self._dimension,
            "_nodes": self._nodes,
            "_degree": self._degree,
            "_edges": self._edges,
        }

    def locate(self, point, verify=True):
        r"""Find a point on the current triangle.

        Solves for :math:`s` and :math:`t` in :math:`B(s, t) = p`.

        This method acts as a (partial) inverse to :meth:`evaluate_cartesian`.

        .. warning::

           A unique solution is only guaranteed if the current triangle is
           valid. This code assumes a valid triangle, but doesn't check.

        .. image:: ../../images/triangle_locate.png
           :align: center

        .. doctest:: triangle-locate

           >>> nodes = np.asfortranarray([
           ...     [0.0,  0.5 , 1.0, 0.25, 0.75, 0.0],
           ...     [0.0, -0.25, 0.0, 0.5 , 0.75, 1.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> point = np.asfortranarray([
           ...     [0.59375],
           ...     [0.25   ],
           ... ])
           >>> s, t = triangle.locate(point)
           >>> s
           0.5
           >>> t
           0.25

        .. testcleanup:: triangle-locate

           import make_images
           make_images.triangle_locate(triangle, point)

        Args:
            point (numpy.ndarray): A (``D x 1``) point on the triangle,
                where :math:`D` is the dimension of the triangle.
            verify (Optional[bool]): Indicates if extra caution should be
                used to verify assumptions about the inputs. Can be
                disabled to speed up execution time. Defaults to :data:`True`.

        Returns:
            Optional[Tuple[float, float]]: The :math:`s` and :math:`t`
            values corresponding to ``point`` or :data:`None` if the point
            is not on the triangle.

        Raises:
            NotImplementedError: If the triangle isn't in :math:`\mathbf{R}^2`.
            ValueError: If the dimension of the ``point`` doesn't match the
                dimension of the current triangle.
        """
        if verify:
            if self._dimension != 2:
                raise NotImplementedError("Only 2D triangles supported.")

            if point.shape != (self._dimension, 1):
                point_dimensions = " x ".join(
                    str(dimension) for dimension in point.shape
                )
                msg = _LOCATE_ERROR_TEMPLATE.format(
                    self._dimension, self._dimension, point, point_dimensions
                )
                raise ValueError(msg)

        return _triangle_intersection.locate_point(
            self._nodes, self._degree, point[0, 0], point[1, 0]
        )

    def intersect(self, other, strategy=_STRATEGY.GEOMETRIC, verify=True):
        """Find the common intersection with another triangle.

        Args:
            other (Triangle): Other triangle to intersect with.
            strategy (Optional[ \
                ~bezier.hazmat.intersection_helpers.IntersectionStrategy]): The
                intersection algorithm to use. Defaults to geometric.
            verify (Optional[bool]): Indicates if extra caution should be
                used to verify assumptions about the algorithm as it
                proceeds. Can be disabled to speed up execution time.
                Defaults to :data:`True`.

        Returns:
            List[Union[~bezier.curved_polygon.CurvedPolygon, \
            ~bezier.triangle.Triangle]]: List of intersections (possibly
            empty).

        Raises:
            TypeError: If ``other`` is not a triangle (and ``verify=True``).
            NotImplementedError: If at least one of the triangles
                isn't two-dimensional (and ``verify=True``).
            ValueError: If ``strategy`` is not a valid
                :class:`.IntersectionStrategy`.
        """
        if verify:
            if not isinstance(other, Triangle):
                raise TypeError(
                    "Can only intersect with another triangle",
                    "Received",
                    other,
                )

            if self._dimension != 2 or other._dimension != 2:
                raise NotImplementedError(
                    "Intersection only implemented in 2D"
                )

        if strategy == _STRATEGY.GEOMETRIC:
            do_intersect = _triangle_intersection.geometric_intersect
        elif strategy == _STRATEGY.ALGEBRAIC:
            do_intersect = _py_triangle_intersection.algebraic_intersect
        else:
            raise ValueError("Unexpected strategy.", strategy)

        edge_infos, contained, all_edge_nodes = do_intersect(
            self._nodes, self._degree, other._nodes, other._degree, verify
        )
        if edge_infos is None:
            if contained:
                return [self]

            else:
                return [other]

        else:
            return [
                _make_intersection(edge_info, all_edge_nodes)
                for edge_info in edge_infos
            ]

    def elevate(self):
        r"""Return a degree-elevated version of the current triangle.

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

        .. image:: ../../images/triangle_elevate.png
           :align: center

        .. doctest:: triangle-elevate
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 1.5, 3.0, 0.75, 2.25, 0.0],
           ...     [0.0, 0.0, 0.0, 1.5 , 2.25, 3.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> elevated = triangle.elevate()
           >>> elevated
           <Triangle (degree=3, dimension=2)>
           >>> elevated.nodes
           array([[0. , 1. , 2. , 3. , 0.5 , 1.5 , 2.5 , 0.5 , 1.5 , 0. ],
                  [0. , 0. , 0. , 0. , 1.  , 1.25, 1.5 , 2.  , 2.5 , 3. ]])

        .. testcleanup:: triangle-elevate

           import make_images
           make_images.triangle_elevate(triangle, elevated)

        Returns:
            Triangle: The degree-elevated triangle.
        """
        _, num_nodes = self._nodes.shape
        # (d + 1)(d + 2)/2 --> (d + 2)(d + 3)/2
        num_new = num_nodes + self._degree + 2
        new_nodes = np.zeros((self._dimension, num_new), order="F")
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
        for k in range(self._degree + 1):
            for j in range(self._degree + 1 - k):
                i = self._degree - j - k
                new_nodes[:, parent_i1] += (i + 1) * self._nodes[:, index]
                new_nodes[:, parent_i2] += (j + 1) * self._nodes[:, index]
                new_nodes[:, parent_i3] += (k + 1) * self._nodes[:, index]
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
        return Triangle(new_nodes, self._degree + 1, copy=False, verify=False)

    def to_symbolic(self):
        """Convert to a SymPy matrix representing :math:`B(s, t)`.

        .. note::

           This method requires SymPy.

        .. doctest:: triangle-to-symbolic

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.5, 1.0, -0.5, 0.0, -1.0],
           ...     [0.0, 0.0, 1.0,  0.0, 0.0,  0.0],
           ...     [0.0, 0.0, 0.0,  0.0, 0.0,  1.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> triangle.to_symbolic()
           Matrix([
           [s - t],
           [ s**2],
           [ t**2]])

        Returns:
            :class:`sympy.Matrix <sympy.matrices.dense.MutableDenseMatrix>`:
            The triangle :math:`B(s, t)`.
        """
        _, _, b_polynomial = _symbolic.triangle_as_polynomial(
            self._nodes, self._degree
        )
        return b_polynomial

    def implicitize(self):
        r"""Implicitize the triangle.

        .. note::

           This method requires SymPy.

        .. doctest:: triangle-implicitize

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.5, 1.0, -0.5, 0.0, -1.0],
           ...     [0.0, 0.0, 1.0,  0.0, 0.0,  0.0],
           ...     [0.0, 0.0, 0.0,  0.0, 0.0,  1.0],
           ... ])
           >>> triangle = bezier.Triangle(nodes, degree=2)
           >>> triangle.implicitize()
           (x**4 - 2*x**2*y - 2*x**2*z + y**2 - 2*y*z + z**2)**2

        Returns:
            :class:`sympy.Expr <sympy.core.expr.Expr>`: The function that
            defines the triangle in :math:`\mathbf{R}^3` via
            :math:`f(x, y, z) = 0`.

        Raises:
            ValueError: If the triangle's dimension is not ``3``.
        """
        if self._dimension != 3:
            raise ValueError(
                "Only a spatial (3D) triangle can be implicitized",
                "Current dimension",
                self._dimension,
            )

        return _symbolic.implicitize_triangle(self._nodes, self._degree)


def _make_intersection(edge_info, all_edge_nodes):
    """Convert a description of edges into a curved polygon.

    .. note::

       This is a helper used only by :meth:`.Triangle.intersect`.

    Args:
        edge_info (Tuple[Tuple[int, float, float], ...]): Information
            describing each edge in the curved polygon by indicating which
            triangle / edge on the triangle and then start and end parameters
            along that edge. (See :func:`.ends_to_curve`.)
        all_edge_nodes (Tuple[numpy.ndarray, ...]): The nodes of three edges
            of the first triangle being intersected followed by the nodes of
            the three edges of the second.

    Returns:
        .CurvedPolygon: The intersection corresponding to ``edge_info``.
    """
    edges = []
    for index, start, end in edge_info:
        nodes = all_edge_nodes[index]
        new_nodes = _curve_helpers.specialize_curve(nodes, start, end)
        degree = new_nodes.shape[1] - 1
        edge = _curve_mod.Curve(new_nodes, degree, copy=False, verify=False)
        edges.append(edge)
    return curved_polygon.CurvedPolygon(
        *edges, metadata=edge_info, verify=False
    )
