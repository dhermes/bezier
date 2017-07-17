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

"""Helper for B |eacute| zier Curves.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

See :doc:`../curve-curve-intersection` for examples using the
:class:`Curve` class to find intersections.

.. testsetup:: *

   import numpy as np
   import bezier

.. autoclass:: IntersectionStrategy
   :members:
"""


import numpy as np

from bezier import _base
from bezier import _curve_helpers
from bezier import _helpers
from bezier import _intersection_helpers
from bezier import _plot_helpers


_REPR_TEMPLATE = (
    '<{} (degree={:d}, dimension={:d}, start={:g}, end={:g})>')
_LOCATE_ERROR_TEMPLATE = (
    'Dimension mismatch: This curve is {:d}-dimensional, so the point should '
    'be a 1x{:d} NumPy array. Instead the point {} has dimensions {}.')
_LINEAR_SUBDIVIDE_LEFT = np.asfortranarray([
    [1.0, 0.0],
    [0.5, 0.5],
])
_LINEAR_SUBDIVIDE_RIGHT = np.asfortranarray([
    [0.5, 0.5],
    [0.0, 1.0],
])
_QUADRATIC_SUBDIVIDE_LEFT = np.asfortranarray([
    [1.0, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.25, 0.5, 0.25],
])
_QUADRATIC_SUBDIVIDE_RIGHT = np.asfortranarray([
    [0.25, 0.5, 0.25],
    [0.0, 0.5, 0.5],
    [0.0, 0.0, 1.0],
])
_CUBIC_SUBDIVIDE_LEFT = np.asfortranarray([
    [1.0, 0.0, 0.0, 0.0],
    [0.5, 0.5, 0.0, 0.0],
    [0.25, 0.5, 0.25, 0.0],
    [0.125, 0.375, 0.375, 0.125],
])
_CUBIC_SUBDIVIDE_RIGHT = np.asfortranarray([
    [0.125, 0.375, 0.375, 0.125],
    [0.0, 0.25, 0.5, 0.25],
    [0.0, 0.0, 0.5, 0.5],
    [0.0, 0.0, 0.0, 1.0],
])


IntersectionStrategy = _intersection_helpers.IntersectionStrategy


class Curve(_base.Base):
    r"""Represents a B |eacute| zier `curve`_.

    .. _curve: https://en.wikipedia.org/wiki/B%C3%A9zier_curve

    We take the traditional definition: a B |eacute| zier curve is a mapping
    from :math:`s \in \left[0, 1\right]` to convex combinations
    of points :math:`v_0, v_1, \ldots, v_n` in some vector space:

    .. math::

       B(s) = \sum_{j = 0}^n \binom{n}{j} s^j (1 - s)^{n - j} \cdot v_j

    .. image:: ../images/curve_constructor.png
       :align: center

    .. doctest:: curve-constructor

       >>> import bezier
       >>> nodes = np.asfortranarray([
       ...     [0.0  , 0.0],
       ...     [0.625, 0.5],
       ...     [1.0  , 0.5],
       ... ])
       >>> curve = bezier.Curve.from_nodes(nodes)
       >>> curve
       <Curve (degree=2, dimension=2)>

    .. testcleanup:: curve-constructor

       import make_images
       make_images.curve_constructor(curve)

    Args:
        nodes (numpy.ndarray): The nodes in the curve. The rows
            represent each node while the columns are the dimension
            of the ambient space.
        degree (int): The degree of the curve. This is assumed to
            correctly correspond to the number of ``nodes``. Use
            :meth:`from_nodes` if the degree has not yet been computed.
        start (Optional[float]): The beginning of the sub-interval
            that this curve represents.
        end (Optional[float]): The end of the sub-interval
            that this curve represents.
        root (Optional[Curve]): The root curve that contains this
            current curve.
        _copy (bool): Flag indicating if the nodes should be copied before
            being stored. Defaults to :data:`True` since callers may
            freely mutate ``nodes`` after passing in.
    """

    __slots__ = (
        '_dimension', '_nodes',  # From base class
        '_degree', '_start', '_end', '_root',  # From constructor
        '_length', '_edge_index', '_next_edge',  # Empty defaults
        '_previous_edge',  # Empty defaults
    )

    def __init__(self, nodes, degree, start=0.0, end=1.0,
                 root=None, _copy=True):
        super(Curve, self).__init__(nodes, _copy=_copy)
        self._degree = degree
        self._start = start
        self._end = end
        if root is None:
            root = self
        self._root = root
        self._length = None
        self._edge_index = None
        self._next_edge = None
        self._previous_edge = None

    @classmethod
    def from_nodes(cls, nodes, start=0.0, end=1.0, root=None, _copy=True):
        """Create a :class:`.Curve` from nodes.

        Computes the ``degree`` based on the shape of ``nodes``.

        Args:
            nodes (numpy.ndarray): The nodes in the curve. The rows
                represent each node while the columns are the dimension
                of the ambient space.
            start (Optional[float]): The beginning of the sub-interval
                that this curve represents.
            end (Optional[float]): The end of the sub-interval
                that this curve represents.
            root (Optional[Curve]): The root curve that contains this
                current curve.
            _copy (bool): Flag indicating if the nodes should be copied before
                being stored. Defaults to :data:`True` since callers may
                freely mutate ``nodes`` after passing in.

        Returns:
            Curve: The constructed curve.
        """
        num_nodes, _ = nodes.shape
        degree = cls._get_degree(num_nodes)
        return cls(nodes, degree, start=start, end=end, root=root, _copy=_copy)

    def __repr__(self):
        """Representation of current object.

        Returns:
            str: Object representation.
        """
        if self._start == 0.0 and self._end == 1.0:
            return super(Curve, self).__repr__()
        else:
            return _REPR_TEMPLATE.format(
                self.__class__.__name__, self._degree,
                self._dimension, self._start, self._end)

    @staticmethod
    def _get_degree(num_nodes):
        """Get the degree of the current curve.

        Args:
            num_nodes (int): The number of nodes provided.

        Returns:
            int: The degree of the current curve.
        """
        return num_nodes - 1

    @property
    def length(self):
        """float: The length of the current curve."""
        if self._length is None:
            self._length = _curve_helpers.compute_length(
                self._nodes, self._degree)
        return self._length

    @property
    def start(self):
        """float: Start of sub-interval this curve represents.

        This value is used to track the current curve in the
        re-parameterization / subdivision process. The curve is still
        defined on the unit interval, but this value illustrates
        how this curve relates to a "parent" curve. For example:

        .. doctest:: curve-start
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.0],
           ...     [1.0, 2.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=1)
           >>> curve
           <Curve (degree=1, dimension=2)>
           >>> left, right = curve.subdivide()
           >>> left
           <Curve (degree=1, dimension=2, start=0, end=0.5)>
           >>> right
           <Curve (degree=1, dimension=2, start=0.5, end=1)>
           >>> _, mid_right = left.subdivide()
           >>> mid_right
           <Curve (degree=1, dimension=2, start=0.25, end=0.5)>
           >>> mid_right.nodes
           array([[ 0.25, 0.5 ],
                  [ 0.5 , 1.  ]])
        """
        return self._start

    @property
    def end(self):
        """float: End of sub-interval this curve represents.

        See :attr:`~Curve.start` for more information.
        """
        return self._end

    @property
    def root(self):
        """Curve: The "root" curve that contains the current curve.

        This indicates that the current curve is a section of the
        "root" curve. For example:

        .. testsetup:: curve-root

           import numpy as np
           import bezier

           nodes = np.asfortranarray([
               [0.0, 0.0],
               [0.75, 0.0],
               [1.0, 1.0],
           ])
           curve = bezier.Curve(nodes, degree=2)

        .. doctest:: curve-root
           :options: +NORMALIZE_WHITESPACE

           >>> _, right = curve.subdivide()
           >>> right
           <Curve (degree=2, dimension=2, start=0.5, end=1)>
           >>> right.root is curve
           True
           >>> right.evaluate(0.0) == curve.evaluate(0.5)
           array([[ True, True]], dtype=bool)
           >>>
           >>> mid_left, _ = right.subdivide()
           >>> mid_left
           <Curve (degree=2, dimension=2, start=0.5, end=0.75)>
           >>> mid_left.root is curve
           True
           >>> mid_left.evaluate(1.0) == curve.evaluate(0.75)
           array([[ True, True]], dtype=bool)
        """
        return self._root

    @property
    def edge_index(self):
        """Optional[int]: The index of the edge among a group of edges.

        .. testsetup:: edge-index

           import numpy as np
           import bezier

           nodes = np.asfortranarray([
               [0.0, 0.0],
               [1.0, 0.0],
               [0.0, 1.0],
           ])
           surface = bezier.Surface(nodes, degree=1)
           _, curve, _ = surface.edges

        .. doctest:: edge-index

           >>> curve.edge_index
           1
           >>> curve.previous_edge
           <Curve (degree=1, dimension=2)>
           >>> curve.previous_edge.edge_index
           0
           >>> curve.next_edge
           <Curve (degree=1, dimension=2)>
           >>> curve.next_edge.edge_index
           2

        This is intended to be used when a :class:`Curve` is created as part
        of a larger structure like a :class:`.Surface` or
        :class:`.CurvedPolygon`.
        """
        return self._edge_index

    @property
    def next_edge(self):
        """Optional[Curve]: An edge that comes after the current one.

        This is intended to be used when a :class:`Curve` is created as part
        of a larger structure like a :class:`.Surface` or
        :class:`.CurvedPolygon`.
        """
        return self._next_edge

    @property
    def previous_edge(self):
        """Optional[Curve]: An edge that comes before the current one.

        This is intended to be used when a :class:`Curve` is created as part
        of a larger structure like a :class:`.Surface` or
        :class:`.CurvedPolygon`.
        """
        return self._previous_edge

    @property
    def __dict__(self):
        """dict: Dictionary of current curve's property namespace.

        This is just a stand-in property for the usual ``__dict__``. This
        class defines ``__slots__`` so by default would not provide a
        ``__dict__``.

        This also means that the current object can't be modified by the
        returned dictionary.
        """
        return {name: getattr(self, name)
                for name in self.__slots__}

    def _copy(self):
        """Make a copy of the current curve.

        Returns:
            .Curve: Copy of current curve.
        """
        root = None if self._root is self else self._root
        result = Curve(
            self._nodes, self._degree, start=self._start, end=self._end,
            root=root, _copy=True)
        # Also copy over any cached computed values. (Ignore the
        # prev/next/index information though: YAGNI.)
        result._length = self._length  # pylint: disable=protected-access
        return result

    def evaluate(self, s):
        r"""Evaluate :math:`B(s)` along the curve.

        This method acts as a (partial) inverse to :meth:`locate`.

        See :meth:`evaluate_multi` for more details.

        .. image:: ../images/curve_evaluate.png
           :align: center

        .. doctest:: curve-eval
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0  , 0.0],
           ...     [0.625, 0.5],
           ...     [1.0  , 0.5],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=2)
           >>> curve.evaluate(0.75)
           array([[ 0.796875, 0.46875 ]])

        .. testcleanup:: curve-eval

           import make_images
           make_images.curve_evaluate(curve)

        Args:
            s (float): Parameter along the curve.

        Returns:
            numpy.ndarray: The point on the curve (as a two dimensional
            NumPy array with a single row).
        """
        return _curve_helpers.evaluate_multi(
            self._nodes, np.asfortranarray([s]))

    def evaluate_multi(self, s_vals):
        r"""Evaluate :math:`B(s)` for multiple points along the curve.

        This is done via a modified Horner's method (vectorized for
        each ``s``-value).

        .. doctest:: curve-eval-multi
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.0, 0.0],
           ...     [1.0, 2.0, 3.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=1)
           >>> curve
           <Curve (degree=1, dimension=3)>
           >>> s_vals = np.linspace(0.0, 1.0, 5)
           >>> curve.evaluate_multi(s_vals)
           array([[ 0.  , 0.  , 0.  ],
                  [ 0.25, 0.5 , 0.75],
                  [ 0.5 , 1.  , 1.5 ],
                  [ 0.75, 1.5 , 2.25],
                  [ 1.  , 2.  , 3.  ]])

        Args:
            s_vals (numpy.ndarray): Parameters along the curve (as a
                1D array).

        Returns:
            numpy.ndarray: The points on the curve. As a two dimensional
            NumPy array, with the rows corresponding to each ``s``
            value and the columns to the dimension.
        """
        return _curve_helpers.evaluate_multi(self._nodes, s_vals)

    def plot(self, num_pts, color=None, alpha=None, ax=None):
        """Plot the current curve.

        Args:
            num_pts (int): Number of points to plot.
            color (Optional[Tuple[float, float, float]]): Color as RGB profile.
            alpha (Optional[float]): The alpha channel for the color.
            ax (Optional[matplotlib.artist.Artist]): matplotlib axis object
                to add plot to.

        Returns:
            matplotlib.artist.Artist: The axis containing the plot. This
            may be a newly created axis.

        Raises:
            NotImplementedError: If the curve's dimension is not ``2``.
        """
        if self._dimension != 2:
            raise NotImplementedError('2D is the only supported dimension',
                                      'Current dimension', self._dimension)

        s_vals = np.linspace(0.0, 1.0, num_pts)
        points = self.evaluate_multi(s_vals)

        if ax is None:
            ax = _plot_helpers.new_axis()

        ax.plot(points[:, 0], points[:, 1], color=color, alpha=alpha)

        return ax

    def subdivide(self):
        r"""Split the curve :math:`B(s)` into a left and right half.

        Takes the interval :math:`\left[0, 1\right]` and splits the curve into
        :math:`B_1 = B\left(\left[0, \frac{1}{2}\right]\right)` and
        :math:`B_2 = B\left(\left[\frac{1}{2}, 1\right]\right)`. In
        order to do this, also reparameterizes the curve, hence the resulting
        left and right halves have new nodes.

        .. image:: ../images/curve_subdivide.png
           :align: center

        .. doctest:: curve-subdivide
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0 , 0.0],
           ...     [1.25, 3.0],
           ...     [2.0 , 1.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=2)
           >>> left, right = curve.subdivide()
           >>> left
           <Curve (degree=2, dimension=2, start=0, end=0.5)>
           >>> left.nodes
           array([[ 0.   , 0.   ],
                  [ 0.625, 1.5  ],
                  [ 1.125, 1.75 ]])
           >>> right
           <Curve (degree=2, dimension=2, start=0.5, end=1)>
           >>> right.nodes
           array([[ 1.125, 1.75 ],
                  [ 1.625, 2.   ],
                  [ 2.   , 1.   ]])

        .. testcleanup:: curve-subdivide

           import make_images
           make_images.curve_subdivide(curve, left, right)

        Returns:
            Tuple[Curve, Curve]: The left and right sub-curves.
        """
        if self._degree == 1:
            left_nodes = _helpers.matrix_product(
                _LINEAR_SUBDIVIDE_LEFT, self._nodes)
            right_nodes = _helpers.matrix_product(
                _LINEAR_SUBDIVIDE_RIGHT, self._nodes)
        elif self._degree == 2:
            left_nodes = _helpers.matrix_product(
                _QUADRATIC_SUBDIVIDE_LEFT, self._nodes)
            right_nodes = _helpers.matrix_product(
                _QUADRATIC_SUBDIVIDE_RIGHT, self._nodes)
        elif self._degree == 3:
            left_nodes = _helpers.matrix_product(
                _CUBIC_SUBDIVIDE_LEFT, self._nodes)
            right_nodes = _helpers.matrix_product(
                _CUBIC_SUBDIVIDE_RIGHT, self._nodes)
        else:
            left_mat, right_mat = _curve_helpers.make_subdivision_matrices(
                self._degree)
            left_nodes = _helpers.matrix_product(left_mat, self._nodes)
            right_nodes = _helpers.matrix_product(right_mat, self._nodes)

        midpoint = 0.5 * (self._start + self._end)
        left = Curve(
            left_nodes, self._degree, start=self._start, end=midpoint,
            root=self._root, _copy=False)
        right = Curve(
            right_nodes, self._degree, start=midpoint, end=self._end,
            root=self._root, _copy=False)
        return left, right

    def intersect(self, other,
                  strategy=IntersectionStrategy.geometric,
                  _verify=True):
        """Find the points of intersection with another curve.

        See :doc:`../curve-curve-intersection` for more details.

        .. image:: ../images/curve_intersect.png
           :align: center

        .. doctest:: curve-intersect
           :options: +NORMALIZE_WHITESPACE

           >>> nodes1 = np.asfortranarray([
           ...     [0.0  , 0.0  ],
           ...     [0.375, 0.75 ],
           ...     [0.75 , 0.375],
           ... ])
           >>> curve1 = bezier.Curve(nodes1, degree=2)
           >>> nodes2 = np.asfortranarray([
           ...     [0.5, 0.0 ],
           ...     [0.5, 0.75],
           ... ])
           >>> curve2 = bezier.Curve(nodes2, degree=1)
           >>> intersections = curve1.intersect(curve2)
           >>> intersections
           array([[ 0.5, 0.5]])

        .. testcleanup:: curve-intersect

           import make_images
           make_images.curve_intersect(curve1, curve2, intersections)

        Args:
            other (Curve): Other curve to intersect with.
            strategy (Optional[~bezier.curve.IntersectionStrategy]): The
                intersection algorithm to use. Defaults to geometric.
            _verify (Optional[bool]): Indicates if extra caution should be
                used to verify assumptions about the input and current
                curve. Can be disabled to speed up execution time.
                Defaults to :data:`True`.

        Returns:
            numpy.ndarray: Array of intersection points (possibly empty).

        Raises:
            TypeError: If ``other`` is not a curve (and ``_verify=True``).
            NotImplementedError: If at least one of the curves
                isn't two-dimensional (and ``_verify=True``).
        """
        if _verify:
            if not isinstance(other, Curve):
                raise TypeError('Can only intersect with another curve',
                                'Received', other)
            if self._dimension != 2 or other._dimension != 2:
                raise NotImplementedError(
                    'Intersection only implemented in 2D')

        candidates = [(self, other)]
        intersections = _intersection_helpers.all_intersections(
            candidates, strategy=strategy)
        if intersections:
            result = np.empty((len(intersections), 2), order='F')
            for index, intersection in enumerate(intersections):
                result[index, :] = intersection.get_point()
            return result
        else:
            return np.zeros((0, 2), order='F')

    def elevate(self):
        r"""Return a degree-elevated version of the current curve.

        Does this by converting the current nodes :math:`v_0, \ldots, v_n`
        to new nodes :math:`w_0, \ldots, w_{n + 1}` where

        .. math::

           \begin{align*}
           w_0 &= v_0 \\
           w_j &= \frac{j}{n + 1} v_{j - 1} + \frac{n + 1 - j}{n + 1} v_j \\
           w_{n + 1} &= v_n
           \end{align*}

        .. image:: ../images/curve_elevate.png
           :align: center

        .. testsetup:: curve-elevate

           import numpy as np
           import bezier

        .. doctest:: curve-elevate
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.0],
           ...     [1.5, 1.5],
           ...     [3.0, 0.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=2)
           >>> elevated = curve.elevate()
           >>> elevated
           <Curve (degree=3, dimension=2)>
           >>> elevated.nodes
           array([[ 0., 0.],
                  [ 1., 1.],
                  [ 2., 1.],
                  [ 3., 0.]])

        .. testcleanup:: curve-elevate

           import make_images
           make_images.curve_elevate(curve, elevated)

        Returns:
            Curve: The degree-elevated curve.
        """
        new_nodes = _curve_helpers.elevate_nodes(
            self._nodes, self._degree, self._dimension)
        root = None if self._root is self else self._root
        return Curve(
            new_nodes, self._degree + 1, start=self._start, end=self._end,
            root=root, _copy=False)

    def reduce_(self):
        r"""Return a degree-reduced version of the current curve.

        .. _pseudo-inverse:
            https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_pseudoinverse

        Does this by converting the current nodes :math:`v_0, \ldots, v_n`
        to new nodes :math:`w_0, \ldots, w_{n - 1}` that correspond to
        reversing the :meth:`elevate` process.

        This uses the `pseudo-inverse`_ of the elevation matrix. For example
        when elevating from degree 2 to 3, the matrix :math:`E_2` is given by

        .. math::

           \mathbf{v} = \left[\begin{array}{c} v_0 \\ v_1 \\ v_2
               \end{array}\right] \longmapsto \left[\begin{array}{c} v_0 \\
               \frac{v_0 + 2 v_1}{3} \\ \frac{2 v_1 + v_2}{3} \\ v_2
               \end{array}\right] = \frac{1}{3}
               \left[\begin{array}{c c} 3 & 0 \\ 2 & 1 \\ 1 & 2
               \\ 0 & 3 \end{array}\right] \mathbf{v}

        and the pseudo-inverse is given by

        .. math::

           R_2 = \left(E_2^T E_2\right)^{-1} E_2^T = \frac{1}{20}
               \left[\begin{array}{c c c c} 19 & 3 & -3 & 1 \\
               -5 & 15 & 15 & -5 \\ 1 & -3 & 3 & 19
               \end{array}\right].

        .. warning::

           Though degree-elevation preserves the start and end nodes, degree
           reduction has no such guarantee. Rather, the nodes produced are
           "best" in the least squares sense (when solving the normal
           equations).

        .. image:: ../images/curve_reduce.png
           :align: center

        .. testsetup:: curve-reduce, curve-reduce-approx

           import numpy as np
           import bezier

        .. doctest:: curve-reduce
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [-3.0, 3.0],
           ...     [ 0.0, 2.0],
           ...     [ 1.0, 3.0],
           ...     [ 0.0, 6.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=3)
           >>> reduced = curve.reduce_()
           >>> reduced
           <Curve (degree=2, dimension=2)>
           >>> reduced.nodes
           array([[-3. , 3. ],
                  [ 1.5, 1.5],
                  [ 0. , 6. ]])

        .. testcleanup:: curve-reduce

           import make_images
           make_images.curve_reduce(curve, reduced)

        In the case that the current curve **is not** degree-elevated.

        .. image:: ../images/curve_reduce_approx.png
           :align: center

        .. doctest:: curve-reduce-approx
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0 , 2.5],
           ...     [1.25, 5.0],
           ...     [3.75, 7.5],
           ...     [5.0 , 2.5],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=3)
           >>> reduced = curve.reduce_()
           >>> reduced
           <Curve (degree=2, dimension=2)>
           >>> reduced.nodes
           array([[-0.125, 2.125],
                  [ 2.5  , 8.125],
                  [ 5.125, 2.875]])

        .. testcleanup:: curve-reduce-approx

           import make_images
           make_images.curve_reduce_approx(curve, reduced)

        Returns:
            Curve: The degree-reduced curve.
        """
        new_nodes = _curve_helpers.reduce_pseudo_inverse(
            self._nodes, self._degree)
        root = None if self._root is self else self._root
        return Curve(
            new_nodes, self._degree - 1, start=self._start, end=self._end,
            root=root, _copy=False)

    def specialize(self, start, end):
        """Specialize the curve to a given sub-interval.

        .. image:: ../images/curve_specialize.png
           :align: center

        .. doctest:: curve-specialize

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.0],
           ...     [0.5, 1.0],
           ...     [1.0, 0.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=2)
           >>> new_curve = curve.specialize(-0.25, 0.75)
           >>> new_curve
           <Curve (degree=2, dimension=2, start=-0.25, end=0.75)>
           >>> new_curve.nodes
           array([[-0.25 , -0.625],
                  [ 0.25 ,  0.875],
                  [ 0.75 ,  0.375]])

        .. testcleanup:: curve-specialize

           import make_images
           make_images.curve_specialize(curve, new_curve)

        This is generalized version of :meth:`subdivide`, and can even
        match the output of that method:

        .. testsetup:: curve-specialize2

           import numpy as np
           import bezier

           nodes = np.asfortranarray([
               [0.0, 0.0],
               [0.5, 1.0],
               [1.0, 0.0],
           ])
           curve = bezier.Curve(nodes, degree=2)

        .. doctest:: curve-specialize2

           >>> left, right = curve.subdivide()
           >>> also_left = curve.specialize(0.0, 0.5)
           >>> np.all(also_left.nodes == left.nodes)
           True
           >>> also_right = curve.specialize(0.5, 1.0)
           >>> np.all(also_right.nodes == right.nodes)
           True

        Args:
            start (float): The start point of the interval we
                are specializing to.
            end (float): The end point of the interval we
                are specializing to.

        Returns:
            Curve: The newly-specialized curve.
        """
        new_nodes, true_start, true_end = _curve_helpers.specialize_curve(
            self._nodes, start, end, self._start, self._end, self._degree)
        return Curve(
            new_nodes, self._degree, start=true_start, end=true_end,
            root=self._root, _copy=False)

    def locate(self, point):
        r"""Find a point on the current curve.

        Solves for :math:`s` in :math:`B(s) = p`.

        This method acts as a (partial) inverse to :meth:`evaluate`.

        .. note::

           A unique solution is only guaranteed if the current curve has no
           self-intersections. This code assumes, but doesn't check, that
           this is true.

        .. image:: ../images/curve_locate.png
           :align: center

        .. doctest:: curve-locate

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.0],
           ...     [1.0, 2.0],
           ...     [3.0, 1.0],
           ...     [4.0, 0.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=3)
           >>> point1 = np.asfortranarray([[3.09375, 0.703125]])
           >>> s = curve.locate(point1)
           >>> s
           0.75
           >>> point2 = np.asfortranarray([[2.0, 0.5]])
           >>> curve.locate(point2) is None
           True

        .. testcleanup:: curve-locate

           import make_images
           make_images.curve_locate(curve, point1, point2)

        Args:
            point (numpy.ndarray): A (``1xD``) point on the curve,
                where :math:`D` is the dimension of the curve.

        Returns:
            Optional[float]: The parameter value (:math:`s`) corresponding
            to ``point`` or :data:`None` if the point is not on
            the ``curve``.

        Raises:
            ValueError: If the dimension of the ``point`` doesn't match the
                dimension of the current curve.
        """
        if point.shape != (1, self._dimension):
            point_dimensions = ' x '.join(
                str(dimension) for dimension in point.shape)
            msg = _LOCATE_ERROR_TEMPLATE.format(
                self._dimension, self._dimension, point, point_dimensions)
            raise ValueError(msg)

        return _curve_helpers.locate_point(self, point)
