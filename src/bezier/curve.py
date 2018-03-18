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

from bezier import _algebraic_intersection
from bezier import _base
from bezier import _curve_helpers
from bezier import _geometric_intersection
from bezier import _intersection_helpers
from bezier import _plot_helpers

_LOCATE_ERROR_TEMPLATE = (
    'Dimension mismatch: This curve is {:d}-dimensional, so the point should '
    'be a {:d} x 1 NumPy array. Instead the point {} has dimensions {}.'
)
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
       ...     [0.0, 0.625, 1.0],
       ...     [0.0, 0.5  , 0.5],
       ... ])
       >>> curve = bezier.Curve(nodes, degree=2)
       >>> curve
       <Curve (degree=2, dimension=2)>

    .. testcleanup:: curve-constructor

       import make_images
       make_images.curve_constructor(curve)

    Args:
        nodes (numpy.ndarray): The nodes in the curve. The columns
            represent each node while the rows are the dimension
            of the ambient space.
        degree (int): The degree of the curve. This is assumed to
            correctly correspond to the number of ``nodes``. Use
            :meth:`from_nodes` if the degree has not yet been computed.
        _copy (bool): Flag indicating if the nodes should be copied before
            being stored. Defaults to :data:`True` since callers may
            freely mutate ``nodes`` after passing in.
    """
    __slots__ = (
        '_dimension',  # From base class
        '_nodes',  # From base class
        '_degree',  # From constructor
        '_length',  # Empty defaults
    )

    def __init__(self, nodes, degree, _copy=True):
        super(Curve, self).__init__(nodes, _copy=_copy)
        self._degree = degree
        self._length = None

    @classmethod
    def from_nodes(cls, nodes, _copy=True):
        """Create a :class:`.Curve` from nodes.

        Computes the ``degree`` based on the shape of ``nodes``.

        Args:
            nodes (numpy.ndarray): The nodes in the curve. The columns
                represent each node while the rows are the dimension
                of the ambient space.
            _copy (bool): Flag indicating if the nodes should be copied before
                being stored. Defaults to :data:`True` since callers may
                freely mutate ``nodes`` after passing in.

        Returns:
            Curve: The constructed curve.
        """
        _, num_nodes = nodes.shape
        degree = cls._get_degree(num_nodes)
        return cls(nodes, degree, _copy=_copy)

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
            self._length = _curve_helpers.compute_length(self._nodes)
        return self._length

    @property
    def __dict__(self):
        """dict: Dictionary of current curve's property namespace.

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
            '_length': self._length,
        }

    def _copy(self):
        """Make a copy of the current curve.

        Returns:
            .Curve: Copy of current curve.
        """
        result = Curve(self._nodes, self._degree, _copy=True)
        # Also copy over any cached computed values.
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
           ...     [0.0, 0.625, 1.0],
           ...     [0.0, 0.5  , 0.5],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=2)
           >>> curve.evaluate(0.75)
           array([[0.796875],
                  [0.46875 ]])

        .. testcleanup:: curve-eval

           import make_images
           make_images.curve_evaluate(curve)

        Args:
            s (float): Parameter along the curve.

        Returns:
            numpy.ndarray: The point on the curve (as a two dimensional
            NumPy array with a single column).
        """
        return _curve_helpers.evaluate_multi(
            self._nodes, np.asfortranarray([s])
        )

    def evaluate_multi(self, s_vals):
        r"""Evaluate :math:`B(s)` for multiple points along the curve.

        This is done via a modified Horner's method (vectorized for
        each ``s``-value).

        .. doctest:: curve-eval-multi
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 1.0],
           ...     [0.0, 2.0],
           ...     [0.0, 3.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=1)
           >>> curve
           <Curve (degree=1, dimension=3)>
           >>> s_vals = np.linspace(0.0, 1.0, 5)
           >>> curve.evaluate_multi(s_vals)
           array([[0.  , 0.25, 0.5 , 0.75, 1.  ],
                  [0.  , 0.5 , 1.  , 1.5 , 2.  ],
                  [0.  , 0.75, 1.5 , 2.25, 3.  ]])

        Args:
            s_vals (numpy.ndarray): Parameters along the curve (as a
                1D array).

        Returns:
            numpy.ndarray: The points on the curve. As a two dimensional
            NumPy array, with the columns corresponding to each ``s``
            value and the rows to the dimension.
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
            raise NotImplementedError(
                '2D is the only supported dimension',
                'Current dimension',
                self._dimension,
            )

        s_vals = np.linspace(0.0, 1.0, num_pts)
        points = self.evaluate_multi(s_vals)
        if ax is None:
            ax = _plot_helpers.new_axis()
        ax.plot(points[0, :], points[1, :], color=color, alpha=alpha)
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
           ...     [0.0, 1.25, 2.0],
           ...     [0.0, 3.0 , 1.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=2)
           >>> left, right = curve.subdivide()
           >>> left.nodes
           array([[0.   , 0.625, 1.125],
                  [0.   , 1.5  , 1.75 ]])
           >>> right.nodes
           array([[1.125, 1.625, 2.   ],
                  [1.75 , 2.   , 1.   ]])

        .. testcleanup:: curve-subdivide

           import make_images
           make_images.curve_subdivide(curve, left, right)

        Returns:
            Tuple[Curve, Curve]: The left and right sub-curves.
        """
        left_nodes, right_nodes = _curve_helpers.subdivide_nodes(self._nodes)
        left = Curve(left_nodes, self._degree, _copy=False)
        right = Curve(right_nodes, self._degree, _copy=False)
        return left, right

    def intersect(
        self, other, strategy=IntersectionStrategy.GEOMETRIC, _verify=True
    ):
        """Find the points of intersection with another curve.

        See :doc:`../curve-curve-intersection` for more details.

        .. image:: ../images/curve_intersect.png
           :align: center

        .. doctest:: curve-intersect
           :options: +NORMALIZE_WHITESPACE

           >>> nodes1 = np.asfortranarray([
           ...     [0.0, 0.375, 0.75 ],
           ...     [0.0, 0.75 , 0.375],
           ... ])
           >>> curve1 = bezier.Curve(nodes1, degree=2)
           >>> nodes2 = np.asfortranarray([
           ...     [0.5, 0.5 ],
           ...     [0.0, 0.75],
           ... ])
           >>> curve2 = bezier.Curve(nodes2, degree=1)
           >>> intersections = curve1.intersect(curve2)
           >>> 3.0 * intersections
           array([[2.],
                  [2.]])
           >>> s_vals = intersections[0, :]
           >>> curve1.evaluate_multi(s_vals)
           array([[0.5],
                  [0.5]])

        .. testcleanup:: curve-intersect

           import make_images
           make_images.curve_intersect(curve1, curve2, s_vals)

        Args:
            other (Curve): Other curve to intersect with.
            strategy (Optional[~bezier.curve.IntersectionStrategy]): The
                intersection algorithm to use. Defaults to geometric.
            _verify (Optional[bool]): Indicates if extra caution should be
                used to verify assumptions about the input and current
                curve. Can be disabled to speed up execution time.
                Defaults to :data:`True`.

        Returns:
            numpy.ndarray: ``2 x N`` array of ``s``- and ``t``-parameters where
            intersections occur (possibly empty).

        Raises:
            TypeError: If ``other`` is not a curve (and ``_verify=True``).
            NotImplementedError: If at least one of the curves
                isn't two-dimensional (and ``_verify=True``).
            ValueError: If ``strategy`` is not a valid
                :attr:`.IntersectionStrategy`.
        """
        if _verify:
            if not isinstance(other, Curve):
                raise TypeError(
                    'Can only intersect with another curve', 'Received', other
                )

            if self._dimension != 2 or other._dimension != 2:
                raise NotImplementedError(
                    'Intersection only implemented in 2D'
                )

        if strategy == IntersectionStrategy.GEOMETRIC:
            all_intersections = _geometric_intersection.all_intersections
        elif strategy == IntersectionStrategy.ALGEBRAIC:
            all_intersections = _algebraic_intersection.all_intersections
        else:
            raise ValueError('Unexpected strategy.', strategy)

        st_vals, _ = all_intersections(self._nodes, other._nodes)
        return st_vals

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
           ...     [0.0, 1.5, 3.0],
           ...     [0.0, 1.5, 0.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=2)
           >>> elevated = curve.elevate()
           >>> elevated
           <Curve (degree=3, dimension=2)>
           >>> elevated.nodes
           array([[0., 1., 2., 3.],
                  [0., 1., 1., 0.]])

        .. testcleanup:: curve-elevate

           import make_images
           make_images.curve_elevate(curve, elevated)

        Returns:
            Curve: The degree-elevated curve.
        """
        new_nodes = _curve_helpers.elevate_nodes(self._nodes)
        return Curve(new_nodes, self._degree + 1, _copy=False)

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

           \mathbf{v} = \left[\begin{array}{c c c} v_0 & v_1 & v_2
               \end{array}\right] \longmapsto \left[\begin{array}{c c c c}
               v_0 & \frac{v_0 + 2 v_1}{3} & \frac{2 v_1 + v_2}{3} & v_2
               \end{array}\right] = \frac{1}{3} \mathbf{v}
               \left[\begin{array}{c c c c} 3 & 1 & 0 & 0 \\
               0 & 2 & 2 & 0 \\ 0 & 0 & 1 & 3 \end{array}\right]

        and the (right) pseudo-inverse is given by

        .. math::

           R_2 = E_2^T \left(E_2 E_2^T\right)^{-1} = \frac{1}{20}
               \left[\begin{array}{c c c} 19 & -5 & 1 \\
               3 & 15 & -3 \\ -3 & 15 & 3 \\ 1 & -5 & 19
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
           ...     [-3.0, 0.0, 1.0, 0.0],
           ...     [ 3.0, 2.0, 3.0, 6.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=3)
           >>> reduced = curve.reduce_()
           >>> reduced
           <Curve (degree=2, dimension=2)>
           >>> reduced.nodes
           array([[-3. ,  1.5,  0. ],
                  [ 3. ,  1.5,  6. ]])

        .. testcleanup:: curve-reduce

           import make_images
           make_images.curve_reduce(curve, reduced)

        In the case that the current curve **is not** degree-elevated.

        .. image:: ../images/curve_reduce_approx.png
           :align: center

        .. doctest:: curve-reduce-approx
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.asfortranarray([
           ...     [0.0, 1.25, 3.75, 5.0],
           ...     [2.5, 5.0 , 7.5 , 2.5],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=3)
           >>> reduced = curve.reduce_()
           >>> reduced
           <Curve (degree=2, dimension=2)>
           >>> reduced.nodes
           array([[-0.125,  2.5  ,  5.125],
                  [ 2.125,  8.125,  2.875]])

        .. testcleanup:: curve-reduce-approx

           import make_images
           make_images.curve_reduce_approx(curve, reduced)

        Returns:
            Curve: The degree-reduced curve.
        """
        new_nodes = _curve_helpers.reduce_pseudo_inverse(self._nodes)
        return Curve(new_nodes, self._degree - 1, _copy=False)

    def specialize(self, start, end):
        """Specialize the curve to a given sub-interval.

        .. image:: ../images/curve_specialize.png
           :align: center

        .. doctest:: curve-specialize

           >>> nodes = np.asfortranarray([
           ...     [0.0, 0.5, 1.0],
           ...     [0.0, 1.0, 0.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=2)
           >>> new_curve = curve.specialize(-0.25, 0.75)
           >>> new_curve.nodes
           array([[-0.25 ,  0.25 ,  0.75 ],
                  [-0.625,  0.875,  0.375]])

        .. testcleanup:: curve-specialize

           import make_images
           make_images.curve_specialize(curve, new_curve)

        This is generalized version of :meth:`subdivide`, and can even
        match the output of that method:

        .. testsetup:: curve-specialize2

           import numpy as np
           import bezier

           nodes = np.asfortranarray([
               [0.0, 0.5, 1.0],
               [0.0, 1.0, 0.0],
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
        new_nodes = _curve_helpers.specialize_curve(self._nodes, start, end)
        return Curve(new_nodes, self._degree, _copy=False)

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
           ...     [0.0, 1.0, 3.0, 4.0],
           ...     [0.0, 2.0, 1.0, 0.0],
           ... ])
           >>> curve = bezier.Curve(nodes, degree=3)
           >>> point1 = np.asfortranarray([
           ...     [3.09375 ],
           ...     [0.703125],
           ... ])
           >>> s = curve.locate(point1)
           >>> s
           0.75
           >>> point2 = np.asfortranarray([
           ...     [2.0],
           ...     [0.5],
           ... ])
           >>> curve.locate(point2) is None
           True

        .. testcleanup:: curve-locate

           import make_images
           make_images.curve_locate(curve, point1, point2)

        Args:
            point (numpy.ndarray): A (``D x 1``) point on the curve,
                where :math:`D` is the dimension of the curve.

        Returns:
            Optional[float]: The parameter value (:math:`s`) corresponding
            to ``point`` or :data:`None` if the point is not on
            the ``curve``.

        Raises:
            ValueError: If the dimension of the ``point`` doesn't match the
                dimension of the current curve.
        """
        if point.shape != (self._dimension, 1):
            point_dimensions = ' x '.join(
                str(dimension) for dimension in point.shape
            )
            msg = _LOCATE_ERROR_TEMPLATE.format(
                self._dimension, self._dimension, point, point_dimensions
            )
            raise ValueError(msg)

        return _curve_helpers.locate_point(self._nodes, point)
