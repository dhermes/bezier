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

.. autofunction:: bezier._intersection_helpers.linearization_error
.. autofunction:: bezier._intersection_helpers.newton_refine
.. autofunction:: bezier._intersection_helpers.segment_intersection
"""


import matplotlib.pyplot as plt
import numpy as np

from bezier import _base
from bezier import _curve_helpers
from bezier import _intersection_helpers


_REPR_TEMPLATE = (
    '<{} (degree={:d}, dimension={:d}, start={:g}, end={:g})>')
_LINEAR_SUBDIVIDE = np.array([
    [1.0, 0.0],
    [0.5, 0.5],
    [0.0, 1.0],
])
_QUADRATIC_SUBDIVIDE = np.array([
    [1.0, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.25, 0.5, 0.25],
    [0.0, 0.5, 0.5],
    [0.0, 0.0, 1.0],
])
_CUBIC_SUBDIVIDE = np.array([
    [1.0, 0.0, 0.0, 0.0],
    [0.5, 0.5, 0.0, 0.0],
    [0.25, 0.5, 0.25, 0.0],
    [0.125, 0.375, 0.375, 0.125],
    [0.0, 0.25, 0.5, 0.25],
    [0.0, 0.0, 0.5, 0.5],
    [0.0, 0.0, 0.0, 1.0],
])


class Curve(_base.Base):
    r"""Represents a B |eacute| zier `curve`_.

    .. _curve: https://en.wikipedia.org/wiki/B%C3%A9zier_curve

    We take the traditional definition: a B |eacute| zier curve is a mapping
    from :math:`s \in \left[0, 1\right]` to convex combinations
    of points :math:`v_0, v_1, \ldots, v_n` in some vector space:

    .. math::

       B(s) = \sum_{j = 0}^n \binom{n}{j} s^j (1 - s)^{n - j} \cdot v_j

    .. doctest:: curve-ctor

       >>> import bezier
       >>> nodes = np.array([
       ...     [0.0  , 0.0],
       ...     [0.625, 0.5],
       ...     [1.0  , 0.5],
       ... ])
       >>> curve = bezier.Curve(nodes)
       >>> curve
       <Curve (degree=2, dimension=2)>

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
    """

    _length = None

    def __init__(self, nodes, start=0.0, end=1.0, root=None, _copy=True):
        super(Curve, self).__init__(nodes, _copy=_copy)
        self._start = start
        self._end = end
        if root is None:
            root = self
        self._root = root

    def __repr__(self):
        """Representation of current object.

        Returns:
            str: Object representation.
        """
        if self.start == 0.0 and self.end == 1.0:
            return super(Curve, self).__repr__()
        else:
            return _REPR_TEMPLATE.format(
                self.__class__.__name__, self.degree,
                self.dimension, self.start, self.end)

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
                self._nodes, self.degree)
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

           >>> nodes = np.array([
           ...     [0.0, 0.0],
           ...     [1.0, 2.0],
           ... ])
           >>> curve = bezier.Curve(nodes)
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

           nodes = np.array([
               [0.0, 0.0],
               [0.75, 0.0],
               [1.0, 1.0],
           ])
           curve = bezier.Curve(nodes)

        .. doctest:: curve-root
           :options: +NORMALIZE_WHITESPACE

           >>> _, right = curve.subdivide()
           >>> right
           <Curve (degree=2, dimension=2, start=0.5, end=1)>
           >>> right.root is curve
           True
           >>> right.evaluate(0.0) == curve.evaluate(0.5)
           array([ True, True], dtype=bool)
           >>>
           >>> mid_left, _ = right.subdivide()
           >>> mid_left
           <Curve (degree=2, dimension=2, start=0.5, end=0.75)>
           >>> mid_left.root is curve
           True
           >>> mid_left.evaluate(1.0) == curve.evaluate(0.75)
           array([ True, True], dtype=bool)
        """
        return self._root

    def evaluate(self, s):
        r"""Evaluate :math:`B(s)` along the curve.

        See :meth:`evaluate_multi` for more details.

        .. doctest:: curve-eval
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.array([
           ...     [0.0  , 0.0],
           ...     [0.625, 0.5],
           ...     [1.0  , 0.5],
           ... ])
           >>> curve = bezier.Curve(nodes)
           >>> curve.evaluate(0.75)
           array([ 0.796875, 0.46875 ])

        Args:
            s (float): Parameter along the curve.

        Returns:
            numpy.ndarray: The point on the curve (as a one dimensional
            NumPy array).
        """
        result = self.evaluate_multi(np.array([s]))
        return result.flatten()

    def evaluate_multi(self, s_vals):
        r"""Evaluate :math:`B(s)` for multiple points along the curve.

        This is done by first evaluating each member of the
        `Bernstein basis`_ at each value in ``s_vals`` and then
        applying those to the control points for the current curve.

        This is done instead of using `de Casteljau's algorithm`_.
        Implementing de Casteljau is problematic because it requires
        a choice between one of two methods:

        * vectorize operations of the form :math:`(1 - s)v + s w`,
          which requires a copy of the curve's control points for
          each value in ``s_vals``
        * avoid vectorization and compute each point in serial

        Instead, we can use vectorized operations to build up the
        Bernstein basis values.

        .. _Bernstein basis:
            https://en.wikipedia.org/wiki/Bernstein_polynomial
        .. _de Casteljau's algorithm:
            https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm

        .. doctest:: curve-eval-multi
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.array([
           ...     [0.0, 0.0, 0.0],
           ...     [1.0, 2.0, 3.0],
           ... ])
           >>> curve = bezier.Curve(nodes)
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
        return _curve_helpers.evaluate_multi(
            self._nodes, self.degree, s_vals)

    def plot(self, num_pts, ax=None, show=False):
        """Plot the current curve.

        Args:
            num_pts (int): Number of points to plot.
            ax (Optional[matplotlib.artist.Artist]): matplotlib axis object
                to add plot to.
            show (Optional[bool]): Flag indicating if the plot should be
                shown.

        Returns:
            matplotlib.artist.Artist: The axis containing the plot. This
            may be a newly created axis.

        Raises:
            NotImplementedError: If the curve's dimension is not ``2``.
        """
        if self.dimension != 2:
            raise NotImplementedError('2D is the only supported dimension',
                                      'Current dimension', self.dimension)

        s_vals = np.linspace(0.0, 1.0, num_pts)
        points = self.evaluate_multi(s_vals)

        if ax is None:
            fig = plt.figure()
            ax = fig.gca()

        ax.plot(points[:, 0], points[:, 1])

        if show:
            plt.show()

        return ax

    def subdivide(self):
        r"""Split the curve :math:`\gamma(s)` into a left and right half.

        Takes the interval :math:`\left[0, 1\right]` and splits the curve into
        :math:`\gamma_1 = \gamma\left(\left[0, \frac{1}{2}\right]\right)` and
        :math:`\gamma_2 = \gamma\left(\left[\frac{1}{2}, 1\right]\right)`. In
        order to do this, also reparameterizes the curve, hence the resulting
        left and right halves have new nodes.

        .. doctest:: curve-subdivide
           :options: +NORMALIZE_WHITESPACE

           >>> nodes = np.array([
           ...     [0.0 , 0.0],
           ...     [1.25, 3.0],
           ...     [2.0 , 1.0],
           ... ])
           >>> curve = bezier.Curve(nodes)
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

        Returns:
            Tuple[Curve, Curve]: The left and right sub-curves.
        """
        if self.degree == 1:
            # pylint: disable=no-member
            new_nodes = _LINEAR_SUBDIVIDE.dot(self._nodes)
            # pylint: enable=no-member
        elif self.degree == 2:
            # pylint: disable=no-member
            new_nodes = _QUADRATIC_SUBDIVIDE.dot(self._nodes)
            # pylint: enable=no-member
        elif self.degree == 3:
            # pylint: disable=no-member
            new_nodes = _CUBIC_SUBDIVIDE.dot(self._nodes)
            # pylint: enable=no-member
        else:
            subdivide_mat = _curve_helpers.make_subdivision_matrix(
                self.degree)
            new_nodes = subdivide_mat.dot(self._nodes)

        left_nodes = new_nodes[:self.degree + 1, :]
        right_nodes = new_nodes[self.degree:, :]

        midpoint = 0.5 * (self.start + self.end)
        left = Curve(left_nodes, start=self.start,
                     end=midpoint, root=self.root, _copy=False)
        right = Curve(right_nodes, start=midpoint,
                      end=self.end, root=self.root, _copy=False)
        return left, right

    def intersect(self, other):
        """Find the points of intersection with another curve.

        .. doctest:: curve-intersect
           :options: +NORMALIZE_WHITESPACE

           >>> nodes1 = np.array([
           ...     [0.0  , 0.0 ],
           ...     [0.375, 0.75 ],
           ...     [0.75 , 0.375],
           ... ])
           >>> curve1 = bezier.Curve(nodes1)
           >>> nodes2 = np.array([
           ...     [0.5, 0.0],
           ...     [0.5, 3.0],
           ... ])
           >>> curve2 = bezier.Curve(nodes2)
           >>> curve1.intersect(curve2)
           array([[ 0.5, 0.5]])

        Args:
            other (Curve): Other curve to intersect with.

        Returns:
            numpy.ndarray: Array of intersection points (possibly empty).

        Raises:
            TypeError: If ``other`` is not a curve.
            NotImplementedError: If both curves aren't two-dimensional.
        """
        if not isinstance(other, Curve):
            raise TypeError('Can only intersect with another curve',
                            'Received', other)
        if self.dimension != 2 or other.dimension != 2:
            raise NotImplementedError(
                'Intersection only implemented in 2D')

        candidates = [(self, other)]
        intersections = _intersection_helpers.all_intersections(candidates)
        if intersections:
            return np.vstack([intersection.point
                              for intersection in intersections])
        else:
            return np.zeros((0, 2))

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

        Returns:
            Curve: The degree-elevated curve.
        """
        new_nodes = _curve_helpers.elevate_nodes(
            self._nodes, self.degree, self.dimension)
        return Curve(new_nodes, _copy=False)
