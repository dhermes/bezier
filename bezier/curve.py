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

.. testsetup:: *

  import numpy as np
  import bezier

.. autofunction:: bezier._intersection_helpers.linearization_error
.. autofunction:: bezier._intersection_helpers.newton_refine
.. autofunction:: bezier._intersection_helpers.segment_intersection
"""


import itertools

import matplotlib.pyplot as plt
import numpy as np
import six

from bezier import _base
from bezier import _curve_helpers
from bezier import _intersection_helpers


_REPR_TEMPLATE = (
    '<{} (degree={:d}, dimension={:d}, start={:g}, end={:g})>')
_FREXP = np.frexp  # pylint: disable=no-member
_MAX_INTERSECT_SUBDIVISIONS = 20
_ERROR_EXPONENT = -26
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


def _make_subdivision_matrix(degree):
    """Make the matrix used to subdivide a curve.

    Args:
        degree (int): The degree of the curve.

    Returns:
        numpy.ndarray: The matrix used to convert the
           nodes into left and right nodes.
    """
    num_rows = 2 * degree + 1
    result = np.zeros((num_rows, degree + 1))
    result[0, 0] = 1.0
    result[-1, -1] = 1.0
    for row in six.moves.xrange(1, degree + 1):
        half_prev = 0.5 * result[row - 1, :row]
        result[row, :row] = half_prev
        result[row, 1:row + 1] += half_prev
        # Populate the complement row as well.
        complement = num_rows - row - 1
        # NOTE: We "should" reverse the results when using
        #       the complement, but they are symmetric so
        #       that would be a waste.
        result[complement, -(row + 1):] = result[row, :row + 1]
    return result


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
    def length(self):  # pylint: disable=missing-returns-doc
        """float: The length of the current curve.

        Raises:
            NotImplementedError: If the length isn't already cached.
        """
        if self._length is None:
            raise NotImplementedError(
                'Length computation not yet implemented.')
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
            subdivide_mat = _make_subdivision_matrix(self.degree)
            new_nodes = subdivide_mat.dot(self._nodes)

        left_nodes = new_nodes[:self.degree + 1, :]
        right_nodes = new_nodes[self.degree:, :]

        root = self.root
        if root is None:
            root = self

        midpoint = 0.5 * (self.start + self.end)
        left = Curve(left_nodes, start=self.start,
                     end=midpoint, root=root, _copy=False)
        right = Curve(right_nodes, start=midpoint,
                      end=self.end, root=root, _copy=False)
        return left, right

    @staticmethod
    def _from_linearized(linearized_pairs):
        """Determine curve-curve intersections from pairs of linearizations.

        Args:
            linearized_pairs (list): List of pairs of
                :class:`._intersection_helpers.Linearization` instances.

        Returns:
            numpy.ndarray: Array of all intersections.
        """
        intersections = []
        for left, right in linearized_pairs:
            s, t = _intersection_helpers.segment_intersection(
                left.start, left.end, right.start, right.end)
            left_curve = left._curve  # pylint: disable=protected-access
            right_curve = right._curve  # pylint: disable=protected-access
            # TODO: Check if s, t are in [0, 1].
            # Now, promote `s` and `t` onto the original curves.
            orig_s = (1 - s) * left_curve.start + s * left_curve.end
            orig_left = left_curve.root
            orig_t = (1 - t) * right_curve.start + t * right_curve.end
            orig_right = right_curve.root
            # Perform one step of Newton iteration to refine the computed
            # values of s and t.
            refined_s, _ = _intersection_helpers.newton_refine(
                orig_s, orig_left, orig_t, orig_right)
            # TODO: Check that (orig_left.evaluate(refined_s) ~=
            #                   orig_right.evaluate(refined_t))
            intersections.append(orig_left.evaluate(refined_s))

        return np.vstack(intersections)

    def intersect(self, other):
        """Find the points of intersection with another curve.

        Args:
            other (Curve): Other curve to intersect with.

        Returns:
            numpy.ndarray: Possible empty array of intersection points.

        Raises:
            TypeError: If ``other`` is not a curve.
            NotImplementedError: If both curves aren't two-dimensional.
            ValueError: If the subdivision iteration does not terminate
                before exhausting the maximum number of subdivisions.
        """
        if not isinstance(other, Curve):
            raise TypeError('Can only intersect with another curve',
                            'Received', other)
        if self.dimension != 2 or other.dimension != 2:
            raise NotImplementedError(
                'Intersection only implemented in 2D')

        candidates = [(self, other)]
        accepted = None
        for _ in six.moves.xrange(_MAX_INTERSECT_SUBDIVISIONS):
            accepted = []
            max_err = 0.0

            for left, right in candidates:
                # pylint: disable=protected-access
                left_nodes = left._nodes
                right_nodes = right._nodes
                # pylint: enable=protected-access
                if not _intersection_helpers.bbox_intersect(
                        left_nodes, right_nodes):
                    continue

                # Attempt to replace the curves with linearizations
                # if they are close enough to lines.
                if isinstance(left, _intersection_helpers.Linearization):
                    is_curve = False
                    # pylint: disable=no-member,protected-access
                    err_left = left._error
                    # pylint: enable=no-member,protected-access
                else:
                    is_curve = True
                    # NOTE: This may be a wasted computation, e.g. if ``left``
                    #       occurs in multiple pairs. However, in practice the
                    #       number of pairs will be small so this cost
                    #       will be low.
                    err_left = _intersection_helpers.linearization_error(left)

                max_err = max(max_err, err_left)
                if is_curve:
                    _, left_exp = _FREXP(err_left)
                    if left_exp <= _ERROR_EXPONENT:
                        left = _intersection_helpers.Linearization(
                            left, error=err_left)

                # Now do the same for the right.
                if isinstance(right, _intersection_helpers.Linearization):
                    is_curve = False
                    # pylint: disable=no-member,protected-access
                    err_right = right._error
                    # pylint: enable=no-member,protected-access
                else:
                    is_curve = True
                    # NOTE: This may also be a wasted computation.
                    err_right = _intersection_helpers.linearization_error(
                        right)

                max_err = max(max_err, err_right)
                if is_curve:
                    _, right_exp = _FREXP(err_right)
                    if right_exp <= _ERROR_EXPONENT:
                        right = _intersection_helpers.Linearization(
                            right, error=err_right)

                # Add the accepted pair.
                accepted.append((left, right))

            # If none of the pairs have been accepted, then there is
            # no intersection.
            if not accepted:
                return np.zeros((0, 2))
            # In the case of ``accepted`` pairs, if the pairs are
            # sufficiently close to their linearizations, we can stop
            # the subdivisions and move on to the next step.
            _, max_exp = _FREXP(max_err)
            if max_exp <= _ERROR_EXPONENT:
                return self._from_linearized(accepted)
            # If we **do** require more subdivisions, we need to update
            # the list of candidates.
            # pylint: disable=redefined-variable-type
            candidates = itertools.chain(*[
                itertools.product(left.subdivide(), right.subdivide())
                for left, right in accepted])
            # pylint: enable=redefined-variable-type

        return ValueError(
            'Curve intersection failed to converge to approximately '
            'linear subdivisions after max iterations.',
            _MAX_INTERSECT_SUBDIVISIONS)
