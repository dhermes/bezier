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
"""


import numpy as np
import six

from bezier import _base


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
      ...     [0.0, 0.0],
      ...     [0.625, 0.5],
      ...     [1.0, 0.5],
      ... ])
      >>> curve = bezier.Curve(nodes)
      >>> curve
      <Curve (degree=2, dimension=2)>

    Args:
        nodes (numpy.ndarray): The nodes in the curve. The rows
            represent each node while the columns are the dimension
            of the ambient space.
    """

    _length = None

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

    def evaluate(self, s):
        r"""Evaluate :math:`B(s)` along the curve.

        Performs `de Casteljau's algorithm`_ to build up :math:`B(s)`.

        .. _de Casteljau's algorithm:
            https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm

        .. doctest:: curve-eval
          :options: +NORMALIZE_WHITESPACE

          >>> nodes = np.array([
          ...     [0.0, 0.0],
          ...     [0.625, 0.5],
          ...     [1.0, 0.5],
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
        # NOTE: This assumes degree > 0, since the constructor requires this.
        weights = np.zeros((self.degree, self.degree + 1))
        eye = np.eye(self.degree)
        weights[:, 1:] += eye * s
        weights[:, :-1] += eye * (1 - s)

        value = weights.dot(self._nodes)
        for stage in six.moves.xrange(1, self.degree):
            value = weights[:-stage, :-stage].dot(value)

        # Here: Value will be 1x2, we just want the 1D point.
        return value.flatten()

    def evaluate_multi(self, s_vals):
        r"""Evaluate :math:`B(s)` for multiple points along the curve.

        Performs `de Casteljau's algorithm`_ to build up :math:`B(s)`.

        .. _de Casteljau's algorithm:
            https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm

        .. note::

            This currently just uses :meth:`evaluate` and so is less
            performant than it could be.

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
          array([[ 0.  ,  0.  ,  0.  ],
                 [ 0.25,  0.5 ,  0.75],
                 [ 0.5 ,  1.  ,  1.5 ],
                 [ 0.75,  1.5 ,  2.25],
                 [ 1.  ,  2.  ,  3.  ]])

        Args:
            s_vals (numpy.ndarray): Parameters along the curve (as a
                1D array).

        Returns:
            numpy.ndarray: The points on the curve. As a two dimensional
            NumPy array, with the rows corresponding to each ``s``
            value and the columns to the dimension.
        """
        num_vals, = s_vals.shape
        result = np.zeros((num_vals, self.dimension))
        for i, s_val in enumerate(s_vals):
            result[i, :] = self.evaluate(s_val)
        return result

    # pylint: disable=missing-param-doc,missing-type-doc,missing-returns-doc
    def plot(self, num_pts, plt, show=False):
        """Plot the current curve.

        Args:
            num_pts (int): Number of points to plot.
            plt (~types.ModuleType): Plotting module (i.e.
                :mod:`plt <matplotlib.pyplot>`) to use for creating
                figures, etc.
            show (Optional[bool]): Flag indicating if the plot should be
                shown.

        Returns:
            ~matplotlib.figure.Figure: The figure created for the plot.

        Raises:
            NotImplementedError: If the curve's dimension is not ``2``.
        """
        if self.dimension != 2:
            raise NotImplementedError('2D is the only supported dimension',
                                      'Current dimension', self.dimension)

        s_vals = np.linspace(0.0, 1.0, num_pts)
        points = self.evaluate_multi(s_vals)

        fig = plt.figure()
        ax = fig.gca()
        ax.plot(points[:, 0], points[:, 1])

        if show:
            plt.show()

        return fig
    # pylint: enable=missing-param-doc,missing-type-doc,missing-returns-doc

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
          ...     [0.0, 0.0],
          ...     [1.25, 3.0],
          ...     [2.0, 1.0],
          ... ])
          >>> curve = bezier.Curve(nodes)
          >>> left, right = curve.subdivide()
          >>> left
          <Curve (degree=2, dimension=2)>
          >>> left.nodes
          array([[ 0.   , 0.   ],
                 [ 0.625, 1.5  ],
                 [ 1.125, 1.75 ]])
          >>> right
          <Curve (degree=2, dimension=2)>
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

        left = new_nodes[:self.degree + 1, :]
        right = new_nodes[self.degree:, :]

        return Curve(left), Curve(right)
