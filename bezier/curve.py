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

"""Helper for Bezier Curves."""


import numpy as np
import six


class Curve(object):
    r"""Represents a `Bezier curve`_.

    .. _Bezier curve: https://en.wikipedia.org/wiki/B%C3%A9zier_curve

    We take the traditional definiton: a `Bezier curve`_ is a mapping from
    :math:`s \in \left[0, 1\right]` to convex combinations
    of points :math:`v_0, v_1, \ldots, v_n` in some vector space:

    .. math::

       B(s) = \sum_{j = 0}^n \binom{n}{j} s^j (1 - s)^{n - j} \cdot v_j

    .. testsetup:: curve-eval

      import numpy as np

    .. doctest:: curve-eval
      :options: +NORMALIZE_WHITESPACE

      >>> import bezier
      >>> nodes = np.array([
      ...     [0.0, 0.0],
      ...     [0.625, 0.5],
      ...     [1.0, 0.5],
      ... ])
      >>> curve = bezier.Curve(nodes)
      >>> curve.evaluate(0.75)
      array([ 0.796875, 0.46875 ])

    Args:
        nodes (numpy.ndarray): The nodes in the curve. The rows
            represent each node while the columns are the dimension
            of the ambient space.
    """

    def __init__(self, nodes):
        rows, cols = nodes.shape
        self._degree = rows - 1
        self._dimension = cols
        self._nodes = nodes

    @property
    def degree(self):
        """int: The degree of the current curve."""
        return self._degree

    @property
    def dimension(self):
        r"""int: The dimension that the curve lives in.

        For example, if the curve is 3D, i.e. if
        :math:`B(s) \in \mathbf{R}^3`, then the dimension is ``3``.
        """
        return self._dimension

    def evaluate(self, s):
        r"""Evaluate :math:`B(s)` along the curve.

        Performs `de Casteljau's algorithm`_ to build up :math:`B(s)`.

        .. _de Casteljau's algorithm:
            https://en.wikipedia.org/wiki/De_Casteljau%27s_algorithm

        Args:
            s (float): Parameter along the curve.

        Returns:
            numpy.ndarray: The point on the curve (as a one dimensional
                NumPy array).
        """
        if self.degree == 0:
            return self._nodes.flatten()

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

            This current just uses :meth:`evaluate` and so is less
            performant than it could be.

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

    def plot(self, num_pts, plt, show=False):
        """Plot the current curve.

        Args:
            num_pts (int): Number of points to plot.
            plt (~types.ModuleType): Plotting module to use for creating
               figures, etc.
            show (bool): (Optional) Flag indicating if the plot should be
                shown.

        Returns:
            matplotlib.figure.Figure: The figure created for the plot.
        """
        s_vals = np.linspace(0.0, 1.0, num_pts)
        points = self.evaluate_multi(s_vals)

        fig = plt.figure()
        ax = fig.gca()
        ax.plot(points[:, 0], points[:, 1])

        if show:
            plt.show()

        return fig
