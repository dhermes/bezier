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


class Curve(object):
    r"""Represents a `Bezier curve`_.

    .. _Bezier curve: https://en.wikipedia.org/wiki/B%C3%A9zier_curve

    We take the traditional definiton: a `Bezier curve`_ is a mapping from
    :math:`s \in \left[0, 1\right]` to convex combinations
    of points :math:`v_0, v_1, \ldots, v_n` in some vector space:

    .. math::

       B(s) = \sum_{j = 0}^n \binom{n}{j} s^j (1 - s)^{n - j} \cdot v_j

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
