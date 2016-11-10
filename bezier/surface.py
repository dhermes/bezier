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
"""


class Surface(object):
    r"""Represents a B |eacute| zier `surface`_.

    .. _surface: https://en.wikipedia.org/wiki/B%C3%A9zier_triangle

    Args:
        nodes (numpy.ndarray): The nodes in the surface. The rows
            represent each node while the columns are the dimension
            of the ambient space.

    Raises:
        ValueError: If the ``nodes`` are not 2D.
    """

    def __init__(self, nodes):
        if nodes.ndim != 2:
            raise ValueError('Nodes must be 2-dimensional, not', nodes.ndim)
        _, cols = nodes.shape
        self._dimension = cols
        self._nodes = nodes
