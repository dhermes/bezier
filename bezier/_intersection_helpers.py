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

"""Private helper methods for intersecting B |eacute| zier shapes.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


import numpy as np


def bbox_intersect(nodes1, nodes2):
    r"""Bounding box intersection predicate.

    Determines if the bounding box of two sets of control points
    intersects in :math:`\mathbf{R}^2` with non-trivial
    intersection (i.e. tangent bounding boxes are insufficient).

    .. note::

        Though we assume (and the code relies on this fact) that
        the nodes are two-dimensional, we don't check it.

    Args:
        nodes1 (numpy.ndarray): Set of control points for a
            B |eacute| zier shape.
        nodes2 (numpy.ndarray): Set of control points for a
            B |eacute| zier shape.

    Returns:
        bool: Predicate indicating if the bounding boxes intersect.
    """
    left1, bottom1 = np.min(nodes1, axis=0)
    right1, top1 = np.max(nodes1, axis=0)
    left2, bottom2 = np.min(nodes2, axis=0)
    right2, top2 = np.max(nodes2, axis=0)

    return (right2 > left1 and right1 > left2 and
            top2 > bottom1 and top1 > bottom2)
