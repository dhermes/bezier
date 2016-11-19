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

"""Private helper methods for B |eacute| zier curves.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


import numpy as np
import six


def evaluate_multi(nodes, degree, s_vals):
    r"""Computes multiple points along a curve.

    Does so by computing the Bernstein basis at each value in ``s_vals``
    rather than using the de Casteljau algorithm.

    Args:
        nodes (numpy.ndarray): The nodes defining a curve.
        degree (int): The degree of the curve (assumed to be one less than
            the number of ``nodes``.
        s_vals (numpy.ndarray): Parameters along the curve (as a
            1D array).

    Returns:
        numpy.ndarray: The evaluated points on the curve as a two dimensional
        NumPy array, with the rows corresponding to each ``s``
        value and the columns to the dimension.
    """
    num_vals, = s_vals.shape

    lambda2 = s_vals[:, np.newaxis]
    lambda1 = 1.0 - lambda2

    weights_next = np.zeros((num_vals, degree + 1))
    weights_curr = np.zeros((num_vals, degree + 1))
    weights_curr[:, 0] = 1.0

    # Increase from degree 0 to ``degree``.
    for curr_deg in six.moves.xrange(degree):
        weights_next[:, :curr_deg + 1] = (
            lambda1 * weights_curr[:, :curr_deg + 1])
        weights_next[:, 1:curr_deg + 2] += (
            lambda2 * weights_curr[:, :curr_deg + 1])
        weights_curr, weights_next = weights_next, weights_curr

    return weights_curr.dot(nodes)
