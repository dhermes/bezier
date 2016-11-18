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


def de_casteljau_multi(nodes, degree, s_vals):
    r"""Performs the de Casteljau algorithm for multiple points along a curve.

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

    # Make a broadcasted copy along an extra axis. We insert
    # the axis in between the #nodes and the dimension since
    # the number of nodes in the result should come before
    # the dimension.
    value = np.repeat(nodes[:, np.newaxis, :], num_vals, axis=1)
    # Put the parameter values on this axis for broadcasting.
    s_vals = s_vals[np.newaxis, :, np.newaxis]
    t_vals = 1.0 - s_vals

    for _ in six.moves.xrange(degree):
        value = t_vals * value[:-1, :, :] + s_vals * value[1:, :, :]

    # Here: Value will be 1x2x(num_vals), we just want the 2D points.
    return value[0, :, :]
