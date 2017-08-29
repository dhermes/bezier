#!python
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

"""Cython "wrapped" interface for `_helpers`."""


from libcpp cimport bool as bool_t
import numpy as np

cimport bezier._helpers


def cross_product(double[::1, :] vec0, double[::1, :] vec1):
    cdef double result

    bezier._helpers.cross_product(
        &vec0[0, 0],
        &vec1[0, 0],
        &result,
    )
    return result


def bbox(double[::1, :] nodes):
    cdef int num_nodes
    cdef double left, right, bottom, top

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)

    bezier._helpers.bbox(
        &num_nodes,
        &nodes[0, 0],
        &left,
        &right,
        &bottom,
        &top,
    )

    return left, right, bottom, top


def wiggle_interval(double value):
    cdef double result
    cdef bool_t success

    bezier._helpers.wiggle_interval(
        &value,
        &result,
        &success,
    )
    return result, success
