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


cdef double EPS = 0.5**40


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


def contains_nd(double[::1, :] nodes, double[::1, :] point):
    cdef int num_nodes, dimension
    cdef bool_t predicate

    num_nodes, dimension = np.shape(nodes)
    if np.shape(point) != (1, dimension):
        msg = 'Point {} was expected to have shape (1, {})'.format(
            np.asarray(point), dimension)
        raise ValueError(msg)

    bezier._helpers.contains_nd(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &point[0, 0],
        &predicate,
    )

    return predicate


def vector_close(double[::1, :] vec1, double[::1, :] vec2, double eps=EPS):
    cdef int num_values

    # NOTE: We don't check that there is 1 row or that
    #       ``np.shape(vec1) == np.shape(vec2)``.
    _, num_values = np.shape(vec1)

    return bezier._helpers.vector_close(
        &num_values,
        &vec1[0, 0],
        &vec2[0, 0],
        &eps,
    )
