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

"""Cython "wrapped" interface for `_intersection_helpers`."""


from libcpp cimport bool as bool_t
import numpy as np

cimport bezier._curve_intersection


BoxIntersectionType_INTERSECTION = (
    bezier._curve_intersection.BoxIntersectionType.INTERSECTION)
BoxIntersectionType_TANGENT = (
    bezier._curve_intersection.BoxIntersectionType.TANGENT)
BoxIntersectionType_DISJOINT = (
    bezier._curve_intersection.BoxIntersectionType.DISJOINT)


def linearization_error(double[::1, :] nodes):
    cdef int num_nodes, dimension
    cdef double error

    num_nodes, dimension = np.shape(nodes)

    bezier._curve_intersection.linearization_error(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &error,
    )
    return error


def segment_intersection(
        double[::1, :] start0, double[::1, :] end0,
        double[::1, :] start1, double[::1, :] end1):
    cdef double s, t
    cdef bool_t success

    bezier._curve_intersection.segment_intersection(
        &start0[0, 0],
        &end0[0, 0],
        &start1[0, 0],
        &end1[0, 0],
        &s,
        &t,
        &success,
    )
    if success:
        return s, t, True
    else:
        return None, None, False


def newton_refine(
        double s, double[::1, :] nodes1, double t, double[::1, :] nodes2):
    cdef int num_nodes1, num_nodes2
    cdef double new_s, new_t

    # NOTE: We don't check that there are 2 columns.
    num_nodes1, _ = np.shape(nodes1)
    # NOTE: We don't check that there are 2 columns.
    num_nodes2, _ = np.shape(nodes2)

    bezier._curve_intersection.newton_refine_intersect(
        &s,
        &num_nodes1,
        &nodes1[0, 0],
        &t,
        &num_nodes2,
        &nodes2[0, 0],
        &new_s,
        &new_t,
    )

    return new_s, new_t


def bbox_intersect(double[::1, :] nodes1, double[::1, :] nodes2):
    cdef int num_nodes1, num_nodes2, enum_val

    # NOTE: We don't check that there are 2 columns.
    num_nodes1, _ = np.shape(nodes1)
    # NOTE: We don't check that there are 2 columns.
    num_nodes2, _ = np.shape(nodes2)

    bezier._curve_intersection.bbox_intersect(
        &num_nodes1,
        &nodes1[0, 0],
        &num_nodes2,
        &nodes2[0, 0],
        &enum_val,
    )

    return enum_val


def parallel_different(
        double[::1, :] start0, double[::1, :] end0,
        double[::1, :] start1, double[::1, :] end1):
    cdef bool_t result

    bezier._curve_intersection.parallel_different(
        &start0[0, 0],
        &end0[0, 0],
        &start1[0, 0],
        &end1[0, 0],
        &result,
    )
    return result


def from_linearized_low_level(
        double error1, double start1, double end1,
        double[::1, :] start_node1, double[::1, :] end_node1,
        double[::1, :] nodes1,
        double error2, double start2, double end2,
        double[::1, :] start_node2, double[::1, :] end_node2,
        double[::1, :] nodes2):
    cdef int num_nodes1, num_nodes2
    cdef double refined_s, refined_t
    cdef bool_t does_intersect
    cdef int py_exc

    # NOTE: We don't check that there are 2 columns.
    num_nodes1, _ = np.shape(nodes1)
    # NOTE: We don't check that there are 2 columns.
    num_nodes2, _ = np.shape(nodes2)

    bezier._curve_intersection.from_linearized(
        &error1,
        &start1,
        &end1,
        &start_node1[0, 0],
        &end_node1[0, 0],
        &num_nodes1,
        &nodes1[0, 0],
        &error2,
        &start2,
        &end2,
        &start_node2[0, 0],
        &end_node2[0, 0],
        &num_nodes2,
        &nodes2[0, 0],
        &refined_s,
        &refined_t,
        &does_intersect,
        &py_exc,
    )

    if py_exc == 1:
        raise NotImplementedError('Line segments parallel.')
    elif py_exc == 2:
        raise ValueError('outside of unit interval')

    return refined_s, refined_t, does_intersect


def bbox_line_intersect(
        double[::1, :] nodes,
        double[::1, :] line_start, double[::1, :] line_end):
    cdef int num_nodes, enum_val

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)

    bezier._curve_intersection.bbox_line_intersect(
        &num_nodes,
        &nodes[0, 0],
        &line_start[0, 0],
        &line_end[0, 0],
        &enum_val,
    )

    return enum_val
