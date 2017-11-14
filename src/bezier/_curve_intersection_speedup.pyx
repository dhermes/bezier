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

"""Cython "wrapped" interface for `curve_intersection.f90`."""


from libcpp cimport bool as bool_t
import numpy as np
from numpy cimport ndarray as ndarray_t

cimport bezier._curve_intersection


BoxIntersectionType_INTERSECTION = (
    bezier._curve_intersection.BoxIntersectionType.INTERSECTION)
BoxIntersectionType_TANGENT = (
    bezier._curve_intersection.BoxIntersectionType.TANGENT)
BoxIntersectionType_DISJOINT = (
    bezier._curve_intersection.BoxIntersectionType.DISJOINT)
cdef double[::1, :] WORKSPACE = np.empty((2, 2), order='F')
TOO_MANY_TEMPLATE = (
    'The number of candidate intersections is too high.\n'
    '{:d} candidate pairs.')
TOO_SMALL_TEMPLATE = (
    'Did not have enough space for intersections. Needed space '
    'for {:d} intersections but only had space for {:d}.')


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
        double[::1, :] root_nodes1,
        double error2, double start2, double end2,
        double[::1, :] start_node2, double[::1, :] end_node2,
        double[::1, :] root_nodes2):
    cdef int num_nodes1, num_nodes2
    cdef double refined_s, refined_t
    cdef bool_t does_intersect
    cdef int py_exc

    # NOTE: We don't check that there are 2 columns.
    num_nodes1, _ = np.shape(root_nodes1)
    # NOTE: We don't check that there are 2 columns.
    num_nodes2, _ = np.shape(root_nodes2)

    bezier._curve_intersection.from_linearized(
        &error1,
        &start1,
        &end1,
        &start_node1[0, 0],
        &end_node1[0, 0],
        &num_nodes1,
        &root_nodes1[0, 0],
        &error2,
        &start2,
        &end2,
        &start_node2[0, 0],
        &end_node2[0, 0],
        &num_nodes2,
        &root_nodes2[0, 0],
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


def reset_workspace(int workspace_size):
    global WORKSPACE
    WORKSPACE = np.empty((2, workspace_size), order='F')


def workspace_size():
    global WORKSPACE
    cdef int intersections_size

    # NOTE: We don't check that there are 2 rows.
    _, intersections_size = np.shape(WORKSPACE)
    return intersections_size


def all_intersections(
        double[::1, :] nodes_first, double[::1, :] nodes_second,
        bool_t allow_resize=True):
    global WORKSPACE
    cdef int num_nodes_first, num_nodes_second
    cdef int intersections_size, num_intersections, status
    cdef ndarray_t[double, ndim=2, mode='fortran'] intersections

    # NOTE: We don't check that there are 2 columns.
    num_nodes_first, _ = np.shape(nodes_first)
    num_nodes_second, _ = np.shape(nodes_second)
    # NOTE: We don't check that there are 2 rows.
    _, intersections_size = np.shape(WORKSPACE)

    bezier._curve_intersection.curve_intersections(
        &num_nodes_first,
        &nodes_first[0, 0],
        &num_nodes_second,
        &nodes_second[0, 0],
        &intersections_size,
        &WORKSPACE[0, 0],
        &num_intersections,
        &status,
    )

    if status == bezier._curve_intersection.AllIntersectionsStatus.SUCCESS:
        intersections = np.empty((num_intersections, 2), order='F')
        intersections[:, :] = WORKSPACE[:, :num_intersections].T
        return intersections
    elif (status ==
              bezier._curve_intersection.AllIntersectionsStatus.NO_CONVERGE):
        # NOTE: This assumes, but does not verify, that the Fortran subroutine
        #       uses ``MAX_INTERSECT_SUBDIVISIONS = 20``.
        raise ValueError(
            'Curve intersection failed to converge to approximately linear '
            'subdivisions after 20 iterations.')
    elif status == bezier._curve_intersection.AllIntersectionsStatus.TOO_SMALL:
        if allow_resize:
            reset_workspace(num_intersections)
            return all_intersections(
                nodes_first, nodes_second, allow_resize=False)
        else:
            msg = TOO_SMALL_TEMPLATE.format(
                intersections_size, num_intersections)
            raise ValueError(msg)
    elif status == bezier._curve_intersection.AllIntersectionsStatus.PARALLEL:
        raise NotImplementedError('Line segments parallel.')
    elif (status ==
              bezier._curve_intersection.AllIntersectionsStatus.WIGGLE_FAIL):
        # NOTE: This branch may not be tested because it's quite difficult to
        #       come up with an example that causes it.
        raise ValueError('outside of unit interval')
    elif status == bezier._curve_intersection.AllIntersectionsStatus.UNKNOWN:
        # NOTE: We exclude this block from testing because it **should**
        #       never occur. It's just a "future-proofing" mechanism of the
        #       Fortran code.
        raise ValueError('Unknown error.')
    else:
        # NOTE: If ``status`` isn't one of the enum values, then it is the
        #       number of candidate intersections.
        raise NotImplementedError(TOO_MANY_TEMPLATE.format(status))


def free_curve_intersections_workspace():
    bezier._curve_intersection.free_curve_intersections_workspace()
