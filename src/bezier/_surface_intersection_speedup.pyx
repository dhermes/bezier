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

"""Cython "wrapped" interface for `_surface_intersection`."""


import numpy as np
from numpy cimport dtype as dtype_t
from numpy cimport ndarray as ndarray_t

cimport bezier._status
cimport bezier._surface
cimport bezier._surface_intersection
from bezier._surface_intersection cimport CurvedPolygonSegment
from bezier._surface_intersection cimport SurfaceContained


cdef int[:] SEGMENT_ENDS_WORKSPACE = np.empty(3, dtype=np.intc)
cdef dtype_t SEGMENT_DTYPE = np.dtype(
    [
        ('start', np.double),
        ('end', np.double),
        ('edge_index', np.intc),
    ],
    align=True,
)
cdef CurvedPolygonSegment[:] SEGMENTS_WORKSPACE = np.empty(
    6, dtype=SEGMENT_DTYPE)
SEGMENT_ENDS_TOO_SMALL = (
    'Did not have enough space for segment ends. Needed space '
    'for {:d} integers but only had space for {:d}.')
SEGMENTS_TOO_SMALL = (
    'Did not have enough space for segments. Needed space '
    'for {:d} `CurvedPolygonSegment`-s but only had space for {:d}.')
# NOTE: ``TOO_MANY_TEMPLATE`` is copy-pasted from
#       ``_curve_intersection_speedup.pyx`` to avoid defining a `.pxd`
#       for either of these modules.
cdef str TOO_MANY_TEMPLATE = (
    'The number of candidate intersections is too high.\n'
    '{:d} candidate pairs.')


def newton_refine(
        double[::1, :] nodes, int degree,
        double x_val, double y_val, double s, double t):
    cdef int num_nodes
    cdef double updated_s, updated_t

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)

    bezier._surface_intersection.newton_refine_surface(
        &num_nodes,
        &nodes[0, 0],
        &degree,
        &x_val,
        &y_val,
        &s,
        &t,
        &updated_s,
        &updated_t,
    )

    return updated_s, updated_t


def locate_point(double[::1, :] nodes, int degree, double x_val, double y_val):
    cdef int num_nodes
    cdef double s_val, t_val

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)

    bezier._surface_intersection.locate_point_surface(
        &num_nodes,
        &nodes[0, 0],
        &degree,
        &x_val,
        &y_val,
        &s_val,
        &t_val,
    )

    if s_val == -1.0:  # LOCATE_MISS
        return None
    else:
        return s_val, t_val


def reset_workspaces(int segment_ends_size=-1, int segments_size=-1):
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    if segment_ends_size != -1:
        SEGMENT_ENDS_WORKSPACE = np.empty(segment_ends_size, dtype=np.intc)
    if segments_size != -1:
        SEGMENTS_WORKSPACE = np.empty(segments_size, dtype=SEGMENT_DTYPE)


def workspace_sizes():
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    cdef int segment_ends_size
    cdef int segments_size

    segment_ends_size, = np.shape(SEGMENT_ENDS_WORKSPACE)
    segments_size, = np.shape(SEGMENTS_WORKSPACE)
    return segment_ends_size, segments_size


def _surface_intersections_success(
        double[::1, :] nodes1, int degree1,
        double[::1, :] nodes2, int degree2,
        int num_intersected):
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    cdef int begin_index, end_index
    cdef size_t i, j
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode='fortran'] edge_nodes1
    cdef ndarray_t[double, ndim=2, mode='fortran'] edge_nodes2
    cdef ndarray_t[double, ndim=2, mode='fortran'] edge_nodes3
    cdef ndarray_t[double, ndim=2, mode='fortran'] edge_nodes4
    cdef ndarray_t[double, ndim=2, mode='fortran'] edge_nodes5
    cdef ndarray_t[double, ndim=2, mode='fortran'] edge_nodes6

    curved_polygons = []
    for i in range(num_intersected):
        if i == 0:
            begin_index = 0
        else:
            begin_index = SEGMENT_ENDS_WORKSPACE[i - 1]
        end_index = SEGMENT_ENDS_WORKSPACE[i]
        # NOTE: We switch from 1-based to 0-based indexing since
        #       the data moves from Fortran to Python.
        triples = tuple(
            (
                SEGMENTS_WORKSPACE[j].edge_index - 1,
                SEGMENTS_WORKSPACE[j].start,
                SEGMENTS_WORKSPACE[j].end,
            )
            for j in range(begin_index, end_index)
        )
        curved_polygons.append(triples)

    if not curved_polygons:
        return [], None, ()

    # NOTE: We compute the nodes for each of the six edges. This is a
    #       "wasted" computation / storage.
    num_nodes, dimension = np.shape(nodes1)
    edge_nodes1 = np.empty((degree1 + 1, dimension), order='F')
    edge_nodes2 = np.empty((degree1 + 1, dimension), order='F')
    edge_nodes3 = np.empty((degree1 + 1, dimension), order='F')
    bezier._surface.compute_edge_nodes(
        &num_nodes,
        &dimension,
        &nodes1[0, 0],
        &degree1,
        &edge_nodes1[0, 0],
        &edge_nodes2[0, 0],
        &edge_nodes3[0, 0],
    )

    num_nodes, dimension = np.shape(nodes2)
    edge_nodes4 = np.empty((degree2 + 1, 2), order='F')
    edge_nodes5 = np.empty((degree2 + 1, 2), order='F')
    edge_nodes6 = np.empty((degree2 + 1, 2), order='F')
    bezier._surface.compute_edge_nodes(
        &num_nodes,
        &dimension,
        &nodes2[0, 0],
        &degree2,
        &edge_nodes4[0, 0],
        &edge_nodes5[0, 0],
        &edge_nodes6[0, 0],
    )

    all_edge_nodes = (
        edge_nodes1,
        edge_nodes2,
        edge_nodes3,
        edge_nodes4,
        edge_nodes5,
        edge_nodes6,
    )
    return curved_polygons, None, all_edge_nodes


def _surface_intersections_resize(
        double[::1, :] nodes1, int degree1,
        double[::1, :] nodes2, int degree2,
        int segment_ends_size, int segments_size,
        int num_intersected, int resizes_allowed):
    global SEGMENT_ENDS_WORKSPACE
    cdef int num_segments

    if resizes_allowed > 0:
        if num_intersected > segment_ends_size:
            reset_workspaces(segment_ends_size=num_intersected)
        else:
            num_segments = SEGMENT_ENDS_WORKSPACE[num_intersected - 1]
            reset_workspaces(segments_size=num_segments)

        return surface_intersections(
            nodes1, degree1, nodes2, degree2,
            resizes_allowed=resizes_allowed - 1)
    else:
        if num_intersected > segment_ends_size:
            msg = SEGMENT_ENDS_TOO_SMALL.format(
                num_intersected, segment_ends_size)
            raise ValueError(msg)
        else:
            num_segments = SEGMENT_ENDS_WORKSPACE[num_intersected - 1]
            msg = SEGMENTS_TOO_SMALL.format(
                num_segments, segments_size)
            raise ValueError(msg)


def surface_intersections(
        double[::1, :] nodes1, int degree1,
        double[::1, :] nodes2, int degree2,
        bint verify=True, int resizes_allowed=2):
    # NOTE: ``verify`` is unused. It is only provided as compatibility
    #       with the Python version.
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    cdef int num_nodes1, num_nodes2
    cdef int segment_ends_size
    cdef int segments_size
    cdef int num_intersected
    cdef SurfaceContained contained
    cdef bezier._status.Status status

    # NOTE: We don't check that there are 2 columns.
    num_nodes1, _ = np.shape(nodes1)
    num_nodes2, _ = np.shape(nodes2)

    segment_ends_size, segments_size = workspace_sizes()

    bezier._surface_intersection.surface_intersections(
        &num_nodes1,
        &nodes1[0, 0],
        &degree1,
        &num_nodes2,
        &nodes2[0, 0],
        &degree2,
        &segment_ends_size,
        &SEGMENT_ENDS_WORKSPACE[0],
        &segments_size,
        &SEGMENTS_WORKSPACE[0],
        &num_intersected,
        &contained,
        &status,
    )

    if status == bezier._status.Status.SUCCESS:
        if contained == bezier._surface_intersection.SurfaceContained.FIRST:
            return None, True, ()
        elif contained == bezier._surface_intersection.SurfaceContained.SECOND:
            return None, False, ()
        else:
            # Assumes, but does not check, that ``contained`` is equal to
            # ``bezier._surface_intersection.NEITHER``.
            return _surface_intersections_success(
                nodes1, degree1, nodes2, degree2, num_intersected)
    elif status == bezier._status.Status.INSUFFICIENT_SPACE:
        return _surface_intersections_resize(
            nodes1, degree1, nodes2, degree2,
            segment_ends_size, segments_size,
            num_intersected, resizes_allowed)
    elif status == bezier._status.Status.NO_CONVERGE:
        # NOTE: This assumes, but does not verify, that this comes from
        #       ``curve_intersections()``.
        raise ValueError(
            'Curve intersection failed to converge to approximately linear '
            'subdivisions after 20 iterations.')
    elif status == bezier._status.Status.PARALLEL:
        raise NotImplementedError('Line segments parallel.')
    elif status == bezier._status.Status.WIGGLE_FAIL:
        # NOTE: This branch may not be tested because it's quite difficult to
        #       come up with an example that causes it.
        raise ValueError('outside of unit interval')
    elif status == bezier._status.Status.EDGE_END:
        # NOTE: This text is identical (or should be) to the exception
        #       in the Python ``classify_intersection()``.
        raise ValueError('Intersection occurs at the end of an edge')
    elif status == bezier._status.Status.BAD_TANGENT:
        # NOTE: This text is identical (or should be) to the exception(s)
        #       in the Python ``classify_tangent_intersection()``.
        raise NotImplementedError(
            'Curves moving in opposite direction but define overlapping arcs.')
    elif status == bezier._status.Status.SAME_CURVATURE:
        # NOTE: This text is identical (or should be) to the exception(s)
        #       in the Python ``classify_tangent_intersection()``.
        raise NotImplementedError('Tangent curves have same curvature.')
    elif status == bezier._status.Status.UNKNOWN:
        # NOTE: This assumes that the Fortran ``interior_combine()`` has
        #       failed to return to the start node after ``MAX_EDGES = 10``
        #       edges have been added. As the status name indicates, this
        #       should never occur, hence will be difficult to test.
        raise RuntimeError('Unexpected number of edges', 11)
    else:
        # NOTE: If ``status`` isn't one of the enum values, then it is the
        #       number of candidate intersections.
        raise NotImplementedError(TOO_MANY_TEMPLATE.format(status))


def free_surface_intersections_workspace():
    bezier._surface_intersection.free_surface_intersections_workspace()


def _type_info():
    return (
        SEGMENT_DTYPE.isnative,
        SEGMENT_DTYPE.itemsize,
        SEGMENT_DTYPE.num,
        sizeof(CurvedPolygonSegment),
    )
