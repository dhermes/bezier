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
cimport bezier._surface_intersection
from bezier._surface_intersection cimport CurvedPolygonSegment


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


def reset_workspaces(int segment_ends_size, int segments_size):
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    SEGMENT_ENDS_WORKSPACE = np.empty(segment_ends_size, dtype=np.intc)
    SEGMENTS_WORKSPACE = np.empty(segments_size, dtype=SEGMENT_DTYPE)


def workspace_sizes():
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    cdef int segment_ends_size
    cdef int segments_size

    segment_ends_size, = np.shape(SEGMENT_ENDS_WORKSPACE)
    segments_size, = np.shape(SEGMENTS_WORKSPACE)
    return segment_ends_size, segments_size


def surface_intersections(
        double[::1, :] nodes1, int degree1,
        double[::1, :] nodes2, int degree2):
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    cdef int num_nodes1, num_nodes2
    cdef int segment_ends_size
    cdef int segments_size
    cdef int num_intersected, contained, status
    cdef int num_segments, i

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

    if status == bezier._status.Status.INSUFFICIENT_SPACE:
        if num_intersected > segment_ends_size:
            return (
                None,
                None,
                num_intersected,
                contained,
                status,
            )
        else:
            return (
                SEGMENT_ENDS_WORKSPACE[:num_intersected],
                None,
                num_intersected,
                contained,
                status,
            )
    else:
        if num_intersected == 0:
            return (
                None,
                None,
                num_intersected,
                contained,
                status,
            )
        else:
            num_segments = SEGMENT_ENDS_WORKSPACE[num_intersected - 1]
            return (
                SEGMENT_ENDS_WORKSPACE[:num_intersected],
                [
                    SEGMENTS_WORKSPACE[i]
                    for i in range(num_segments)
                ],
                num_intersected,
                contained,
                status,
            )
