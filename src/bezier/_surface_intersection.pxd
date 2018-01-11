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

"""Cython wrapper for ``surface_intersection.f90``."""


from bezier._status cimport Status


cdef extern from "bezier/surface_intersection.h":
    cpdef enum SurfaceContained:
        NEITHER = 0
        FIRST = 1
        SECOND = 2

    ctypedef struct CurvedPolygonSegment:
        double start
        double end
        int edge_index

    void newton_refine_surface(
        int *num_nodes, double *nodes, int *degree,
        double *x_val, double *y_val, double *s, double *t,
        double *updated_s, double *updated_t)
    void locate_point_surface(
        int *num_nodes, double *nodes, int *degree,
        double *x_val, double *y_val, double *s_val, double *t_val)
    void surface_intersections(
        int *num_nodes1, double *nodes1, int *degree1,
        int *num_nodes2, double *nodes2, int *degree2,
        int *segment_ends_size, int *segment_ends,
        int *segments_size, CurvedPolygonSegment *segments,
        int *num_intersected, SurfaceContained *contained, Status *status)
    void free_surface_intersections_workspace()
