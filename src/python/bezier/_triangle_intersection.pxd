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

"""Cython wrapper for ``triangle_intersection.f90``."""


from bezier._status cimport Status


cdef extern from "bezier/triangle_intersection.h":
    cdef enum TriangleContained:
        NEITHER = 0
        FIRST = 1
        SECOND = 2

    ctypedef struct CurvedPolygonSegment:
        double start
        double end
        int edge_index

    void newton_refine_triangle "BEZ_newton_refine_triangle" (
        const int* num_nodes, const double* nodes,
        const int* degree, const double* x_val, const double* y_val,
        const double* s, const double* t, double* updated_s, double* updated_t)
    void locate_point_triangle "BEZ_locate_point_triangle" (
        const int* num_nodes, const double* nodes,
        const int* degree, const double* x_val, const double* y_val, double* s_val,
        double* t_val)
    void triangle_intersections "BEZ_triangle_intersections" (
        const int* num_nodes1, const double* nodes1,
        const int* degree1, const int* num_nodes2, const double* nodes2,
        const int* degree2, const int* segment_ends_size, int* segment_ends,
        int* segments_size, CurvedPolygonSegment* segments, int* num_intersected,
        TriangleContained* contained, Status* status)
    void free_triangle_intersections_workspace "BEZ_free_triangle_intersections_workspace" ()
