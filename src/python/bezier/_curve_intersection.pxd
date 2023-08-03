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

"""Cython wrapper for ``curve_intersection.f90``."""


from libcpp cimport bool as bool_t

from bezier._status cimport Status


cdef extern from "bezier/curve_intersection.h":
    cdef enum BoxIntersectionType:
        INTERSECTION = 0
        TANGENT = 1
        DISJOINT = 2

    void newton_refine_curve_intersect "BEZ_newton_refine_curve_intersect" (
        const double* s, const int* num_nodes1,
        const double* nodes1, const double* t, const int* num_nodes2,
        const double* nodes2, double* new_s, double* new_t, Status* status)
    void bbox_intersect "BEZ_bbox_intersect" (
        const int* num_nodes1, const double* nodes1,
        const int* num_nodes2, const double* nodes2, BoxIntersectionType* enum_)
    void curve_intersections "BEZ_curve_intersections" (
        const int* num_nodes_first,
        const double* nodes_first, const int* num_nodes_second,
        const double* nodes_second, const int* intersections_size,
        double* intersections, const int* num_intersections, bool_t* coincident,
        Status* status)
    void free_curve_intersections_workspace "BEZ_free_curve_intersections_workspace" ()
