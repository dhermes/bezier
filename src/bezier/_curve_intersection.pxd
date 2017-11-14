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


cdef extern from "bezier/curve_intersection.h":
    cpdef enum BoxIntersectionType:
        INTERSECTION = 0
        TANGENT = 1
        DISJOINT = 2

    cpdef enum AllIntersectionsStatus:
        SUCCESS = 0
        NO_CONVERGE = 1
        TOO_SMALL = 2
        PARALLEL = 3
        WIGGLE_FAIL = 4
        UNKNOWN = 5

    void linearization_error(
        int *num_nodes, int *dimension, double *nodes, double *error)
    void segment_intersection(
        double *start0, double *end0, double *start1, double *end1,
        double *s, double *t, bool_t *success)
    void newton_refine_intersect(
        double *s, int *num_nodes1, double *nodes1,
        double *t, int *num_nodes2, double *nodes2,
        double *new_s, double *new_t)
    void bbox_intersect(
        int *num_nodes1, double *nodes1,
        int *num_nodes2, double *nodes2, int *enum_)
    void parallel_different(
        double *start0, double *end0,
        double *start1, double *end1, bool_t *result)
    void from_linearized(
        double *error1, double *start1, double *end1, double *start_node1,
        double *end_node1, int *num_nodes1, double *root_nodes1,
        double *error2, double *start2, double *end2, double *start_node2,
        double *end_node2, int *num_nodes2, double *root_nodes2,
        double *refined_s, double *refined_t,
        bool_t *does_intersect, int *py_exc)
    void bbox_line_intersect(
        int *num_nodes, double *nodes, double *line_start, double *line_end,
        int *enum_)
    void curve_intersections(
        int *num_nodes_first, double *nodes_first,
        int *num_nodes_second, double *nodes_second,
        int *intersections_size, double *intersections,
        int *num_intersections, int *status)
    void free_curve_intersections_workspace()
