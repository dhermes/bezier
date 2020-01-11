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

"""Cython wrapper for ``helpers.f90``."""


from libcpp cimport bool as bool_t


cdef extern from "bezier/helpers.h":
    void cross_product "BEZ_cross_product" (
        const double* vec0, const double* vec1, double* result)
    void bbox "BEZ_bbox" (
        const int* num_nodes, const double* nodes, double* left,
        double* right, double* bottom, double* top)
    void wiggle_interval "BEZ_wiggle_interval" (
        const double* value, double* result, bool_t* success)
    void contains_nd "BEZ_contains_nd" (
        const int* num_nodes, const int* dimension,
        const double* nodes, const double* point, bool_t* predicate)
    bool_t vector_close "BEZ_vector_close" (
        int* num_values, const double* vec1, const double* vec2,
        const double* eps)
    bool_t in_interval "BEZ_in_interval" (
        const double* value, const double* start, const double* end)
    void simple_convex_hull "BEZ_simple_convex_hull" (
        int* num_points, const double* points, int* polygon_size, double* polygon)
    void polygon_collide "BEZ_polygon_collide" (
        const int* polygon_size1, const double* polygon1,
        const int* polygon_size2, const double* polygon2, bool_t* collision)
