// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     https://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef BEZIER_HELPERS_H
#define BEZIER_HELPERS_H

#include <stdbool.h>

#if defined(__cplusplus)
extern "C" {
#endif

void BEZ_cross_product(const double* vec0, const double* vec1, double* result);
void BEZ_bbox(const int* num_nodes, const double* nodes, double* left,
    double* right, double* bottom, double* top);
void BEZ_wiggle_interval(const double* value, double* result, bool* success);
void BEZ_contains_nd(const int* num_nodes, const int* dimension,
    const double* nodes, const double* point, bool* predicate);
bool BEZ_vector_close(int* num_values, const double* vec1, const double* vec2,
    const double* eps);
bool BEZ_in_interval(
    const double* value, const double* start, const double* end);
void BEZ_simple_convex_hull(
    int* num_points, const double* points, int* polygon_size, double* polygon);
void BEZ_polygon_collide(const int* polygon_size1, const double* polygon1,
    const int* polygon_size2, const double* polygon2, bool* collision);

#if defined(__cplusplus)
}
#endif

#if defined(__cplusplus)
#include <vector>

namespace bezier {
double cross_product(
    const std::vector<double>& vec0, const std::vector<double>& vec1)
{
    double result;
    BEZ_cross_product(vec0.data(), vec1.data(), &result);
    return result;
}
}
#endif

#endif /* BEZIER_HELPERS_H */
