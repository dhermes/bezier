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

#ifndef BEZIER_CURVE_INTERSECTION_H
#define BEZIER_CURVE_INTERSECTION_H

#include "bezier/status.h"
#include <stdbool.h>

#if defined(__cplusplus)
extern "C" {
#endif

typedef enum BoxIntersectionType {
    INTERSECTION = 0,
    TANGENT = 1,
    DISJOINT = 2,
} BoxIntersectionType;

void BEZ_newton_refine_curve_intersect(const double* s, const int* num_nodes1,
    const double* nodes1, const double* t, const int* num_nodes2,
    const double* nodes2, double* new_s, double* new_t, Status* status);
void BEZ_bbox_intersect(const int* num_nodes1, const double* nodes1,
    const int* num_nodes2, const double* nodes2, BoxIntersectionType* enum_);
void BEZ_curve_intersections(const int* num_nodes_first,
    const double* nodes_first, const int* num_nodes_second,
    const double* nodes_second, const int* intersections_size,
    double* intersections, const int* num_intersections, bool* coincident,
    Status* status);
void BEZ_free_curve_intersections_workspace(void);

#if defined(__cplusplus)
}
#endif

#endif /* BEZIER_CURVE_INTERSECTION_H */
