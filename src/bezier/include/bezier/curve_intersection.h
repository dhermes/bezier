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

#include "bezier/_bool_patch.h"
#include "bezier/status.h"

#if defined (__cplusplus)
extern "C" {
#endif

typedef enum BoxIntersectionType {
  INTERSECTION = 0,
  TANGENT = 1,
  DISJOINT = 2,
} BoxIntersectionType;

void newton_refine_curve_intersect(
    double *s, int *num_nodes1, double *nodes1,
    double *t, int *num_nodes2, double *nodes2,
    double *new_s, double *new_t, Status *status);
void bbox_intersect(
    int *num_nodes1, double *nodes1,
    int *num_nodes2, double *nodes2, BoxIntersectionType *enum_);
void curve_intersections(
    int *num_nodes_first, double *nodes_first,
    int *num_nodes_second, double *nodes_second,
    int *intersections_size, double *intersections,
    int *num_intersections, bool *coincident, Status *status);
void set_max_candidates(int *num_candidates);
void get_max_candidates(int *num_candidates);
void free_curve_intersections_workspace(void);

#if defined (__cplusplus)
}
#endif

#endif /* BEZIER_CURVE_INTERSECTION_H */
