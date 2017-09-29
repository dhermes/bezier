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

#if defined (__cplusplus)
extern "C" {
#elif !defined (_MSC_VER)
#include <stdbool.h>
#endif

#if !defined (__cplusplus) && defined (_MSC_VER)
#  if !defined (FALSE)
#    define FALSE (0)
#  endif
#  if !defined (TRUE)
#    define TRUE (!FALSE)
#  endif
#  if _MSC_VER >= 1800
#    include <stdbool.h>
#  else
#    define bool int
#    define true TRUE
#    define false FALSE
#  endif
#endif

enum BoxIntersectionType {
  INTERSECTION = 0,
  TANGENT = 1,
  DISJOINT = 2,
};

void linearization_error(
    int *num_nodes, int *dimension, double *nodes, double *error);
void segment_intersection(
    double *start0, double *end0, double *start1, double *end1,
    double *s, double *t, bool *success);
void newton_refine_intersect(
    double *s, int *degree1, double *nodes1,
    double *t, int *degree2, double *nodes2,
    double *new_s, double *new_t);
void bbox_intersect(
    int *num_nodes1, double *nodes1,
    int *num_nodes2, double *nodes2, int *enum_);
void parallel_different(
    double *start0, double *end0,
    double *start1, double *end1, bool *result);
void from_linearized(
    double *error1, double *start1, double *end1,
    double *start_node1, double *end_node1, int *degree1, double *nodes1,
    double *error2, double *start2, double *end2,
    double *start_node2, double *end_node2, int *degree2, double *nodes2,
    double *refined_s, double *refined_t,
    bool *does_intersect, int *py_exc);
void bbox_line_intersect(
    int *num_nodes, double *nodes, double *line_start, double *line_end,
    int *enum_);

#if defined (__cplusplus)
}
#endif

#endif /* BEZIER_CURVE_INTERSECTION_H */
