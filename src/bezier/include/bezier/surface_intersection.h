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

#ifndef BEZIER_SURFACE_INTERSECTION_H
#define BEZIER_SURFACE_INTERSECTION_H

#include "bezier/status.h"

#if defined (__cplusplus)
extern "C" {
#endif

typedef enum SurfaceContained {
  NEITHER = 0,
  FIRST = 1,
  SECOND = 2,
} SurfaceContained;

typedef struct CurvedPolygonSegment {
  double start;
  double end;
  int edge_index;
} CurvedPolygonSegment;

void newton_refine_surface(
    int *num_nodes, double *nodes, int *degree,
    double *x_val, double *y_val, double *s, double *t,
    double *updated_s, double *updated_t);
void locate_point_surface(
    int *num_nodes, double *nodes, int *degree,
    double *x_val, double *y_val, double *s_val, double *t_val);
void surface_intersections(
    int *num_nodes1, double *nodes1, int *degree1,
    int *num_nodes2, double *nodes2, int *degree2,
    int *segment_ends_size, int *segment_ends,
    int *segments_size, CurvedPolygonSegment *segments,
    int *num_intersected, SurfaceContained *contained, Status *status);
void free_surface_intersections_workspace(void);

#if defined (__cplusplus)
}
#endif

#endif /* BEZIER_SURFACE_INTERSECTION_H */
