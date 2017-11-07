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

#if defined (__cplusplus)
extern "C" {
#endif

void newton_refine_surface(
    int *num_nodes, double *nodes, int *degree,
    double *x_val, double *y_val, double *s, double *t,
    double *updated_s, double *updated_t);
void locate_point_surface(
    int *num_nodes, double *nodes, int *degree,
    double *x_val, double *y_val, double *s_val, double *t_val);

#if defined (__cplusplus)
}
#endif

#endif /* BEZIER_SURFACE_INTERSECTION_H */
