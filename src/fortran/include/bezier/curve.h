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

#ifndef BEZIER_CURVE_H
#define BEZIER_CURVE_H

#include <stdbool.h>

#if defined(__cplusplus)
extern "C" {
#endif

void BEZ_evaluate_curve_barycentric(int* num_nodes, int* dimension,
    double* nodes, int* num_vals, double* lambda1, double* lambda2,
    double* evaluated);
void BEZ_evaluate_multi(int* num_nodes, int* dimension, double* nodes,
    int* num_vals, double* s_vals, double* evaluated);
void BEZ_specialize_curve(int* num_nodes, int* dimension, double* nodes,
    double* start, double* end, double* new_nodes);
void BEZ_evaluate_hodograph(double* s, int* num_nodes, int* dimension,
    double* nodes, double* hodograph);
void BEZ_subdivide_nodes_curve(int* num_nodes, int* dimension, double* nodes,
    double* left_nodes, double* right_nodes);
void BEZ_newton_refine_curve(int* num_nodes, int* dimension, double* nodes,
    double* point, double* s, double* updated_s);
void BEZ_locate_point_curve(int* num_nodes, int* dimension, double* nodes,
    double* point, double* s_approx);
void BEZ_elevate_nodes_curve(
    int* num_nodes, int* dimension, double* nodes, double* elevated);
void BEZ_get_curvature(int* num_nodes, double* nodes, double* tangent_vec,
    double* s, double* curvature);
void BEZ_reduce_pseudo_inverse(int* num_nodes, int* dimension, double* nodes,
    double* reduced, bool* not_implemented);
void BEZ_full_reduce(int* num_nodes, int* dimension, double* nodes,
    int* num_reduced_nodes, double* reduced, bool* not_implemented);
void BEZ_compute_length(int* num_nodes, int* dimension, double* nodes,
    double* length, int* error_val);

#if defined(__cplusplus)
}
#endif

#endif /* BEZIER_CURVE_H */
