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

void BEZ_evaluate_curve_barycentric(const int* num_nodes, const int* dimension,
    const double* nodes, const int* num_vals, const double* lambda1,
    const double* lambda2, double* evaluated);
void BEZ_evaluate_multi(const int* num_nodes, const int* dimension,
    const double* nodes, const int* num_vals, const double* s_vals,
    double* evaluated);
void BEZ_specialize_curve(const int* num_nodes, const int* dimension,
    const double* nodes, const double* start, const double* end,
    double* new_nodes);
void BEZ_evaluate_hodograph(const double* s, const int* num_nodes,
    const int* dimension, const double* nodes, double* hodograph);
void BEZ_subdivide_nodes_curve(const int* num_nodes, const int* dimension,
    const double* nodes, double* left_nodes, double* right_nodes);
void BEZ_newton_refine_curve(const int* num_nodes, const int* dimension,
    const double* nodes, const double* point, const double* s,
    double* updated_s);
void BEZ_locate_point_curve(const int* num_nodes, const int* dimension,
    const double* nodes, const double* point, double* s_approx);
void BEZ_elevate_nodes_curve(const int* num_nodes, const int* dimension,
    const double* nodes, double* elevated);
void BEZ_get_curvature(const int* num_nodes, const double* nodes,
    const double* tangent_vec, const double* s, double* curvature);
void BEZ_reduce_pseudo_inverse(const int* num_nodes, const int* dimension,
    const double* nodes, double* reduced, bool* not_implemented);
void BEZ_full_reduce(const int* num_nodes, const int* dimension,
    const double* nodes, int* num_reduced_nodes, double* reduced,
    bool* not_implemented);
void BEZ_compute_length(const int* num_nodes, const int* dimension,
    const double* nodes, double* length, int* error_val);

#if defined(__cplusplus)
}
#endif

#endif /* BEZIER_CURVE_H */
