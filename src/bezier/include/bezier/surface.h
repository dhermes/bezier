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

#ifndef BEZIER_SURFACE_H
#define BEZIER_SURFACE_H

#if defined (__cplusplus)
extern "C" {
#endif

void de_casteljau_one_round(
    int *num_nodes, int *dimension, double *nodes, int *degree,
    double *lambda1, double *lambda2, double *lambda3, double *new_nodes);
void evaluate_barycentric(
    int *num_nodes, int *dimension, double *nodes, int *degree,
    double *lambda1, double *lambda2, double *lambda3, double *point);
void evaluate_barycentric_multi(
    int *num_nodes, int *dimension, double *nodes, int *degree,
    int *num_vals, double *param_vals, double *evaluated);
void evaluate_cartesian_multi(
    int *num_nodes, int *dimension, double *nodes, int *degree,
    int *num_vals, double *param_vals, double *evaluated);
void jacobian_both(
    int *num_nodes, int *dimension, double *nodes,
    int *degree, double *new_nodes);
void jacobian_det(
    int *num_nodes, double *nodes, int *degree,
    int *num_vals, double *param_vals, double *evaluated);
void specialize_surface(
    int *num_nodes, int *dimension, double *nodes, int *degree,
    double *weights_a, double *weights_b, double *weights_c,
    double *specialized);
void subdivide_nodes_surface(
    int *num_nodes, int *dimension, double *nodes, int *degree,
    double *nodes_a, double *nodes_b, double *nodes_c, double *nodes_d);
void compute_edge_nodes(
    int *num_nodes, int *dimension, double *nodes, int *degree,
    double *nodes1, double *nodes2, double *nodes3);

#if defined (__cplusplus)
}
#endif

#endif /* BEZIER_SURFACE_H */
