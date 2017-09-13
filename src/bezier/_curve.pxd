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

"""Cython wrapper for ``curve.f90``."""


cdef extern from "bezier/curve.h":
    void evaluate_curve_barycentric(
        int *degree, int *dimension, double *nodes, int *num_vals,
        double *lambda1, double *lambda2, double *evaluated)
    void evaluate_multi(
        int *degree, int *dimension, double *nodes,
        int *num_vals, double *s_vals, double *evaluated)
    void specialize_curve_generic(
        int *degree, int *dimension, double *nodes,
        double *start, double *end, double *new_nodes)
    void specialize_curve_quadratic(
        int *dimension, double *nodes, double *start,
        double *end, double *new_nodes)
    void specialize_curve(
        int *degree, int *dimension, double *nodes, double *start, double *end,
        double *curve_start, double *curve_end, double *new_nodes,
        double *true_start, double *true_end)
    void evaluate_hodograph(
        double *s, int *degree, int *dimension, double *nodes,
        double *hodograph)
    void subdivide_nodes(
        int *num_nodes, int *dimension, double *nodes,
        double *left_nodes, double *right_nodes)
    void newton_refine(
        int *num_nodes, int *dimension, double *nodes,
        double *point, double *s, double *updated_s)
