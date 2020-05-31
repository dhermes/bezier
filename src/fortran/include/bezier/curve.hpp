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

#ifndef BEZIER_CURVE_HPP
#define BEZIER_CURVE_HPP

#include "bezier/curve.h"
#include "xtensor/xtensor.hpp"

namespace bezier {

template <size_t N, size_t D, size_t k>
xt::xtensor_fixed<double, xt::xshape<D, k>, xt::layout_type::column_major>
evaluate_curve_barycentric(const xt::xtensor_fixed<double, xt::xshape<D, N>,
                               xt::layout_type::column_major>& nodes,
    const xt::xtensor_fixed<double, xt::xshape<k>,
        xt::layout_type::column_major>& lambda1,
    const xt::xtensor_fixed<double, xt::xshape<k>,
        xt::layout_type::column_major>& lambda2)
{
    xt::xtensor_fixed<double, xt::xshape<D, k>, xt::layout_type::column_major>
        evaluated;
    int num_nodes = N;
    int dimension = D;
    int num_vals = k;
    BEZ_evaluate_curve_barycentric(&num_nodes, &dimension, nodes.data(),
        &num_vals, lambda1.data(), lambda2.data(), evaluated.data());
    return evaluated;
}

template <size_t N, size_t D, size_t k>
xt::xtensor_fixed<double, xt::xshape<D, k>, xt::layout_type::column_major>
evaluate_multi(const xt::xtensor_fixed<double, xt::xshape<D, N>,
                   xt::layout_type::column_major>& nodes,
    const xt::xtensor_fixed<double, xt::xshape<k>,
        xt::layout_type::column_major>& s_vals)
{
    xt::xtensor_fixed<double, xt::xshape<D, k>, xt::layout_type::column_major>
        evaluated;
    int num_nodes = N;
    int dimension = D;
    int num_vals = k;
    BEZ_evaluate_multi(&num_nodes, &dimension, nodes.data(), &num_vals,
        s_vals.data(), evaluated.data());
    return evaluated;
}

template <size_t N, size_t D>
xt::xtensor_fixed<double, xt::xshape<D, N>, xt::layout_type::column_major>
specialize_curve(const xt::xtensor_fixed<double, xt::xshape<D, N>,
                     xt::layout_type::column_major>& nodes,
    const double& start, const double& end)
{
    xt::xtensor_fixed<double, xt::xshape<D, N>, xt::layout_type::column_major>
        new_nodes;
    int num_nodes = N;
    int dimension = D;
    BEZ_specialize_curve(
        &num_nodes, &dimension, nodes.data(), &start, &end, new_nodes.data());
    return new_nodes;
}

template <size_t N, size_t D>
xt::xtensor_fixed<double, xt::xshape<D>, xt::layout_type::column_major>
evaluate_hodograph(const double& s,
    const xt::xtensor_fixed<double, xt::xshape<D, N>,
        xt::layout_type::column_major>& nodes)
{
    xt::xtensor_fixed<double, xt::xshape<D>, xt::layout_type::column_major>
        hodograph;
    int num_nodes = N;
    int dimension = D;
    BEZ_evaluate_hodograph(
        &s, &num_nodes, &dimension, nodes.data(), hodograph.data());
    return hodograph;
}

template <size_t N, size_t D>
void subdivide_nodes_curve(const xt::xtensor_fixed<double, xt::xshape<D, N>,
                               xt::layout_type::column_major>& nodes,
    xt::xtensor_fixed<double, xt::xshape<D, N>, xt::layout_type::column_major>&
        left_nodes,
    xt::xtensor_fixed<double, xt::xshape<D, N>, xt::layout_type::column_major>&
        right_nodes)
{
    int num_nodes = N;
    int dimension = D;
    BEZ_subdivide_nodes_curve(&num_nodes, &dimension, nodes.data(),
        left_nodes.data(), right_nodes.data());
}

template <size_t N, size_t D>
double newton_refine_curve(const xt::xtensor_fixed<double, xt::xshape<D, N>,
                               xt::layout_type::column_major>& nodes,
    const xt::xtensor_fixed<double, xt::xshape<D>,
        xt::layout_type::column_major>& point,
    const double& s)
{
    int num_nodes = N;
    int dimension = D;
    double updated_s;
    BEZ_newton_refine_curve(
        &num_nodes, &dimension, nodes.data(), point.data(), &s, &updated_s);
    return updated_s;
}

template <size_t N, size_t D>
double locate_point_curve(const xt::xtensor_fixed<double, xt::xshape<D, N>,
                              xt::layout_type::column_major>& nodes,
    const xt::xtensor_fixed<double, xt::xshape<D>,
        xt::layout_type::column_major>& point)
{
    int num_nodes = N;
    int dimension = D;
    double s_approx;
    BEZ_locate_point_curve(
        &num_nodes, &dimension, nodes.data(), point.data(), &s_approx);
    return s_approx;
}

template <size_t N, size_t D>
xt::xtensor_fixed<double, xt::xshape<D, N + 1>, xt::layout_type::column_major>
elevate_nodes_curve(const xt::xtensor_fixed<double, xt::xshape<D, N>,
    xt::layout_type::column_major>& nodes)
{
    xt::xtensor_fixed<double, xt::xshape<D, N + 1>,
        xt::layout_type::column_major>
        elevated;

    int num_nodes = N;
    int dimension = D;
    BEZ_elevate_nodes_curve(
        &num_nodes, &dimension, nodes.data(), elevated.data());
    return elevated;
}

template <size_t N>
double get_curvature(const xt::xtensor_fixed<double, xt::xshape<2, N>,
                         xt::layout_type::column_major>& nodes,
    const std::array<double, 2>& tangent_vec, const double& s)
{
    int num_nodes = N;
    double curvature;
    BEZ_get_curvature(
        &num_nodes, nodes.data(), tangent_vec.data(), &s, &curvature);
    return curvature;
}

// TODO: reduce_pseudo_inverse
// TODO: full_reduce
// TODO: compute_length

}

#endif /* BEZIER_CURVE_HPP */
