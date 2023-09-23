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

#include <tuple>

#include "bezier/curve.h"
#include "bezier/helpers.hpp"
#include "xtensor/xtensor.hpp"

namespace bezier {

template <std::size_t N, std::size_t D, std::size_t k>
Matrix<D, k> evaluate_curve_barycentric(const Matrix<D, N>& nodes,
    const Vector<k>& lambda1, const Vector<k>& lambda2)
{
    Matrix<D, k> evaluated;

    int num_nodes = N;
    int dimension = D;
    int num_vals = k;
    BEZ_evaluate_curve_barycentric(&num_nodes, &dimension, nodes.data(),
        &num_vals, lambda1.data(), lambda2.data(), evaluated.data());
    return evaluated;
}

template <std::size_t N, std::size_t D, std::size_t k>
Matrix<D, k> evaluate_multi(const Matrix<D, N>& nodes, const Vector<k>& s_vals)
{
    Matrix<D, k> evaluated;

    int num_nodes = N;
    int dimension = D;
    int num_vals = k;
    BEZ_evaluate_multi(&num_nodes, &dimension, nodes.data(), &num_vals,
        s_vals.data(), evaluated.data());
    return evaluated;
}

template <std::size_t N, std::size_t D>
Matrix<D, N> specialize_curve(
    const Matrix<D, N>& nodes, const double& start, const double& end)
{
    Matrix<D, N> new_nodes;

    int num_nodes = N;
    int dimension = D;
    BEZ_specialize_curve(
        &num_nodes, &dimension, nodes.data(), &start, &end, new_nodes.data());
    return new_nodes;
}

template <std::size_t N, std::size_t D>
Vector<D> evaluate_hodograph(const double& s, const Matrix<D, N>& nodes)
{
    Vector<D> hodograph;

    int num_nodes = N;
    int dimension = D;
    BEZ_evaluate_hodograph(
        &s, &num_nodes, &dimension, nodes.data(), hodograph.data());
    return hodograph;
}

template <std::size_t N, std::size_t D>
void subdivide_nodes_curve(const Matrix<D, N>& nodes, Matrix<D, N>& left_nodes,
    Matrix<D, N>& right_nodes)
{
    int num_nodes = N;
    int dimension = D;
    BEZ_subdivide_nodes_curve(&num_nodes, &dimension, nodes.data(),
        left_nodes.data(), right_nodes.data());
}

template <std::size_t N, std::size_t D>
double newton_refine_curve(
    const Matrix<D, N>& nodes, const Vector<D>& point, const double& s)
{
    int num_nodes = N;
    int dimension = D;
    double updated_s;
    BEZ_newton_refine_curve(
        &num_nodes, &dimension, nodes.data(), point.data(), &s, &updated_s);
    return updated_s;
}

template <std::size_t N, std::size_t D>
double locate_point_curve(const Matrix<D, N>& nodes, const Vector<D>& point)
{
    int num_nodes = N;
    int dimension = D;
    double s_approx;
    BEZ_locate_point_curve(
        &num_nodes, &dimension, nodes.data(), point.data(), &s_approx);
    return s_approx;
}

template <std::size_t N, std::size_t D>
Matrix<D, N + 1> elevate_nodes_curve(const Matrix<D, N>& nodes)
{
    Matrix<D, N + 1> elevated;

    int num_nodes = N;
    int dimension = D;
    BEZ_elevate_nodes_curve(
        &num_nodes, &dimension, nodes.data(), elevated.data());
    return elevated;
}

template <std::size_t N>
double get_curvature(
    const Matrix<2, N>& nodes, const Vector<2>& tangent_vec, const double& s)
{
    int num_nodes = N;
    double curvature;
    BEZ_get_curvature(
        &num_nodes, nodes.data(), tangent_vec.data(), &s, &curvature);
    return curvature;
}

template <std::size_t N, std::size_t D>
bool reduce_pseudo_inverse(
    const Matrix<D, N>& nodes, Matrix<D, N - 1>& reduced)
{
    int num_nodes = N;
    int dimension = D;
    bool not_implemented;
    BEZ_reduce_pseudo_inverse(&num_nodes, &dimension, nodes.data(),
        reduced.data(), &not_implemented);
    return not_implemented;
}

template <std::size_t N, std::size_t D>
std::tuple<int, bool> full_reduce(
    const Matrix<D, N>& nodes, Matrix<D, N>& reduced)
{
    int num_nodes = N;
    int dimension = D;
    int num_reduced_nodes;
    bool not_implemented;
    BEZ_full_reduce(&num_nodes, &dimension, nodes.data(), &num_reduced_nodes,
        reduced.data(), &not_implemented);
    return std::make_tuple(num_reduced_nodes, not_implemented);
}

template <std::size_t N, std::size_t D>
std::tuple<double, int> compute_length(const Matrix<D, N>& nodes)
{
    int num_nodes = N;
    int dimension = D;
    double length;
    int error_val;
    BEZ_compute_length(
        &num_nodes, &dimension, nodes.data(), &length, &error_val);
    return std::make_tuple(length, error_val);
}

}

#endif /* BEZIER_CURVE_HPP */
