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

#ifndef BEZIER_HELPERS_HPP
#define BEZIER_HELPERS_HPP

#include <array>
#include <tuple>
#include <vector>

#include "bezier/helpers.h"
#include "xtensor/xtensor.hpp"

namespace bezier {

template <size_t M, size_t N>
using Matrix = xt::xtensor_fixed<double, xt::xshape<M, N>,
    xt::layout_type::column_major>;
template <size_t N> using Vector = std::array<double, N>;

double cross_product(
    const std::vector<double>& vec0, const std::vector<double>& vec1)
{
    double result;
    BEZ_cross_product(vec0.data(), vec1.data(), &result);
    return result;
}

template <size_t N>
Vector<4> bbox(const xt::xtensor_fixed<double, xt::xshape<2, N>,
    xt::layout_type::column_major>& nodes)
{
    double left, right, bottom, top;
    int num_nodes = N;
    BEZ_bbox(&num_nodes, nodes.data(), &left, &right, &bottom, &top);
    return { left, right, bottom, top };
}

std::tuple<double, bool> wiggle_interval(const double& value)
{
    double result;
    bool success;
    BEZ_wiggle_interval(&value, &result, &success);
    return std::make_tuple(result, success);
}

template <size_t D, size_t N>
bool contains_nd(const xt::xtensor_fixed<double, xt::xshape<D, N>,
                     xt::layout_type::column_major>& nodes,
    const Vector<D>& point)
{
    int dimension = D;
    int num_nodes = N;
    bool predicate;
    BEZ_contains_nd(
        &num_nodes, &dimension, nodes.data(), point.data(), &predicate);
    return predicate;
}

template <size_t D>
bool vector_close(
    const Vector<D>& vec1, const Vector<D>& vec2, const double& eps)
{
    int num_values = D;
    return BEZ_vector_close(&num_values, vec1.data(), vec2.data(), &eps);
}

bool in_interval(const double& value, const double& start, const double& end)
{
    return BEZ_in_interval(&value, &start, &end);
}

template <size_t N>
int simple_convex_hull(const xt::xtensor_fixed<double, xt::xshape<2, N>,
                           xt::layout_type::column_major>& points,
    xt::xtensor_fixed<double, xt::xshape<2, N>, xt::layout_type::column_major>&
        polygon)
{
    int num_points = N;
    int polygon_size;
    BEZ_simple_convex_hull(
        &num_points, points.data(), &polygon_size, polygon.data());
    return polygon_size;
}

template <size_t N1, size_t N2>
bool polygon_collide(const xt::xtensor_fixed<double, xt::xshape<2, N1>,
                         xt::layout_type::column_major>& polygon1,
    const xt::xtensor_fixed<double, xt::xshape<2, N2>,
        xt::layout_type::column_major>& polygon2)
{
    int polygon_size1 = N1;
    int polygon_size2 = N2;
    bool collision;
    BEZ_polygon_collide(&polygon_size1, polygon1.data(), &polygon_size2,
        polygon2.data(), &collision);
    return collision;
}

}

#endif /* BEZIER_HELPERS_HPP */
