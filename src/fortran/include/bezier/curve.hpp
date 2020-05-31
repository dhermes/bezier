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

}

#endif /* BEZIER_CURVE_HPP */
