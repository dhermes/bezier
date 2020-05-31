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

#include <iostream>

#include "bezier.hpp"
#include "xtensor/xfixed.hpp"
#include "xtensor/xio.hpp"

int main(int argc, char* argv[])
{
    // Inputs.
    xt::xtensor_fixed<double, xt::xshape<2, 3>, xt::layout_type::column_major>
        nodes {
            { 0.0, 2.0, 3.0 },
            { 1.0, 1.0, 3.0 },
        };
    xt::xtensor_fixed<double, xt::xshape<4>, xt::layout_type::column_major>
        lambda1 { 0.25, 0.5, 0.0, 1.0 };
    xt::xtensor_fixed<double, xt::xshape<4>, xt::layout_type::column_major>
        lambda2 { 0.75, 0.25, 0.5, 0.25 };

    // Outputs.
    auto evaluated
        = bezier::evaluate_curve_barycentric(nodes, lambda1, lambda2);

    std::cout << "Evaluated:" << std::endl;
    std::cout << evaluated << std::endl;

    return 0;
}
