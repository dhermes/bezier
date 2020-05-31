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
    xt::xtensor_fixed<double, xt::xshape<2, 4>, xt::layout_type::column_major>
        nodes {
            { 1.0, 1.0, 2.0, 2.0 },
            { 0.0, 1.0, 0.0, 1.0 },
        };
    xt::xtensor_fixed<double, xt::xshape<3>, xt::layout_type::column_major>
        s_vals { 0.0, 0.5, 1.0 };

    // Outputs.
    auto evaluated = bezier::evaluate_multi(nodes, s_vals);

    std::cout << "Evaluated:" << std::endl;
    std::cout << evaluated << std::endl;

    return 0;
}
