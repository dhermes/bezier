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

int main(int argc, char* argv[])
{
    // Inputs.
    double s = 0.125;
    xt::xtensor_fixed<double, xt::xshape<2, 4>, xt::layout_type::column_major>
        nodes {
            { 1.0, 1.0, 2.0, 2.0 },
            { 0.0, 1.0, 0.0, 1.0 },
        };

    // Outputs.
    std::array<double, 2> hodograph = bezier::evaluate_hodograph(s, nodes);

    std::cout << "Hodograph:" << std::endl;
    std::cout << hodograph[0] << std::endl;
    std::cout << hodograph[1] << std::endl;

    return 0;
}
