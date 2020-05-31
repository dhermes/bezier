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
            { -3.0, 0.0, 1.0, 0.0 },
            { 3.0, 2.0, 3.0, 6.0 },
        };

    // Outputs.
    xt::xtensor_fixed<double, xt::xshape<2, 3>, xt::layout_type::column_major>
        reduced;
    bool not_implemented;

    not_implemented = bezier::reduce_pseudo_inverse(nodes, reduced);

    std::cout << "Reduced:" << std::endl;
    std::cout << reduced << std::endl;
    std::cout << "Not implemented: " << (not_implemented ? "TRUE" : "FALSE")
              << std::endl;

    return 0;
}
