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
#include "xtensor/xview.hpp"

int main(int argc, char* argv[])
{
    // Inputs.
    bezier::Matrix<2, 5> nodes {
        { 1.0, 1.25, 1.5, 1.75, 2.0 },
        { 3.0, 3.5, 4.0, 4.5, 5.0 },
    };

    // Outputs.
    int num_reduced_nodes;
    bezier::Matrix<2, 5> reduced;
    bool not_implemented;

    std::tuple<int, bool> status_pair = bezier::full_reduce(nodes, reduced);
    num_reduced_nodes = std::get<0>(status_pair);
    not_implemented = std::get<1>(status_pair);
    auto reduced_view
        = xt::view(reduced, xt::all(), xt::range(0, num_reduced_nodes));

    std::cout << "Number of reduced nodes: " << num_reduced_nodes << std::endl;
    std::cout << "Reduced:" << std::endl;
    std::cout << reduced_view << std::endl;
    std::cout << "Not implemented: " << (not_implemented ? "TRUE" : "FALSE")
              << std::endl;

    return 0;
}
