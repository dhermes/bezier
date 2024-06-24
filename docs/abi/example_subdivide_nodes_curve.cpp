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
    bezier::Matrix<2, 3> nodes {
        { 0.0, 1.25, 2.0 },
        { 0.0, 3.0, 1.0 },
    };

    // Outputs.
    bezier::Matrix<2, 3> left_nodes;
    bezier::Matrix<2, 3> right_nodes;

    bezier::subdivide_nodes_curve(nodes, left_nodes, right_nodes);

    std::cout << "Left Nodes:" << std::endl;
    std::cout << left_nodes << std::endl;
    std::cout << "Right Nodes:" << std::endl;
    std::cout << right_nodes << std::endl;

    return 0;
}
