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
    bezier::Matrix<2, 5> nodes {
        { 1.0, 0.75, 0.5, 0.25, 0.0 },
        { 0.0, 2.0, -2.0, 2.0, 0.0 },
    };
    bezier::Vector<2> tangent_vec { -1.0, 0.0 };
    double s = 0.5;

    // Outputs.
    double curvature = bezier::get_curvature(nodes, tangent_vec, s);

    std::cout << "Curvature: " << curvature << std::endl;

    return 0;
}
