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
    bezier::Matrix<2, 2> nodes {
        { 0.0, 3.0 },
        { 0.0, 4.0 },
    };

    // Outputs.
    double length;
    int error_val;

    std::tuple<int, bool> status_pair = bezier::compute_length(nodes);
    length = std::get<0>(status_pair);
    error_val = std::get<1>(status_pair);

    std::cout << "Length: " << length << std::endl;
    std::cout << "Error value: " << error_val << std::endl;

    return 0;
}
