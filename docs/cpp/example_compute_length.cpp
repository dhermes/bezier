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
    bezier::Matrix<2, 2> nodes {
        { 0.0, 3.0 },
        { 0.0, 4.0 },
    };
    std::tuple<int, bool> status_pair = bezier::compute_length(nodes);

    std::cout << "Length: " << std::get<0>(status_pair) << std::endl;
    std::cout << "Error value: " << std::get<1>(status_pair) << std::endl;

    return 0;
}
