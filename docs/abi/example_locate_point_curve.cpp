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
            { 0.0, -1.0, 1.0, -0.75 },
            { 2.0, 0.0, 1.0, 1.625 },
        };
    xt::xtensor_fixed<double, xt::xshape<2>, xt::layout_type::column_major>
        point1 { -0.09375, 0.828125 };
    xt::xtensor_fixed<double, xt::xshape<2>, xt::layout_type::column_major>
        point2 { 0.0, 1.5 };
    xt::xtensor_fixed<double, xt::xshape<2>, xt::layout_type::column_major>
        point3 { -0.25, 1.375 };

    // Outputs.
    double s_approx;

    s_approx = bezier::locate_point_curve(nodes, point1);
    std::cout << "When B(s) = " << point1 << "; s = " << s_approx << std::endl;
    s_approx = bezier::locate_point_curve(nodes, point2);
    std::cout << "When B(s) = " << point2 << "; s = " << s_approx << std::endl;
    s_approx = bezier::locate_point_curve(nodes, point3);
    std::cout << "When B(s) = " << point3 << "; s = " << s_approx << std::endl;

    return 0;
}
