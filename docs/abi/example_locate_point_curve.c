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

#include "bezier.h"
#include <stdio.h>

int main(void)
{
    // Inputs.
    int num_nodes = 4;
    int dimension = 2;
    double nodes[8] = { 0.0, 2.0, -1.0, 0.0, 1.0, 1.0, -0.75, 1.625 };
    double point1[2] = { -0.09375, 0.828125 };
    double point2[2] = { 0.0, 1.5 };
    double point3[2] = { -0.25, 1.375 };

    // Outputs.
    double s_approx;

    locate_point_curve(&num_nodes, &dimension, nodes, point1, &s_approx);
    printf("When B(s) = [% f, %f]; s = % f\n", point1[0], point1[1], s_approx);
    locate_point_curve(&num_nodes, &dimension, nodes, point2, &s_approx);
    printf("When B(s) = [% f, %f]; s = % f\n", point2[0], point2[1], s_approx);
    locate_point_curve(&num_nodes, &dimension, nodes, point3, &s_approx);
    printf("When B(s) = [% f, %f]; s = % f\n", point3[0], point3[1], s_approx);

    return 0;
}
