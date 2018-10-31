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
    int num_nodes = 3;
    int dimension = 2;
    double nodes[6] = { 0.0, 0.0, 1.25, 3.0, 2.0, 1.0 };

    // Outputs.
    double left_nodes[6];
    double right_nodes[6];

    subdivide_nodes_curve(
        &num_nodes, &dimension, nodes, left_nodes, right_nodes);
    printf("Left Nodes:\n");
    printf("%f, %f, %f\n", left_nodes[0], left_nodes[2], left_nodes[4]);
    printf("%f, %f, %f\n", left_nodes[1], left_nodes[3], left_nodes[5]);
    printf("Right Nodes:\n");
    printf("%f, %f, %f\n", right_nodes[0], right_nodes[2], right_nodes[4]);
    printf("%f, %f, %f\n", right_nodes[1], right_nodes[3], right_nodes[5]);

    return 0;
}
