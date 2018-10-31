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
    double nodes[6] = { 0.0, 0.0, 0.5, 1.0, 1.0, 0.0 };
    double start = -0.25;
    double end = 0.75;

    // Outputs.
    double new_nodes[6];

    specialize_curve(&num_nodes, &dimension, nodes, &start, &end, new_nodes);
    printf("New Nodes:\n");
    printf("%f, %f, %f\n", new_nodes[0], new_nodes[2], new_nodes[4]);
    printf("%f, %f, %f\n", new_nodes[1], new_nodes[3], new_nodes[5]);

    return 0;
}
