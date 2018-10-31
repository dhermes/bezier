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
    double nodes[6] = { 0.0, 0.0, 1.5, 1.5, 3.0, 0.0 };

    // Outputs.
    double elevated[8];

    elevate_nodes_curve(&num_nodes, &dimension, nodes, elevated);
    printf("Elevated:\n");
    printf("%f, %f, %f, %f\n", elevated[0], elevated[2], elevated[4],
        elevated[6]);
    printf("%f, %f, %f, %f\n", elevated[1], elevated[3], elevated[5],
        elevated[7]);

    return 0;
}
