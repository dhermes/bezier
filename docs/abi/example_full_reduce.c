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
    int num_nodes = 5;
    int dimension = 2;
    double nodes[10] = { 1.0, 3.0, 1.25, 3.5, 1.5, 4.0, 1.75, 4.5, 2.0, 5.0 };

    // Outputs.
    int num_reduced_nodes;
    double reduced[10];
    bool not_implemented;

    full_reduce(&num_nodes, &dimension, nodes, &num_reduced_nodes, reduced,
        &not_implemented);
    printf("Number of reduced nodes: %d\n", num_reduced_nodes);
    printf("Reduced:\n");
    printf("%f, %f\n", reduced[0], reduced[2]);
    printf("%f, %f\n", reduced[1], reduced[3]);
    printf("Not implemented: %s\n", not_implemented ? "TRUE" : "FALSE");

    return 0;
}
