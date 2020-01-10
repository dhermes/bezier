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
    double nodes[8] = { -3.0, 3.0, 0.0, 2.0, 1.0, 3.0, 0.0, 6.0 };

    // Outputs.
    double reduced[6];
    bool not_implemented;

    BEZ_reduce_pseudo_inverse(
        &num_nodes, &dimension, nodes, reduced, &not_implemented);
    printf("Reduced:\n");
    printf("% f, %f, %f\n", reduced[0], reduced[2], reduced[4]);
    printf("% f, %f, %f\n", reduced[1], reduced[3], reduced[5]);
    printf("Not implemented: %s\n", not_implemented ? "TRUE" : "FALSE");

    return 0;
}
