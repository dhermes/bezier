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
    double nodes[8] = { 1.0, 0.0, 1.0, 1.0, 2.0, 0.0, 2.0, 1.0 };
    int num_vals = 3;
    double s_vals[3] = { 0.0, 0.5, 1.0 };

    // Outputs.
    double evaluated[6];

    BEZ_evaluate_multi(
        &num_nodes, &dimension, nodes, &num_vals, s_vals, evaluated);
    printf("Evaluated:\n");
    printf("%f, %f, %f\n", evaluated[0], evaluated[2], evaluated[4]);
    printf("%f, %f, %f\n", evaluated[1], evaluated[3], evaluated[5]);

    return 0;
}
