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
    double nodes[10] = { 1.0, 0.0, 0.75, 2.0, 0.5, -2.0, 0.25, 2.0, 0.0, 0.0 };
    double tangent_vec[2] = { -1.0, 0.0 };
    double s = 0.5;

    // Outputs.
    double curvature;

    get_curvature(&num_nodes, nodes, tangent_vec, &s, &curvature);
    printf("Curvature: %f\n", curvature);

    return 0;
}
