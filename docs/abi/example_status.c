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
    int num_nodes1 = 3;
    double nodes1[6] = { 0.0, 0.0, 2.0, 4.0, 4.0, 0.0 };
    double s = 7.0 / 8.0;
    int num_nodes2 = 2;
    double nodes2[4] = { 2.0, 0.0, 0.0, 3.0 };
    double t = -2.0;

    // Outputs.
    double new_s, new_t;
    Status status;

    newton_refine_curve_intersect(&s, &num_nodes1, nodes1, &t, &num_nodes2,
        nodes2, &new_s, &new_t, &status);
    if (status == SINGULAR) {
        printf("Jacobian is singular.\n");
    }

    return 0;
}
