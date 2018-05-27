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

#include <stdio.h>
#include "bezier.h"

int main(void)
{
  // Inputs.
  int num_nodes = 3;
  int dimension = 2;
  double nodes[6] = {0.0, 1.0, 2.0, 1.0, 3.0, 3.0};
  int num_vals = 4;
  double lambda1[4] = {0.25, 0.5, 0.0, 1.0};
  double lambda2[4] = {0.75, 0.25, 0.5, 0.25};

  // Outputs.
  double evaluated[8];

  evaluate_curve_barycentric(
      &num_nodes, &dimension, nodes,
      &num_vals, lambda1, lambda2, evaluated);
  printf("Evaluated:\n");
  printf(
      "%f, %f, %f, %f\n",
      evaluated[0], evaluated[2], evaluated[4], evaluated[6]);
  printf(
      "%f, %f, %f, %f\n",
      evaluated[1], evaluated[3], evaluated[5], evaluated[7]);

  return 0;
}
