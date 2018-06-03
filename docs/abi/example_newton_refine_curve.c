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
  double nodes[6] = {0.0, 0.0, 1.0, 2.0, 3.0, 1.0};
  double point[2] = {0.5625, 0.8125};
  double s = 0.75;

  // Outputs.
  double updated_s;

  newton_refine_curve(&num_nodes, &dimension, nodes, point, &s, &updated_s);
  printf("Updated s: %f\n", updated_s);

  return 0;
}
