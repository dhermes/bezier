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
  int num_nodes = 2;
  int dimension = 2;
  double nodes[4] = {0.0, 0.0, 3.0, 4.0};

  // Outputs.
  double length;
  int error_val;

  compute_length(&num_nodes, &dimension, nodes, &length, &error_val);
  printf("Length: %f\n", length);

  return 0;
}
