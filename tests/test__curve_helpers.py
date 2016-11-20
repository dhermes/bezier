# See the License for the specific language governing permissions and
# limitations under the License.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest

import numpy as np


class Test_evaluate_multi(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, s_vals):
        from bezier import _curve_helpers

        return _curve_helpers.evaluate_multi(nodes, degree, s_vals)

    def test_it(self):
        import six

        s_vals = np.array([0.0, 0.25, 0.75, 1.0])
        degree = 2
        # B(s) = [s(4 - s), 2s(2s - 1)]
        nodes = np.array([
            [0.0, 0.0],
            [2.0, -1.0],
            [3.0, 2.0],
        ])

        result = self._call_function_under_test(nodes, degree, s_vals)
        self.assertEqual(result.shape, (4, 2))

        for index in six.moves.xrange(4):
            s_val = s_vals[index]
            self.assertEqual(result[index, 0], s_val * (4.0 - s_val))
            self.assertEqual(result[index, 1],
                             2.0 * s_val * (2.0 * s_val - 1.0))
