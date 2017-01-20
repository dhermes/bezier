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


class Test_evaluate(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, x, y):
        from bezier import _implicitization

        return _implicitization.evaluate(nodes, x, y)

    def test_point(self):
        nodes = np.asfortranarray([[1.0, 1.0]])
        with self.assertRaises(ValueError):
            self._call_function_under_test(nodes, 0.0, 0.0)

    def test_linear(self):
        # f(x, y) = -4 x + y + 3
        nodes = np.asfortranarray([
            [1.0, 1.0],
            [2.0, 5.0],
        ])

        result0 = self._call_function_under_test(nodes, 0.0, 0.0)
        self.assertEqual(result0, 3.0)
        result1 = self._call_function_under_test(nodes, 0.0, 1.0)
        self.assertEqual(result1, 4.0)
        result2 = self._call_function_under_test(nodes, 1.0, 0.0)
        self.assertEqual(result2, -1.0)
        result3 = self._call_function_under_test(nodes, 1.0, 1.0)
        self.assertEqual(result3, 0.0)

        # f(x, y) = (-12 x + 8 y - 5) / 32
        nodes = np.asfortranarray([
            [0.0, 0.625],
            [0.25, 1.0],
        ])

        result0 = self._call_function_under_test(nodes, 0.0, 0.0)
        self.assertEqual(result0, -5.0 / 32)
        result1 = self._call_function_under_test(nodes, 0.0, 1.0)
        self.assertEqual(result1, 3.0 / 32)
        result2 = self._call_function_under_test(nodes, 1.0, 0.0)
        self.assertEqual(result2, -17.0 / 32)
        result3 = self._call_function_under_test(nodes, 1.0, 1.0)
        self.assertEqual(result3, -9.0 / 32)

        # f(x, y) = -x
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.0, 1.0],
        ])

        vals = np.linspace(0.0, 1.0, 9)
        for xv in vals:
            for yv in vals:
                result = self._call_function_under_test(nodes, xv, yv)
                self.assertEqual(result, -xv)

    def test_quadratic(self):
        # f(x, y) = x^2 + 4 x - 4 y
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
        ])

        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(nodes, 0.0, 0.0)
