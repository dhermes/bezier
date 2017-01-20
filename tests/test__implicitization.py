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

from tests import utils


class Test_evaluate(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, x_val, y_val):
        from bezier import _implicitization

        return _implicitization.evaluate(nodes, x_val, y_val)

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
        for x_val in vals:
            for y_val in vals:
                result = self._call_function_under_test(nodes, x_val, y_val)
                self.assertEqual(result, -x_val)

    def test_quadratic(self):
        # f(x, y) = x^2 + 4 x - 4 y
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
        ])

        values = [
            self._call_function_under_test(nodes, 0.0, 0.0),
            self._call_function_under_test(nodes, 1.0, 0.0),
            self._call_function_under_test(nodes, 2.0, 0.0),
            self._call_function_under_test(nodes, 0.0, 1.0),
            self._call_function_under_test(nodes, 1.0, 1.0),
            self._call_function_under_test(nodes, 0.0, 2.0),
        ]
        expected = [0.0, 5.0, 12.0, -4.0, 1.0, -8.0]
        self.assertEqual(values, expected)

        # f(x, y) = (x - y)^2 - y
        nodes = np.asfortranarray([
            [0.75, 0.25],
            [-0.25, -0.25],
            [-0.25, 0.25],
        ])
        xy_vals = utils.get_random_nodes(
            shape=(50, 2), seed=7930932, num_bits=8)
        values = []
        expected = []
        for x_val, y_val in xy_vals:
            values.append(self._call_function_under_test(nodes, x_val, y_val))
            expected.append((x_val - y_val) * (x_val - y_val) - y_val)
        self.assertEqual(values, expected)

    def test_cubic(self):
        # f(x, y) = 13824 (x^3 - 24 y^2)
        nodes = np.asfortranarray([
            [6.0, -3.0],
            [-2.0, 3.0],
            [-2.0, -3.0],
            [6.0, 3.0],
        ])

        local_eps = 0.5**26  # sqrt(machine precision)

        xy_vals = utils.get_random_nodes(
            shape=(50, 2), seed=238382, num_bits=8)
        for x_val, y_val in xy_vals:
            result = self._call_function_under_test(nodes, x_val, y_val)
            expected = 13824.0 * (x_val * x_val * x_val - 24.0 * y_val * y_val)
            self.assertAlmostEqual(result, expected, delta=local_eps)

    def test_quartic(self):
        # f(x, y) = -28 x^4 + 56 x^3 - 36 x^2 + 8 x - y
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.25, 2.0],
            [0.5, -2.0],
            [0.75, 2.0],
            [1.0, 0.0],
        ])

        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(nodes, 0.0, 0.0)
