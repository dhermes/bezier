# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest
import unittest.mock

import numpy as np
try:
    import scipy.linalg.lapack as SCIPY_LAPACK
except ImportError:  # pragma: NO COVER
    SCIPY_LAPACK = None
import six  # noqa: I202

from tests import utils as base_utils
from tests.unit import utils


FLOAT64 = np.float64  # pylint: disable=no-member
SPACING = np.spacing  # pylint: disable=no-member
LOCAL_EPS = 0.5**25  # 2 * sqrt(machine precision)


class Test__evaluate3(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, x_val, y_val):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._evaluate3(nodes, x_val, y_val)

    @staticmethod
    def _compute_expected(x_val, y_val):
        return 289.0 * (
            ((17.0 * x_val - 81.0) * x_val + 135.0) * x_val - 27.0 * y_val)

    def test_it(self):
        # f(x, y) = 289(17 x^3 - 81 x^2 + 135 x - 27 y)
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 5.0],
            [2.0, 1.0],
            [3.0, 5.0],
        ])

        xy_vals = utils.get_random_nodes(
            shape=(50, 2), seed=81390, num_bits=8)
        for x_val, y_val in xy_vals:
            result = self._call_function_under_test(nodes, x_val, y_val)
            expected = self._compute_expected(x_val, y_val)
            self.assertAlmostEqual(result, expected, delta=LOCAL_EPS)


class Test_evaluate(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, x_val, y_val):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.evaluate(nodes, x_val, y_val)

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

        xy_vals = utils.get_random_nodes(
            shape=(50, 2), seed=238382, num_bits=8)
        for x_val, y_val in xy_vals:
            result = self._call_function_under_test(nodes, x_val, y_val)
            expected = 13824.0 * (x_val * x_val * x_val - 24.0 * y_val * y_val)
            self.assertAlmostEqual(result, expected, delta=LOCAL_EPS)

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


class Test_eval_intersection_polynomial(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2, t):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.eval_intersection_polynomial(
            nodes1, nodes2, t)

    def test_degrees_1_1(self):
        # f1(x, y) = (8 y - 3) / 8
        nodes1 = np.asfortranarray([
            [0.0, 0.375],
            [1.0, 0.375],
        ])
        # x2(t), y2(t) = 1 / 2, 3 s / 4
        nodes2 = np.asfortranarray([
            [0.5, 0.0],
            [0.5, 0.75],
        ])

        values = [
            self._call_function_under_test(nodes1, nodes2, 0.0),
            self._call_function_under_test(nodes1, nodes2, 0.5),
            self._call_function_under_test(nodes1, nodes2, 1.0),
        ]
        # f1(x2(t), y2(t)) = 3 (2 t - 1) / 8
        expected = [-0.375, 0.0, 0.375]
        self.assertEqual(values, expected)

    def test_degrees_1_2(self):
        # f1(x, y) = 2 (4 x + 3 y - 24)
        nodes1 = np.asfortranarray([
            [0.0, 8.0],
            [6.0, 0.0],
        ])
        # x2(t), y2(t) = 9 t, 18 t (1 - t)
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [4.5, 9.0],
            [9.0, 0.0],
        ])

        values = [
            self._call_function_under_test(nodes1, nodes2, 0.0),
            self._call_function_under_test(nodes1, nodes2, 0.25),
            self._call_function_under_test(nodes1, nodes2, 0.5),
            self._call_function_under_test(nodes1, nodes2, 0.75),
            self._call_function_under_test(nodes1, nodes2, 1.0),
        ]
        # f1(x2(t), y2(t)) = 12 (4 - 3 t) (3 t - 1)
        expected = [-48.0, -9.75, 15.0, 26.25, 24.0]
        self.assertEqual(values, expected)


class Test__to_power_basis11(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._to_power_basis11(nodes1, nodes2)

    def test_it(self):
        # f1(x, y) = -(12 x - 8 y + 5) / 32
        nodes1 = np.asfortranarray([
            [0.0, 0.625],
            [0.25, 1.0],
        ])
        # x2(t), y2(t) = (2 - 3 t) / 4, (3 t + 2) / 4
        nodes2 = np.asfortranarray([
            [0.5, 0.5],
            [-0.25, 1.25],
        ])
        # f1(x2(t), y2(t)) = (15 t - 7) / 32
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray([-7.0, 15.0]) / 32.0
        self.assertEqual(result, expected)


class Test__to_power_basis12(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._to_power_basis12(nodes1, nodes2)

    def test_it(self):
        # f1(x, y) = (2 y - 1) / 2
        nodes1 = np.asfortranarray([
            [0.0, 0.5],
            [1.0, 0.5],
        ])
        # x2(t), y2(t) = t, 2 t (1 - t)
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ])
        # f1(x2(t), y2(t)) = -(2 t - 1)^2 / 2
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray([-0.5, 2.0, -2.0])
        self.assertEqual(result, expected)


class Test__to_power_basis13(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._to_power_basis13(nodes1, nodes2)

    def test_it(self):
        # f1(x, y) = -(152 x + 112 y - 967) / 64
        nodes1 = np.asfortranarray([
            [2.5625, 5.15625],
            [0.8125, 7.53125],
        ])
        # x2(t), y2(t) = 3 (14 t + 1) / 8, 18 t^3 - 27 t^2 + 3 t + 7
        nodes2 = np.asfortranarray([
            [0.375, 7.0],
            [2.125, 8.0],
            [3.875, 0.0],
            [5.625, 1.0],
        ])
        # f1(x2(t), y2(t)) = -63 (t - 1) (4 t - 1)^2 / 32
        result = self._call_function_under_test(nodes1, nodes2)
        # The 1-3 method avoids a division by 3.0
        expected = 3.0 * (-63.0 / 32.0) * np.asfortranarray(
            [-1.0, 9.0, -24.0, 16.0])
        self.assertEqual(result, expected)


class Test__to_power_basis_degree4(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._to_power_basis_degree4(nodes1, nodes2)

    def test_degrees_2_2(self):
        # f1(x, y) = (x^2 - 4 x y + 4 y^2 - y) / 16
        nodes1 = np.asfortranarray([
            [0.375, 0.0625],
            [-0.125, -0.0625],
            [-0.125, 0.0625],
        ])
        # x2(t), y2(t) = (2 t - 3) (2 t - 1) / 4, (2 t - 1)^2 / 4
        nodes2 = np.asfortranarray([
            [0.75, 0.25],
            [-0.25, -0.25],
            [-0.25, 0.25],
        ])
        # f1(x2(t), y2(t)) = (2 t - 1)^3 (2 t + 3) / 256
        result = self._call_function_under_test(nodes1, nodes2)
        # The function avoids a division by 3.0
        expected = (3.0 / 256.0) * np.asfortranarray(
            [-3.0, 16.0, -24.0, 0.0, 16.0])
        self.assertEqual(result, expected)

    def test_degrees_1_4(self):
        # f1(x, y) = 4 (y - 3)
        nodes1 = np.asfortranarray([
            [0.0, 3.0],
            [4.0, 3.0],
        ])
        # x2(t), y2(t) = 4 s, 2 s (s - 1) (5 s^2 - 5 s - 2)
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
            [3.0, 1.0],
            [4.0, 0.0],
        ])
        # f1(x2(t), y2(t)) = 4 (10 t^4 - 20 t^3 + 6 t^2 + 4 t - 3)
        result = self._call_function_under_test(nodes1, nodes2)
        # The function avoids a division by 3.0
        expected = 3.0 * np.asfortranarray([-12.0, 16.0, 24.0, -80.0, 40.0])
        self.assertEqual(result, expected)


class Test__to_power_basis23(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._to_power_basis23(nodes1, nodes2)

    def test_it(self):
        # f1(x, y) = 4 (4 x^2 - 12 x - 4 y + 11)
        nodes1 = np.asfortranarray([
            [0.5, 1.5],
            [1.5, -0.5],
            [2.5, 1.5],
        ])
        # x2(t), y2(t) = 3 t, t (4 t^2 - 6 t + 3)
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 0.0],
            [3.0, 1.0],
        ])
        # f1(x2(t), y2(t)) = 4 (2 s - 1)^2 (4 s - 11)
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray([44, -192, 240, -64, 0.0, 0.0, 0.0])
        self.assertLess(np.abs(result - expected).max(), LOCAL_EPS)


class Test__to_power_basis_degree8(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._to_power_basis_degree8(nodes1, nodes2)

    def test_degrees_2_4(self):
        # f1(x, y) = 2 (9 x - 2 y^2 - 6 y)
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.5],
            [4.0, 3.0],
        ])
        # x2(t), y2(t) = 4 t, -t (7 t^3 - 6 t - 4)
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
            [3.0, 6.0],
            [4.0, 3.0],
        ])
        # f1(x2(t), y2(t)) = (
        #     4 t (t - 1) (49 t^6 + 49 t^5 - 35 t^4 -
        #                  91 t^3 - 76 t^2 - 28 t + 6))
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray(
            [0.0, -24.0, 136.0, 192.0, 60.0, -224.0, -336.0, 0.0, 196.0])
        self.assertLess(np.abs(result - expected).max(), LOCAL_EPS)


class Test__to_power_basis33(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._to_power_basis33(nodes1, nodes2)

    def test_it(self):
        # f1(x, y) = x^3 - 9 x^2 + 27 x - 27 y
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 1.0],
            [3.0, 1.0],
        ])
        # x2(t), y2(t) = (1 - t) (t^2 + t + 1), 3 (1 - t)
        nodes2 = np.asfortranarray([
            [1.0, 3.0],
            [1.0, 2.0],
            [1.0, 1.0],
            [0.0, 0.0],
        ])
        # f1(x2(t), y2(t)) = (
        #     -(s - 1)^2 (s + 2) (s^6 + 3 s^4 + 4 s^3 + 9 s^2 + 6 s + 31)
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray(
            [-62, 81, 0, -12, 0, 0, -6, 0, 0, -1], dtype=FLOAT64)
        self.assertLess(np.abs(result - expected).max(), LOCAL_EPS)


class Test_to_power_basis(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.to_power_basis(nodes1, nodes2)

    def test_degrees_1_1(self):
        # f1(x, y) = -x
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [0.0, 1.0],
        ])
        # x2(t), y2(t) = 1, t
        nodes2 = np.asfortranarray([
            [1.0, 0.0],
            [1.0, 1.0],
        ])
        # f1(x2(t), y2(t)) = -1
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray([-1.0, 0.0])
        self.assertEqual(result, expected)

    def test_degrees_1_2(self):
        # f1(x, y) = 2 (4 x + 3 y - 24)
        nodes1 = np.asfortranarray([
            [0.0, 8.0],
            [6.0, 0.0],
        ])
        # x2(t), y2(t) = 9 t, 18 t (1 - t)
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [4.5, 9.0],
            [9.0, 0.0],
        ])
        # f1(x2(t), y2(t)) = 12 (4 - 3 t) (3 t - 1)
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray([-48.0, 180.0, -108.0])
        self.assertEqual(result, expected)

    def test_degrees_1_3(self):
        # f1(x, y) = -(2 x - y - 1) / 2
        nodes1 = np.asfortranarray([
            [0.5, 0.0],
            [1.0, 1.0],
        ])
        # x2(t), y2(t) = -t (2 t^2 - 3 t - 3) / 4, -(3 t^3 - 3 t - 1) / 2
        nodes2 = np.asfortranarray([
            [0.0, 0.5],
            [0.25, 1.0],
            [0.75, 1.5],
            [1.0, 0.5],
        ])
        # f1(x2(t), y2(t)) = -(t^3 + 3 t^2 - 3) / 4
        result = self._call_function_under_test(nodes1, nodes2)
        # The 1-3 method avoids a division by 3.0
        expected = 3.0 * np.asfortranarray([0.75, 0.0, -0.75, -0.25])
        self.assertEqual(result, expected)

    def test_degrees_1_4(self):
        # f1(x, y) = 4 y - 3 x
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [4.0, 3.0],
        ])
        # x2(t), y2(t) = 4 t, -t (7 t^3 - 6 t - 4)
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
            [3.0, 6.0],
            [4.0, 3.0],
        ])
        # f1(x2(t), y2(t)) = 4 t (1 - t) (7 t^2 + 7 t + 1)
        result = self._call_function_under_test(nodes1, nodes2)
        # The function avoids a division by 3.0
        expected = 3.0 * np.asfortranarray([0.0, 4.0, 24.0, 0.0, -28.0])
        self.assertEqual(result, expected)

    def test_degrees_2_2(self):
        # f1(x, y) = (x^2 - 2 x y + 4 x + y^2 + 4 y - 5) / 4
        nodes1 = np.asfortranarray([
            [1.0, 0.0],
            [0.75, 0.75],
            [0.0, 1.0],
        ])
        # x2(t), y2(t) = (t - 1) (5 t - 8) / 8, -(t - 1) (7 t + 8) / 8
        nodes2 = np.asfortranarray([
            [1.0, 1.0],
            [0.1875, 0.9375],
            [0.0, 0.0],
        ])
        # f1(x2(t), y2(t)) = (9 t^4 - 18 t^3 + 5 t^2 - 28 t + 12) / 16
        result = self._call_function_under_test(nodes1, nodes2)
        # The function avoids a division by 3.0
        expected = (3.0 / 16.0) * np.asfortranarray(
            [12.0, -28.0, 5.0, -18.0, 9.0])
        self.assertEqual(result, expected)

    def test_degrees_2_3(self):
        # f1(x, y) = 81 (2 x^2 - 2 x - y + 1) / 128
        nodes1 = np.asfortranarray([
            [0.25, 0.625],
            [0.625, 0.25],
            [1.0, 1.0],
        ])
        # x2(t), y2(t) = -t (2 t^2 - 3 t - 3) / 4, -(3 t^3 - 3 t - 1) / 2
        nodes2 = np.asfortranarray([
            [0.0, 0.5],
            [0.25, 1.0],
            [0.75, 1.5],
            [1.0, 0.5],
        ])
        # f1(x2(t), y2(t)) = (
        #     81 (t + 1) (4 t^5 - 16 t^4 + 13 t^3 + 25 t^2 - 28 t + 4) / 1024)
        result = self._call_function_under_test(nodes1, nodes2)
        expected = (81.0 / 1024.0) * np.asfortranarray(
            [4, -24, -3, 38, -3, -12, 4])
        self.assertTrue(
            np.allclose(result, expected, atol=0.0, rtol=LOCAL_EPS))

    def test_degrees_2_4(self):
        # f1(x, y) = 2*(9*x - 2*y**2 - 6*y)
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.5],
            [4.0, 3.0],
        ])
        # x2(t), y2(t) = 4*t, 2*t*(t - 1)*(5*t**2 - 5*t - 2)]])
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
            [3.0, 1.0],
            [4.0, 0.0],
        ])
        # f1(x2(t), y2(t)) = (
        #     8 t (50 t^7 - 200 t^6 + 260 t^5 - 80 t^4 -
        #          47 t^3 - 6 t^2 + 17 t - 3))
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray(
            [0.0, -24.0, 136.0, -48.0, -376.0, -640.0, 2080.0, -1600.0, 400.0])
        self.assertLess(np.abs(result - expected).max(), LOCAL_EPS)

    def test_degrees_3_3(self):
        # f1(x, y) = -(13824 x^3 + 3456 x^2 y - 55296 x^2 +
        #              288 x y^2 - 2088 x y + 39816 x + 8 y^3 +
        #              129 y^2 + 6846 y - 6983) / 512
        nodes1 = np.asfortranarray([
            [0.0, 1.0],
            [0.375, -1.0],
            [0.625, 0.0],
            [1.0, 1.0],
        ])
        # x2(t), y2(t) = -t (2 t^2 - 3 t - 3) / 4, -(3 t^3 - 3 t - 1) / 2
        nodes2 = np.asfortranarray([
            [0.0, 0.5],
            [0.25, 1.0],
            [0.75, 1.5],
            [1.0, 0.5],
        ])
        # f1(x2(t), y2(t)) = (
        #     13500 t^9 - 48600 t^8 + 1620 t^7 + 170451 t^6 - 171072 t^5 -
        #     146394 t^4 + 331686 t^3 + 10827 t^2 - 158418 t + 14107) / 2048
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray([
            14107, -158418, 10827, 331686, -146394,
            -171072, 170451, 1620, -48600, 13500]) / 2048.0
        self.assertTrue(
            np.allclose(result, expected, atol=0.0, rtol=LOCAL_EPS))

    def test_unsupported(self):
        nodes_yes1 = np.zeros((2, 2), order='F')
        nodes_yes2 = np.zeros((3, 2), order='F')
        nodes_yes3 = np.zeros((4, 2), order='F')
        nodes_no = np.zeros((6, 2), order='F')

        # Just make sure we fall through **all** of the implicit
        # ``else`` branches.
        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(nodes_yes1, nodes_no)

        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(nodes_yes2, nodes_no)

        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(nodes_yes3, nodes_no)

        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(nodes_no, nodes_no)


class Test_polynomial_norm(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(coeffs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.polynomial_norm(coeffs)

    def test_it(self):
        coeffs = np.asfortranarray([2.0, 1.0, 3.0])
        result = self._call_function_under_test(coeffs)
        expected = np.sqrt(409.0 / 30.0)
        self.assertAlmostEqual(result, expected, delta=LOCAL_EPS)


class Test_roots_in_unit_interval(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(coeffs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.roots_in_unit_interval(coeffs)

    def test_it(self):
        # P(t) = (t^2 + 1) (2 t - 1) (3 t - 1) (3 t - 2)
        coeffs = np.asfortranarray([-2.0, 13.0, -29.0, 31.0, -27.0, 18.0])
        all_roots = self._call_function_under_test(coeffs)
        all_roots = np.sort(all_roots)
        self.assertEqual(all_roots.shape, (3,))
        for index in (0, 1, 2):
            expected = (index + 2.0) / 6.0
            self.assertAlmostEqual(
                all_roots[index], expected, delta=LOCAL_EPS)


class Test__strip_leading_zeros(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(coeffs, **kwargs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._strip_leading_zeros(coeffs, **kwargs)

    def test_default_threshold(self):
        coeffs = np.asfortranarray([0.0, 1.0, 1.5])
        result = self._call_function_under_test(coeffs)
        self.assertIs(result, coeffs)

    def test_custom_threshold(self):
        coeffs = np.asfortranarray([2.0, 0.0, 0.0, 0.5**10])
        result = self._call_function_under_test(coeffs)
        self.assertIs(result, coeffs)

        result = self._call_function_under_test(coeffs, threshold=0.5**9)
        self.assertEqual(result, coeffs[:1])


class Test__check_non_simple(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(coeffs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._check_non_simple(coeffs)

    def test_extra_zeros(self):
        coeffs = np.asfortranarray([2.0, -3.0, 1.0, 0.0])
        # Just make sure no exception was thrown.
        self.assertIsNone(self._call_function_under_test(coeffs))

    def test_line(self):
        coeffs = np.asfortranarray([4.0, 1.0])
        # Just make sure no exception was thrown.
        self.assertIsNone(self._call_function_under_test(coeffs))

    def test_double_root(self):
        # f(t) = (t - 2)^2
        coeffs = np.asfortranarray([4.0, -4.0, 1.0])
        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(coeffs)

        # f(t) = (t + 2)^2 (3 t + 2) (4 t + 19)
        coeffs = np.asfortranarray([152.0, 412.0, 346.0, 113.0, 12.0])
        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(coeffs)

    def test_scale_invariant(self):
        # f(t) = (t - 1) (t - 2) (t - 3) (t - 4)
        coeffs = np.asfortranarray([24.0, -50.0, 35.0, -10.0, 1.0])
        # Make sure no exception was thrown.
        self.assertIsNone(self._call_function_under_test(coeffs))

        for exponent in (-30, -20, -10, 10, 20, 30):
            new_coeffs = 0.5**exponent * coeffs
            self.assertIsNone(self._call_function_under_test(new_coeffs))


class Test__resolve_and_add(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(
            nodes1, s_val, final_s, nodes2, t_val, final_t):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._resolve_and_add(
            nodes1, s_val, final_s, nodes2, t_val, final_t)

    def _helper(self, s, t):
        nodes1 = unittest.mock.sentinel.nodes1
        nodes2 = unittest.mock.sentinel.nodes2
        final_s = []
        final_t = []

        patch = unittest.mock.patch(
            'bezier._intersection_helpers.newton_refine',
            return_value=(s, t))
        with patch as mocked:
            self._call_function_under_test(
                nodes1, s, final_s, nodes2, t, final_t)
            mocked.assert_called_once_with(s, nodes1, t, nodes2)

        return final_s, final_t

    def test_unchanged(self):
        s = 0.5
        t = 0.25
        final_s, final_t = self._helper(s, t)
        self.assertEqual(final_s, [s])
        self.assertEqual(final_t, [t])

    def test_to_zero(self):
        s = -0.5**60
        t = 0.5
        final_s, final_t = self._helper(s, t)
        self.assertEqual(final_s, [0.0])
        self.assertEqual(final_t, [t])

    def test_still_negative(self):
        s = 0.125
        t = -0.5**20
        final_s, final_t = self._helper(s, t)
        self.assertEqual(final_s, [])
        self.assertEqual(final_t, [])


class Test_intersect_curves(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.intersect_curves(nodes1, nodes2)

    def test_degrees_1_1(self):
        # f1(x, y) = (8 y - 3) / 8
        nodes1 = np.asfortranarray([
            [0.0, 0.375],
            [1.0, 0.375],
        ])
        # x2(t), y2(t) = 1 / 2, 3 s / 4
        nodes2 = np.asfortranarray([
            [0.5, 0.0],
            [0.5, 0.75],
        ])
        # f1(x2(t), y2(t)) = 3 (2 t - 1) / 8
        result = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray([
            [0.5, 0.5],
        ])
        self.assertEqual(result, expected)

    def _degrees_1_2_helper(self, swapped=False):
        # f1(x, y) = 2 (4 x + 3 y - 24)
        nodes1 = np.asfortranarray([
            [0.0, 8.0],
            [6.0, 0.0],
        ])
        # x2(t), y2(t) = 9 t, 18 t (1 - t)
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [4.5, 9.0],
            [9.0, 0.0],
        ])

        # f1(x2(t), y2(t)) = 12 (4 - 3 t) (3 t - 1)
        if swapped:
            col1, col2 = 1, 0
            result = self._call_function_under_test(nodes2, nodes1)
        else:
            col1, col2 = 0, 1
            result = self._call_function_under_test(nodes1, nodes2)

        self.assertEqual(result.shape, (1, 2))
        self.assertAlmostEqual(result[0, col1], 0.5, delta=LOCAL_EPS)
        self.assertAlmostEqual(result[0, col2], 1.0 / 3.0, delta=LOCAL_EPS)

    def test_degrees_1_2(self):
        self._degrees_1_2_helper()

    def test_degrees_2_1(self):
        self._degrees_1_2_helper(swapped=True)

    def test_t_val_without_s_val(self):
        # f1(x, y) = 81 (2 x^2 - 2 x - y + 1) / 128
        nodes1 = np.asfortranarray([
            [0.25, 0.625],
            [0.625, 0.25],
            [1.0, 1.0],
        ])
        # x2(t), y2(t) = -t (2 t^2 - 3 t - 3) / 4, -(3 t^3 - 3 t - 1) / 2
        nodes2 = np.asfortranarray([
            [0.0, 0.5],
            [0.25, 1.0],
            [0.75, 1.5],
            [1.0, 0.5],
        ])

        # f1(x2(t), y2(t)) = (
        #     81 (t + 1) (4 t^5 - 16 t^4 + 13 t^3 + 25 t^2 - 28 t + 4) / 1024)
        # NOTE: This polynomial has two roots inside the unit interval
        #       t1 = 0.17072782629577715 and t2 = 0.8734528541508780397 but
        #       they correspond to s-values s1 = -0.1367750984247189222649977
        #       and s2 = 0.8587897065534014787757 so only s2, t2 is returned.
        result = self._call_function_under_test(nodes1, nodes2)
        self.assertEqual(result.shape, (1, 2))
        # 486 s^5 - 3726 s^4 + 13905 s^3 - 18405 s^2 + 6213 s + 1231
        self.assertAlmostEqual(
            result[0, 0], float.fromhex('0x1.b7b348cf939b9p-1'),
            delta=LOCAL_EPS)
        # 4 t^5 - 16 t^4 + 13 t^3 + 25 t^2 - 28 t + 4
        self.assertAlmostEqual(
            result[0, 1], float.fromhex('0x1.bf3536665a0cdp-1'),
            delta=LOCAL_EPS)

    def test_coincident(self):
        # f1(x, y) = -4 (4 x^2 + 4 x y - 16 x + y^2 + 8 y)
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 2.0],
            [4.0, 0.0],
        ])
        # x2(t), y2(t) = 3 (t + 1) (3 t + 7) / 8, -3 (t + 1) (3 t - 1) / 4
        nodes2 = np.asfortranarray([
            [2.625, 0.75],
            [4.5, 0.0],
            [7.5, -3.0],
        ])

        # f1(x2(t), y2(t)) = 0
        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(nodes1, nodes2)


class Test_normalize_polynomial(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(coeffs, **kwargs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.normalize_polynomial(coeffs, **kwargs)

    def test_nonzero(self):
        coeffs = np.asfortranarray([2.0])
        result = self._call_function_under_test(coeffs)
        self.assertIs(result, coeffs)
        # Check changed in place.
        expected = np.asfortranarray([1.0])
        self.assertEqual(result, expected)

        coeffs = np.asfortranarray([-2.0, 6.0])
        result = self._call_function_under_test(coeffs)
        self.assertIs(result, coeffs)
        # Check changed in place.
        expected = np.asfortranarray([-1.0, 3.0])
        self.assertEqual(result, expected)

    def test_actual_zero(self):
        for num_vals in (1, 2, 3, 4, 5):
            coeffs = np.zeros((num_vals,), order='F')
            result = self._call_function_under_test(coeffs)
            self.assertIsNot(result, coeffs)
            self.assertEqual(result, coeffs)

    def test_almost_zero(self):
        shape = (4,)
        coeffs = 0.5**42 * np.random.random(shape)
        result = self._call_function_under_test(coeffs)
        self.assertIsNot(result, coeffs)
        self.assertEqual(result, np.zeros(shape, order='F'))

        coeffs = 0.5**10 * np.random.random(shape)
        result = self._call_function_under_test(coeffs, threshold=0.5**8)
        self.assertIsNot(result, coeffs)
        self.assertEqual(result, np.zeros(shape, order='F'))


class Test__get_sigma_coeffs(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(coeffs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._get_sigma_coeffs(coeffs)

    def test_all_zero(self):
        for num_zeros in (1, 2, 3, 4):
            coeffs = np.zeros((num_zeros,), order='F')
            result = self._call_function_under_test(coeffs)
            sigma_coeffs, degree, effective_degree = result
            self.assertIsNone(sigma_coeffs)
            self.assertEqual(degree, 0)
            self.assertEqual(effective_degree, 0)

    def test_linear(self):
        # s
        coeffs = np.asfortranarray([0.0, 1.0])
        result = self._call_function_under_test(coeffs)
        sigma_coeffs, degree, effective_degree = result
        self.assertEqual(sigma_coeffs, np.asfortranarray([0.0]))
        self.assertEqual(degree, 1)
        self.assertEqual(effective_degree, 1)

    def test_linear_drop_degree(self):
        # 4 (1 - s)
        coeffs = np.asfortranarray([4.0, 0.0])
        result = self._call_function_under_test(coeffs)
        sigma_coeffs, degree, effective_degree = result
        self.assertIsNone(sigma_coeffs)
        self.assertEqual(degree, 1)
        self.assertEqual(effective_degree, 0)

    def test_quadratic(self):
        # 2 s (s + 1)
        coeffs = np.asfortranarray([0.0, 1.0, 4.0])
        result = self._call_function_under_test(coeffs)
        sigma_coeffs, degree, effective_degree = result
        expected = np.asfortranarray([0.0, 0.5])
        self.assertEqual(sigma_coeffs, expected)
        self.assertEqual(degree, 2)
        self.assertEqual(effective_degree, 2)

    def test_quadratic_drop_degree(self):
        # 2 (s - 2) (s - 1)
        coeffs = np.asfortranarray([4.0, 1.0, 0.0])
        result = self._call_function_under_test(coeffs)
        sigma_coeffs, degree, effective_degree = result
        self.assertEqual(sigma_coeffs, np.asfortranarray([2.0]))
        self.assertEqual(degree, 2)
        self.assertEqual(effective_degree, 1)

    def test_cubic(self):
        # -(s - 17) (s - 5) (s - 2)
        coeffs = np.asfortranarray([170.0, 127.0, 92.0, 64.0])
        result = self._call_function_under_test(coeffs)
        sigma_coeffs, degree, effective_degree = result
        expected = np.asfortranarray([2.65625, 5.953125, 4.3125])
        self.assertEqual(sigma_coeffs, expected)
        self.assertEqual(degree, 3)
        self.assertEqual(effective_degree, 3)

    def test_cubic_drop_degree(self):
        # 3 (1 - s)^2 (3 s + 1)
        coeffs = np.asfortranarray([3.0, 4.0, 0.0, 0.0])
        result = self._call_function_under_test(coeffs)
        sigma_coeffs, degree, effective_degree = result
        self.assertEqual(sigma_coeffs, np.asfortranarray([0.25]))
        self.assertEqual(degree, 3)
        self.assertEqual(effective_degree, 1)

    def test_quartic(self):
        # 4 s^3 (5 - s)
        coeffs = np.asfortranarray([0.0, 0.0, 0.0, 5.0, 16.0])
        result = self._call_function_under_test(coeffs)
        sigma_coeffs, degree, effective_degree = result
        expected = np.asfortranarray([0.0, 0.0, 0.0, 1.25])
        self.assertEqual(sigma_coeffs, expected)
        self.assertEqual(degree, 4)
        self.assertEqual(effective_degree, 4)


class Test_bernstein_companion(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(coeffs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.bernstein_companion(coeffs)

    def test_all_zero(self):
        for num_zeros in (1, 2, 3, 4):
            coeffs = np.zeros((num_zeros,), order='F')
            result = self._call_function_under_test(coeffs)
            companion, degree, effective_degree = result
            self.assertEqual(companion.shape, (0, 0))
            self.assertEqual(degree, 0)
            self.assertEqual(effective_degree, 0)

    def test_quadratic(self):
        # 2 s (s + 1)
        coeffs = np.asfortranarray([0.0, 1.0, 4.0])
        companion, degree, effective_degree = self._call_function_under_test(
            coeffs)
        expected = np.asfortranarray([
            [-0.5, 0.0],
            [1.0, 0.0],
        ])
        self.assertEqual(companion, expected)
        self.assertEqual(degree, 2)
        self.assertEqual(effective_degree, 2)


class Test_bezier_roots(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(coeffs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.bezier_roots(coeffs)

    def test_all_zero(self):
        for num_zeros in (1, 2, 3, 4):
            coeffs = np.zeros((num_zeros,), order='F')
            roots = self._call_function_under_test(coeffs)
            self.assertEqual(roots.shape, (0,))

    def test_linear(self):
        # s
        coeffs = np.asfortranarray([0.0, 1.0])
        roots = self._call_function_under_test(coeffs)
        self.assertEqual(roots, np.asfortranarray([0.0]))

    def test_linear_drop_degree(self):
        # 4 (1 - s)
        coeffs = np.asfortranarray([4.0, 0.0])
        roots = self._call_function_under_test(coeffs)
        self.assertEqual(roots, np.asfortranarray([1.0]))

    def test_linear_as_elevated(self):
        # 2
        coeffs = np.asfortranarray([2.0, 2.0])
        roots = self._call_function_under_test(coeffs)
        self.assertEqual(roots.shape, (0,))

    def test_quadratic(self):
        # 2 s (s + 1)
        coeffs = np.asfortranarray([0.0, 1.0, 4.0])
        roots = self._call_function_under_test(coeffs)

        roots = np.sort(roots)
        expected = np.asfortranarray([-1.0, 0.0])
        self.assertEqual(roots, expected)

    def test_quadratic_complex(self):
        # s^2 + 1
        coeffs = np.asfortranarray([1.0, 1.0, 2.0])
        roots = self._call_function_under_test(coeffs)

        roots = np.sort(roots)
        expected = np.asfortranarray([-1.0j, 1.0j])
        ulp_errs = np.abs((roots - expected) / SPACING(expected.imag))
        self.assertEqual(ulp_errs.shape, (2,))
        self.assertLess(ulp_errs[0], 2)
        self.assertEqual(ulp_errs[0], ulp_errs[1])

    def test_quadratic_drop_degree(self):
        # 2 (s - 2) (s - 1)
        coeffs = np.asfortranarray([4.0, 1.0, 0.0])
        roots = self._call_function_under_test(coeffs)
        self.assertEqual(roots, np.asfortranarray([2.0, 1.0]))

    def test_quadratic_as_elevated(self):
        # 2 (s + 1)
        coeffs = np.asfortranarray([2.0, 3.0, 4.0])
        roots = self._call_function_under_test(coeffs)
        self.assertEqual(roots, np.asfortranarray([-1.0]))

    def test_cubic(self):
        # -(s - 17) (s - 5) (s - 2)
        coeffs = np.asfortranarray([170.0, 127.0, 92.0, 64.0])
        roots = self._call_function_under_test(coeffs)

        roots = np.sort(roots)
        expected = np.asfortranarray([2.0, 5.0, 17.0])
        ulp_errs = np.abs((roots - expected) / SPACING(expected))
        self.assertEqual(ulp_errs.shape, (3,))
        self.assertLess(ulp_errs[0], 16)
        self.assertLess(ulp_errs[1], 512)
        self.assertLess(ulp_errs[2], 1024)

    def test_cubic_drop_degree(self):
        # 6 (1 - s)^2 (s + 1)
        coeffs = np.asfortranarray([6.0, 4.0, 0.0, 0.0])
        roots = self._call_function_under_test(coeffs)
        self.assertEqual(roots, np.asfortranarray([-1.0, 1.0, 1.0]))

    def test_cubic_as_elevated(self):
        # 3 (s - 5) (4 s - 3)
        coeffs = np.asfortranarray([45.0, 22.0, 3.0, -12.0])
        roots = self._call_function_under_test(coeffs)

        roots = np.sort(roots)
        expected = np.asfortranarray([0.75, 5.0])
        ulp_errs = np.abs((roots - expected) / SPACING(expected))
        self.assertEqual(ulp_errs.shape, (2,))
        self.assertEqual(ulp_errs[0], 0.0)
        self.assertLess(ulp_errs[1], 64)

    def test_quartic(self):
        # 4 s^3 (5 - s)
        coeffs = np.asfortranarray([0.0, 0.0, 0.0, 5.0, 16.0])
        roots = self._call_function_under_test(coeffs)

        roots = np.sort(roots)
        expected = np.asfortranarray([0.0, 0.0, 0.0, 5.0])
        self.assertEqual(roots, expected)

    def test_quartic_as_twice_elevated(self):
        # 6 (s - 3) (s + 7)
        coeffs = np.asfortranarray([-126.0, -120.0, -113.0, -105.0, -96.0])
        roots = self._call_function_under_test(coeffs)

        roots = np.sort(roots)
        expected = np.asfortranarray([-7.0, 3.0])
        ulp_errs = np.abs((roots - expected) / SPACING(expected))
        self.assertEqual(ulp_errs.shape, (2,))
        self.assertLess(ulp_errs[0], 16384)
        self.assertLess(ulp_errs[1], 256)


class Test_lu_companion(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(top_row, value):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.lu_companion(top_row, value)

    def test_linear(self):
        # t + 3
        top_row = np.asfortranarray([-3.0])
        lu_mat, one_norm = self._call_function_under_test(top_row, -2.0)
        expected = np.asfortranarray([
            [-1.0],
        ])
        self.assertEqual(expected, lu_mat)
        self.assertEqual(one_norm, 1.0)

    def _check_lu(self, lu_mat, expected_a):
        from bezier import _helpers

        rows, cols = lu_mat.shape
        self.assertEqual(rows, cols)

        l_mat = np.asfortranarray(np.tril(lu_mat, -1) + _helpers.eye(rows))
        u_mat = np.triu(lu_mat)
        a_mat = _helpers.matrix_product(l_mat, u_mat)
        self.assertEqual(a_mat, expected_a)

    def test_quadratic(self):
        # t^2 + 5 t - 4
        top_row = np.asfortranarray([-5.0, 4.0])
        value = 2.0
        lu_mat, one_norm = self._call_function_under_test(top_row, value)
        expected = np.asfortranarray([
            [1.0, -value],
            [-7.0, -10.0],
        ])
        self.assertEqual(expected, lu_mat)
        self.assertEqual(one_norm, 8.0)

        expected_a = np.asfortranarray([
            [1.0, -value],
            [-7.0, 4.0],
        ])
        self._check_lu(lu_mat, expected_a)

    def test_cubic(self):
        # t^3 + 3 t^2 - t + 2
        top_row = np.asfortranarray([-3.0, 1.0, -2.0])
        value = 3.0
        lu_mat, one_norm = self._call_function_under_test(top_row, value)
        expected = np.asfortranarray([
            [1.0, -value, 0.0],
            [0.0, 1.0, -value],
            [-6.0, -17.0, -53.0],
        ])

        self.assertEqual(expected, lu_mat)
        self.assertEqual(one_norm, 7.0)

        expected_a = np.asfortranarray([
            [1.0, -value, 0.0],
            [0.0, 1.0, -value],
            [-6.0, 1.0, -2.0],
        ])
        self._check_lu(lu_mat, expected_a)

    def test_quartic(self):
        # t^4 + t^3 + t^2 - 5 t - 2
        top_row = np.asfortranarray([-1.0, -1.0, 5.0, 2.0])
        value = 4.0
        lu_mat, one_norm = self._call_function_under_test(top_row, value)
        expected = np.asfortranarray([
            [1.0, -value, 0.0, 0.0],
            [0.0, 1.0, -value, 0.0],
            [0.0, 0.0, 1.0, -value],
            [-5.0, -21.0, -79.0, -314.0],
        ])
        self.assertEqual(expected, lu_mat)
        self.assertEqual(one_norm, 10.0)

        expected_a = np.asfortranarray([
            [1.0, -value, 0.0, 0.0],
            [0.0, 1.0, -value, 0.0],
            [0.0, 0.0, 1.0, -value],
            [-5.0, -1.0, 5.0, 2.0],
        ])
        self._check_lu(lu_mat, expected_a)


class Test__reciprocal_condition_number(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(lu_mat, one_norm):
        from bezier import _algebraic_intersection

        return _algebraic_intersection._reciprocal_condition_number(
            lu_mat, one_norm)

    @unittest.mock.patch(
        'bezier._algebraic_intersection._scipy_lapack', new=None)
    def test_without_scipy(self):
        lu_mat = np.zeros((2, 2), order='F')
        one_norm = 0.0
        with self.assertRaises(OSError):
            self._call_function_under_test(lu_mat, one_norm)

    @unittest.mock.patch('bezier._algebraic_intersection._scipy_lapack')
    def test_dgecon_failure(self, _scipy_lapack):
        rcond = 0.5
        info = -1
        _scipy_lapack.dgecon.return_value = rcond, info

        one_norm = 1.0
        with self.assertRaises(RuntimeError):
            self._call_function_under_test(
                unittest.mock.sentinel.lu_mat, one_norm)

        _scipy_lapack.dgecon.assert_called_once_with(
            unittest.mock.sentinel.lu_mat, one_norm)

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_singular(self):
        # A = [[ 2.   3. ]
        #      [ 1.   1.5]]
        one_norm = 4.5
        lu_mat = np.asfortranarray([
            [2.0, 3.0],
            [0.5, 0.0],
        ])
        rcond = self._call_function_under_test(lu_mat, one_norm)
        self.assertEqual(rcond, 0.0)

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_invertible(self):
        # A = [[  4.,  10.,   0.],
        #      [  0.,  16.,  32.],
        #      [  2.,  -3.,   0.]]
        one_norm = 32.0
        lu_mat = np.asfortranarray([
            [4.0, 10.0, 0.0],
            [0.0, 16.0, 32.0],
            [0.5, -0.5, 16.0],
        ])
        rcond = self._call_function_under_test(lu_mat, one_norm)
        self.assertEqual(rcond, 0.0625)


class Test_bezier_value_check(utils.NumPyTestCase):

    CUBIC_COEFFS = np.asfortranarray([170.0, 127.0, 92.0, 64.0])

    @staticmethod
    def _call_function_under_test(coeffs, s_val, **kwargs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.bezier_value_check(
            coeffs, s_val, **kwargs)

    def test_check_one_failure(self):
        # 3 (s^2 + 3) (2 - s)
        coeffs = np.asfortranarray([18, 15, 14, 12])
        # 3 (1^2 + 3) (2 - 1) = 12
        rhs_val = 10.0
        self.assertFalse(
            self._call_function_under_test(coeffs, 1.0, rhs_val=rhs_val))

    def test_check_one_success(self):
        # 2 (1 - s) (s + 4)
        coeffs = np.asfortranarray([8.0, 5.0, 0.0])
        self.assertTrue(self._call_function_under_test(coeffs, 1.0))

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_constant(self):
        input_s = (-9.0, -2.0, 0.0, 0.5, 1.0, 1.5, 4.0)
        values = (1.0, 2.0, 4.5, 0.0)
        for index, value in enumerate(values):
            count = index + 1
            coeffs = value * np.ones((count,), order='F')
            for s_val in input_s:
                is_root = self._call_function_under_test(
                    coeffs, s_val, rhs_val=value)
                self.assertTrue(is_root)

                is_root = self._call_function_under_test(
                    coeffs, s_val, rhs_val=value + 1.0)
                self.assertFalse(is_root)

    def test_all_zero(self):
        for num_zeros in (1, 2, 3, 4):
            coeffs = np.zeros((num_zeros,), order='F')
            is_root = self._call_function_under_test(coeffs, 0.0)
            self.assertTrue(is_root)

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_linear(self):
        # 8 s + 1
        coeffs = np.asfortranarray([1.0, 9.0])
        for s_val in (-2.0, -0.125, 0.0, 0.75, 1.0, 9.5):
            rhs_val = 8.0 * s_val + 1.0
            is_root = self._call_function_under_test(
                coeffs, s_val, rhs_val=rhs_val)
            self.assertTrue(is_root)

            is_root = self._call_function_under_test(
                coeffs, s_val, rhs_val=rhs_val + 2.0)
            self.assertFalse(is_root)

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_quadratic(self):
        # 2 s (s + 1)
        coeffs = np.asfortranarray([0.0, 1.0, 4.0])
        s_vals = (-2.0, -1.0, -0.5, 0.0, 0.25, 1.0)
        check_values = [
            self._call_function_under_test(coeffs, s_val)
            for s_val in s_vals
        ]
        self.assertEqual(
            check_values, [False, True, False, True, False, False])

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_quadratic_complex(self):
        # s^2 + 1
        coeffs = np.asfortranarray([1.0, 1.0, 2.0])
        rhs_val = 10.0
        s_vals = (-4.5, -3.0, -1.0, 0.0, 1.0, 3.0, 5.5)
        check_values = [
            self._call_function_under_test(coeffs, s_val, rhs_val=rhs_val)
            for s_val in s_vals
        ]
        self.assertEqual(
            check_values, [False, True, False, False, False, True, False])

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_quadratic_drop_degree(self):
        # 2 (s - 2) (s - 1)
        coeffs = np.asfortranarray([4.0, 1.0, 0.0])
        s_vals = (-2.25, -1.0, 0.0, 1.0, 1.25, 2.0, 4.5)
        check_values = [
            self._call_function_under_test(coeffs, s_val)
            for s_val in s_vals
        ]
        self.assertEqual(
            check_values, [False, False, False, True, False, True, False])

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_quadratic_as_elevated(self):
        # 2 (s + 1)
        coeffs = np.asfortranarray([2.0, 3.0, 4.0])
        s_vals = (-2.0, -1.0, -0.5, 0.0, 0.25, 1.0)
        check_values = [
            self._call_function_under_test(coeffs, s_val)
            for s_val in s_vals
        ]
        self.assertEqual(
            check_values, [False, True, False, False, False, False])

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_cubic(self):
        # -(s - 17) (s - 5) (s - 2)
        s_vals = (0.0, 1.0, 2.0, 3.0, 5.0, 6.0, 16.0, 17.0, 18.0)
        check_values = [
            self._call_function_under_test(self.CUBIC_COEFFS, s_val)
            for s_val in s_vals
        ]
        self.assertEqual(
            check_values,
            [False, False, True, False, True, False, False, True, False],
        )

    def _cubic_wiggle(self, root, max_ulps):
        # -(s - 17) (s - 5) (s - 2)
        coeffs = self.CUBIC_COEFFS
        eps = SPACING(root)
        for delta in six.moves.xrange(-32, 32 + 1):
            s_val = root + delta * eps
            self.assertTrue(self._call_function_under_test(coeffs, s_val))

        s_val1 = root + max_ulps * eps
        s_val2 = root - max_ulps * eps
        s_val3 = root + (max_ulps + 1) * eps
        s_val4 = root - (max_ulps + 1) * eps
        self.assertTrue(self._call_function_under_test(coeffs, s_val1))
        self.assertTrue(self._call_function_under_test(coeffs, s_val2))
        check3 = self._call_function_under_test(coeffs, s_val3)
        check4 = self._call_function_under_test(coeffs, s_val4)
        self.assertFalse(check3 and check4)

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_cubic_wiggle(self):
        self._cubic_wiggle(2.0, 308)
        # This wiggle room is especially egregious.
        if (base_utils.IS_LINUX and
                not base_utils.IS_64_BIT):  # pragma: NO COVER
            self._cubic_wiggle(5.0, 8125)
            self._cubic_wiggle(17.0, 22504)
        else:
            self._cubic_wiggle(5.0, 8129)
            self._cubic_wiggle(17.0, 22503)

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_cubic_drop_degree(self):
        # 3 (1 - s)^2 (3 s + 1)
        coeffs = np.asfortranarray([3.0, 4.0, 0.0, 0.0])
        s_vals = (-2.5, -1.0, -1.0 / 3.0, 0.0, 0.5, 1.0, 1.25)
        check_values = [
            self._call_function_under_test(coeffs, s_val)
            for s_val in s_vals
        ]
        self.assertEqual(
            check_values, [False, False, True, False, False, True, False])

    @unittest.skipIf(SCIPY_LAPACK is None, 'SciPy not installed')
    def test_quartic(self):
        # 4 s^3 (5 - s)
        coeffs = np.asfortranarray([0.0, 0.0, 0.0, 5.0, 16.0])
        s_vals = (-1.0, 0.0, 1.0, 3.5, 5.0, 7.25)
        check_values = [
            self._call_function_under_test(coeffs, s_val)
            for s_val in s_vals
        ]
        self.assertEqual(
            check_values, [False, True, False, False, True, False])


class Test_poly_to_power_basis(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(bezier_coeffs):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.poly_to_power_basis(bezier_coeffs)

    def test_point(self):
        bezier_coeffs = np.asfortranarray([1.0])
        expected = bezier_coeffs
        self.assertEqual(
            self._call_function_under_test(bezier_coeffs), expected)

    def test_line(self):
        bezier_coeffs = np.asfortranarray([1.0, 2.0])
        expected = np.asfortranarray([1.0, 1.0])
        self.assertEqual(
            self._call_function_under_test(bezier_coeffs), expected)

    def test_quadratic(self):
        bezier_coeffs = np.asfortranarray([-1.0, 0.0, 3.0])
        expected = np.asfortranarray([-1.0, 2.0, 2.0])
        self.assertEqual(
            self._call_function_under_test(bezier_coeffs), expected)

    def test_cubic(self):
        bezier_coeffs = np.asfortranarray([1.0, 2.0, 1.0, 2.0])
        expected = np.asfortranarray([1.0, 3.0, -6.0, 4.0])
        self.assertEqual(
            self._call_function_under_test(bezier_coeffs), expected)

    def test_unsupported(self):
        bezier_coeffs = np.asfortranarray([1.0, 0.0, 0.0, 2.0, 1.0])
        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(bezier_coeffs)


class Test_locate_point(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, x_val, y_val):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.locate_point(nodes, x_val, y_val)

    def test_reduced(self):
        nodes = np.asfortranarray([
            [0.0, 3.0],
            [1.0, 2.0],
            [3.0, 0.0],
        ])
        result = self._call_function_under_test(nodes, 1.25, 1.75)
        self.assertEqual(result, 0.5)

    def test_elevated_swap(self):
        nodes = np.asfortranarray([
            [0.0, 1.0],
            [1.0, 2.0],
            [1.0, 3.0],
            [2.0, 4.0],
        ])
        result = self._call_function_under_test(nodes, 1.40625, 3.25)
        self.assertEqual(result, 0.75)

    def test_with_point(self):
        nodes = np.asfortranarray([
            [0.0, 1.0],
            [0.0, 5.0],
        ])
        result = self._call_function_under_test(nodes, 0.0, 1.5)
        self.assertEqual(result, 0.125)

    def test_no_solution_first(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [3.0, 0.0],
        ])
        result = self._call_function_under_test(nodes, 4.0, 0.25)
        self.assertIsNone(result)

    def test_no_solution_second(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [3.0, 0.0],
        ])
        result = self._call_function_under_test(nodes, 0.5, 1.0)
        self.assertIsNone(result)


class Test_all_intersections(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(curve_first, curve_second):
        from bezier import _algebraic_intersection

        return _algebraic_intersection.all_intersections(
            curve_first, curve_second)

    def test_no_intersections(self):
        import bezier

        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve1 = bezier.Curve(nodes1, degree=1, _copy=False)

        nodes2 = np.asfortranarray([
            [3.0, 3.0],
            [4.0, 3.0],
        ])
        curve2 = bezier.Curve(nodes2, degree=1, _copy=False)

        intersections = self._call_function_under_test(curve1, curve2)
        self.assertEqual(intersections.shape, (0, 2))

    def test_success(self):
        import bezier

        # NOTE: ``nodes1`` is a specialization of [0, 0], [1/2, 1], [1, 1]
        #       onto the interval [1/4, 1] and ``nodes`` is a specialization
        #       of [0, 1], [1/2, 1], [1, 0] onto the interval [0, 3/4].
        #       We expect them to intersect at s = 1/3, t = 2/3, which is
        #       the point [1/2, 3/4].
        nodes1 = np.asfortranarray([
            [0.25, 0.4375],
            [0.625, 1.0],
            [1.0, 1.0],
        ])
        curve1 = bezier.Curve(nodes1, degree=2, _copy=False)

        nodes2 = np.asfortranarray([
            [0.0, 1.0],
            [0.375, 1.0],
            [0.75, 0.4375],
        ])
        curve2 = bezier.Curve(nodes2, degree=2, _copy=False)
        s_val = 1.0 / 3.0
        t_val = 2.0 / 3.0

        intersections = self._call_function_under_test(curve1, curve2)
        expected = np.asfortranarray([[s_val, t_val]])
        self.assertEqual(intersections, expected)
