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


FLOAT64 = np.float64  # pylint: disable=no-member
LOCAL_EPS = 0.5**26  # sqrt(machine precision)


class Test__evaluate3(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, x_val, y_val):
        from bezier import _implicitization

        return _implicitization._evaluate3(nodes, x_val, y_val)

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
        from bezier import _implicitization

        return _implicitization.eval_intersection_polynomial(
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
        from bezier import _implicitization

        return _implicitization._to_power_basis11(nodes1, nodes2)

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
        expected = np.array([-7.0, 15.0]) / 32.0
        self.assertEqual(result, expected)


class Test__to_power_basis12(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _implicitization

        return _implicitization._to_power_basis12(nodes1, nodes2)

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
        expected = np.array([-0.5, 2.0, -2.0])
        self.assertEqual(result, expected)


class Test__to_power_basis13(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _implicitization

        return _implicitization._to_power_basis13(nodes1, nodes2)

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
        expected = 3.0 * (-63.0) * np.array([-1.0, 9.0, -24.0, 16.0]) / 32.0
        self.assertEqual(result, expected)


class Test__to_power_basis22(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _implicitization

        return _implicitization._to_power_basis22(nodes1, nodes2)

    def test_it(self):
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
        # The 2-2 method avoids a division by 3.0
        expected = 3.0 * np.array([-3.0, 16.0, -24.0, 0.0, 16.0]) / 256.0
        self.assertEqual(result, expected)


class Test__to_power_basis23(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _implicitization

        return _implicitization._to_power_basis23(nodes1, nodes2)

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
        expected = np.array([44, -192, 240, -64, 0.0, 0.0, 0.0])
        self.assertLess(np.abs(result - expected).max(), LOCAL_EPS)


class Test__to_power_basis33(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _implicitization

        return _implicitization._to_power_basis33(nodes1, nodes2)

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
        expected = np.array(
            [-62, 81, 0, -12, 0, 0, -6, 0, 0, -1], dtype=FLOAT64)
        self.assertLess(np.abs(result - expected).max(), LOCAL_EPS)


class Test_to_power_basis(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _implicitization

        return _implicitization.to_power_basis(nodes1, nodes2)

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
        expected = np.array([-1.0, 0.0])
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
        expected = np.array([-48.0, 180.0, -108.0])
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
        expected = 3.0 * np.array([0.75, 0.0, -0.75, -0.25])
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
        # The 2-2 method avoids a division by 3.0
        expected = 3.0 * np.array([12.0, -28.0, 5.0, -18.0, 9.0]) / 16.0
        self.assertEqual(result, expected)

    def test_degrees_2_3(self):
        # f1(x, y) = 81 (2 x^2 - 2 x - y + 1) / 128
        nodes1 = np.asfortranarray([
            [0.25, 0.625],
            [0.625, 0.25],
            [1.0, 1.0],
        ])
        # x2(t), y2(t) = -t (2 t^2 - 3 t - 3) / 4, y2(t) = -(3 t^3 - 3 t - 1) / 2
        nodes2 = np.asfortranarray([
            [0.0, 0.5],
            [0.25, 1.0],
            [0.75, 1.5],
            [1.0, 0.5],
        ])
        # f1(x2(t), y2(t)) = (
        #     81 (t + 1) (4 t^5 - 16 t^4 + 13 t^3 + 25 t^2 - 28 t + 4) / 1024)
        result = self._call_function_under_test(nodes1, nodes2)
        expected = 81.0 * np.array([4, -24, -3, 38, -3, -12, 4]) / 1024.0
        self.assertTrue(np.allclose(result, expected, atol=0.0, rtol=LOCAL_EPS))

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
        expected = np.array([
            14107, -158418, 10827, 331686, -146394,
            -171072, 170451, 1620, -48600, 13500]) / 2048.0
        self.assertTrue(np.allclose(result, expected, atol=0.0, rtol=LOCAL_EPS))

    def test_unsupported(self):
        nodes_yes1 = np.zeros((2, 2), order='F')
        nodes_yes2 = np.zeros((3, 2), order='F')
        nodes_yes3 = np.zeros((4, 2), order='F')
        nodes_no = np.zeros((5, 2), order='F')

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
