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

import mock
import numpy as np
try:
    import scipy.integrate as SCIPY_INT
except ImportError:  # pragma: NO COVER
    SCIPY_INT = None

from tests import utils


class Test_make_subdivision_matrix(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(degree):
        from bezier import _curve_helpers

        return _curve_helpers.make_subdivision_matrix(degree)

    def _helper(self, degree, expected):
        result = self._call_function_under_test(degree)
        self.assertEqual(result, expected)

    def test_linear(self):
        from bezier import curve

        self._helper(1, curve._LINEAR_SUBDIVIDE)

    def test_quadratic(self):
        from bezier import curve

        self._helper(2, curve._QUADRATIC_SUBDIVIDE)

    def test_cubic(self):
        from bezier import curve

        self._helper(3, curve._CUBIC_SUBDIVIDE)


# pylint: disable=no-member
class Base_evaluate_multi(object):  # pylint: disable=invalid-name

    def test_linear(self):
        num_vals = 129
        s_vals = np.linspace(0.0, 1.0, num_vals)
        # B(s) = [s + 1, 1 - 2 s, 3 s - 7]
        nodes = np.array([
            [1.0, 1.0, -7.0],
            [2.0, -1.0, -4.0],
        ])

        result = self._call_function_under_test(nodes, s_vals)

        expected = np.empty((num_vals, 3))
        expected[:, 0] = 1.0 + s_vals
        expected[:, 1] = 1.0 - 2.0 * s_vals
        expected[:, 2] = -7.0 + 3.0 * s_vals

        self.assertEqual(result, expected)

    def test_quadratic(self):
        num_vals = 65
        s_vals = np.linspace(0.0, 1.0, num_vals)
        # B(s) = [s(4 - s), 2s(2s - 1)]
        nodes = np.array([
            [0.0, 0.0],
            [2.0, -1.0],
            [3.0, 2.0],
        ])

        result = self._call_function_under_test(nodes, s_vals)

        expected = np.empty((num_vals, 2))
        expected[:, 0] = s_vals * (4.0 - s_vals)
        expected[:, 1] = 2.0 * s_vals * (2.0 * s_vals - 1.0)

        self.assertEqual(result, expected)
# pylint: enable=no-member


class Test__evaluate_multi(Base_evaluate_multi, utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, s_vals):
        from bezier import _curve_helpers

        return _curve_helpers._evaluate_multi(nodes, s_vals)


@unittest.skipIf(utils.WITHOUT_SPEEDUPS, 'No speedups available')
class Test_speedup_evaluate_multi(Base_evaluate_multi, utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, s_vals):
        from bezier import _speedup

        return _speedup.speedup.evaluate_multi(nodes, s_vals)


class Test__vec_size(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, s_val):
        from bezier import _curve_helpers

        return _curve_helpers._vec_size(nodes, s_val)

    def test_linear(self):
        nodes = np.array([
            [0.0, 0.0],
            [3.0, -4.0],
        ])
        size = self._call_function_under_test(nodes, 0.25)
        self.assertEqual(size, 0.25 * 5.0)

    def test_quadratic(self):
        nodes = np.array([
            [0.0, 0.0],
            [2.0, 3.0],
            [1.0, 6.0],
        ])
        size = self._call_function_under_test(nodes, 0.5)
        self.assertEqual(size, 3.25)


class Test_compute_length(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _curve_helpers

        return _curve_helpers.compute_length(nodes, degree)

    def test_linear(self):
        nodes = np.array([
            [0.0, 0.0],
            [3.0, 4.0],
        ])
        length = self._call_function_under_test(nodes, 1)
        self.assertEqual(length, 5.0)

    @unittest.skipIf(SCIPY_INT is None, 'SciPy not installed')
    def test_quadratic(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
            [2.0, 0.0],
        ])
        length = self._call_function_under_test(nodes, 2)
        # 2 INT_0^1 SQRT(16 s^2  - 16 s + 5) ds = SQRT(5) + sinh^{-1}(2)/2
        arcs2 = np.arcsinh(2.0)  # pylint: disable=no-member
        expected = np.sqrt(5.0) + 0.5 * arcs2
        self.assertLess(abs(length - expected), 1e-15)

    def test_without_scipy(self):
        nodes = np.zeros((5, 2))
        with mock.patch('bezier._curve_helpers._scipy_int', new=None):
            with self.assertRaises(OSError):
                self._call_function_under_test(nodes, 4)


class Test_elevate_nodes(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, dimension):
        from bezier import _curve_helpers

        return _curve_helpers.elevate_nodes(nodes, degree, dimension)

    def test_linear(self):
        nodes = np.array([
            [0.0, 0.0],
            [2.0, 4.0],
        ])
        result = self._call_function_under_test(nodes, 1, 2)
        expected = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
            [2.0, 4.0],
        ])
        self.assertEqual(result, expected)

    def test_quadratic(self):
        nodes = np.array([
            [0.0, 0.5, 0.75],
            [3.0, 0.5, 3.0],
            [6.0, 0.5, 2.25],
        ])
        result = self._call_function_under_test(nodes, 2, 3)
        expected = np.array([
            [0.0, 0.5, 0.75],
            [2.0, 0.5, 2.25],
            [4.0, 0.5, 2.75],
            [6.0, 0.5, 2.25],
        ])
        self.assertEqual(result, expected)


class Test_de_casteljau_one_round(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, lambda1, lambda2):
        from bezier import _curve_helpers

        return _curve_helpers.de_casteljau_one_round(
            nodes, lambda1, lambda2)

    def test_it(self):
        nodes = np.array([
            [0.0, 1.0],
            [3.0, 5.0],
        ])
        result = self._call_function_under_test(nodes, 0.25, 0.75)
        self.assertEqual(result, np.array([[2.25, 4.0]]))


class Test_specialize_curve(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(
            nodes, degree, start, end, curve_start, curve_end):
        from bezier import _curve_helpers

        return _curve_helpers.specialize_curve(
            nodes, degree, start, end, curve_start, curve_end)

    def test_it(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        result, true_start, true_end = self._call_function_under_test(
            nodes, 1, 0.25, 0.75, 0.125, 0.25)
        expected = np.array([
            [0.25, 0.25],
            [0.75, 0.75],
        ])
        self.assertEqual(result, expected)
        self.assertEqual(true_start, 0.15625)
        self.assertEqual(true_end, 0.21875)

    def test_againt_subdivision(self):
        import bezier

        nodes = np.array([
            [0.0, 1.0],
            [1.0, 6.0],
            [3.0, 5.0],
        ])
        curve = bezier.Curve(nodes, 2)
        left, right = curve.subdivide()

        left_nodes, true_start, true_end = self._call_function_under_test(
            nodes, 2, 0.0, 0.5, 0.0, 1.0)
        self.assertEqual(left.nodes, left_nodes)
        self.assertEqual(true_start, left.start)
        self.assertEqual(true_end, left.end)

        right_nodes, true_start, true_end = self._call_function_under_test(
            nodes, 2, 0.5, 1.0, 0.0, 1.0)
        self.assertEqual(right.nodes, right_nodes)
        self.assertEqual(true_start, right.start)
        self.assertEqual(true_end, right.end)


class Test_evaluate_hodograph(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, s):
        from bezier import _curve_helpers

        return _curve_helpers.evaluate_hodograph(nodes, degree, s)

    def test_line(self):
        degree = 1
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])

        first_deriv1 = self._call_function_under_test(nodes, degree, 0.25)
        self.assertEqual(first_deriv1.shape, (1, 2))
        self.assertEqual(first_deriv1, nodes[[1], :] - nodes[[0], :])
        # Make sure it is the same elsewhere since
        # the derivative curve is degree 0.
        first_deriv2 = self._call_function_under_test(nodes, degree, 0.75)
        self.assertEqual(first_deriv1, first_deriv2)

    def test_quadratic(self):
        degree = 2
        nodes = np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.25, 0.25],
        ])
        # This defines the curve
        #  B(s) = [s(s + 4)/4, s(8 - 7s)/4]
        # B'(s) = [(2 + s)/2, (4 - 7s)/2]

        for s_val in (0.0, 0.25, 0.5, 0.625, 0.875):
            first_deriv = self._call_function_under_test(nodes, degree, s_val)
            self.assertEqual(first_deriv.shape, (1, 2))
            self.assertEqual(first_deriv[0, 0], (2.0 + s_val) / 2.0)
            self.assertEqual(first_deriv[0, 1], (4.0 - 7.0 * s_val) / 2.0)

    def test_cubic(self):
        degree = 3
        nodes = np.array([
            [0.0, 0.0],
            [0.25, 1.0],
            [0.75, 0.5],
            [1.25, 1.0],
        ])
        # This defines the curve
        #  B(s) = [s(3 + 3s - s^2)/4, s(5s^2 - 9s + 6)/2]
        # B'(s) = [3(1 + 2s - s^2)/4, 3(5s^2 - 6s + 2)/2]
        for s_val in (0.125, 0.5, 0.75, 1.0, 1.125):
            first_deriv = self._call_function_under_test(nodes, degree, s_val)
            self.assertEqual(first_deriv.shape, (1, 2))
            x_prime = 3.0 * (1.0 + 2.0 * s_val - s_val * s_val) / 4.0
            self.assertEqual(first_deriv[0, 0], x_prime)
            y_prime = 3.0 * (5.0 * s_val * s_val - 6.0 * s_val + 2.0) / 2.0
            self.assertEqual(first_deriv[0, 1], y_prime)


class Test_get_curvature(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, tangent_vec, s):
        from bezier import _curve_helpers

        return _curve_helpers.get_curvature(
            nodes, degree, tangent_vec, s)

    @staticmethod
    def _get_tangent_vec(nodes, degree, s):
        from bezier import _curve_helpers

        return _curve_helpers.evaluate_hodograph(nodes, degree, s)

    def test_line(self):
        s = 0.5
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        tangent_vec = self._get_tangent_vec(nodes, 1, s)
        result = self._call_function_under_test(nodes, 1, tangent_vec, s)
        self.assertEqual(result, 0.0)

    def test_elevated_line(self):
        s = 0.25
        nodes = np.array([
            [0.0, 0.0],
            [0.5, 0.5],
            [1.0, 1.0],
        ])
        tangent_vec = self._get_tangent_vec(nodes, 2, s)
        result = self._call_function_under_test(nodes, 2, tangent_vec, s)
        self.assertEqual(result, 0.0)

    def test_quadratic(self):
        s = 0.5
        nodes = np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ])
        tangent_vec = self._get_tangent_vec(nodes, 2, s)
        result = self._call_function_under_test(nodes, 2, tangent_vec, s)
        self.assertEqual(result, -4.0)


class Test_newton_refine(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(curve, point, s):
        from bezier import _curve_helpers

        return _curve_helpers.newton_refine(curve, point, s)

    def test_it(self):
        import bezier

        curve = bezier.Curve.from_nodes(np.array([
            [0.0, 0.0, 0.0],
            [1.0, -1.0, 1.0],
            [3.0, 2.0, 2.0],
            [2.0, 2.0, 4.0],
        ]))
        point = curve.evaluate_multi(np.array([0.5]))
        new_s = self._call_function_under_test(curve, point, 0.25)
        self.assertEqual(110.0 * new_s, 57.0)


class Test_locate_point(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(curve, point):
        from bezier import _curve_helpers

        return _curve_helpers.locate_point(curve, point)

    def test_it(self):
        import bezier

        curve = bezier.Curve.from_nodes(np.array([
            [0.0, 0.0, 0.0],
            [3.0, 0.0, -1.0],
            [1.0, 1.0, 3.0],
        ]))
        point = curve.evaluate_multi(np.array([0.125]))
        result = self._call_function_under_test(curve, point)
        self.assertEqual(result, 0.125)

    def test_no_match(self):
        import bezier

        curve = bezier.Curve.from_nodes(np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ]))
        point = np.array([[0.5, 2.0]])
        self.assertIsNone(self._call_function_under_test(curve, point))

    def test_failure_on_invalid(self):
        import bezier

        nodes = np.array([
            [0.0, 2.0],
            [-1.0, 0.0],
            [1.0, 1.0],
            [-0.75, 1.625],
        ])
        curve = bezier.Curve(nodes, 3)
        point = np.array([[-0.25, 1.375]])
        with self.assertRaises(ValueError):
            self._call_function_under_test(curve, point)
