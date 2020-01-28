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
    import scipy.integrate as SCIPY_INT
except ImportError:  # pragma: NO COVER
    SCIPY_INT = None

from tests import utils as base_utils
from tests.unit import utils


FLOAT64 = np.float64  # pylint: disable=no-member
SPACING = np.spacing  # pylint: disable=no-member


class Test_make_subdivision_matrices(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(degree):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.make_subdivision_matrices(degree)

    def _helper(self, degree, expected_l, expected_r):
        left, right = self._call_function_under_test(degree)
        self.assertEqual(left, expected_l)
        self.assertEqual(right, expected_r)

    def test_linear(self):
        from bezier import _py_curve_helpers

        self._helper(
            1,
            _py_curve_helpers._LINEAR_SUBDIVIDE_LEFT,
            _py_curve_helpers._LINEAR_SUBDIVIDE_RIGHT,
        )

    def test_quadratic(self):
        from bezier import _py_curve_helpers

        self._helper(
            2,
            _py_curve_helpers._QUADRATIC_SUBDIVIDE_LEFT,
            _py_curve_helpers._QUADRATIC_SUBDIVIDE_RIGHT,
        )

    def test_cubic(self):
        from bezier import _py_curve_helpers

        self._helper(
            3,
            _py_curve_helpers._CUBIC_SUBDIVIDE_LEFT,
            _py_curve_helpers._CUBIC_SUBDIVIDE_RIGHT,
        )

    def test_quartic(self):
        expected_l = np.asfortranarray(
            [
                [1.0, 1.0, 1.0, 1.0, 1.0],
                [0.0, 1.0, 2.0, 3.0, 4.0],
                [0.0, 0.0, 1.0, 3.0, 6.0],
                [0.0, 0.0, 0.0, 1.0, 4.0],
                [0.0, 0.0, 0.0, 0.0, 1.0],
            ]
        )
        expected_r = np.asfortranarray(
            [
                [1.0, 0.0, 0.0, 0.0, 0.0],
                [4.0, 1.0, 0.0, 0.0, 0.0],
                [6.0, 3.0, 1.0, 0.0, 0.0],
                [4.0, 3.0, 2.0, 1.0, 0.0],
                [1.0, 1.0, 1.0, 1.0, 1.0],
            ]
        )
        col_scaling = np.asfortranarray([[1.0, 2.0, 4.0, 8.0, 16.0]])
        expected_l /= col_scaling
        expected_r /= col_scaling[:, ::-1]
        self._helper(4, expected_l, expected_r)


class Test_subdivide_nodes(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.subdivide_nodes(nodes)

    def _helper(self, nodes, expected_l, expected_r):
        left, right = self._call_function_under_test(nodes)
        self.assertEqual(left, expected_l)
        self.assertEqual(right, expected_r)

    def _points_check(self, nodes, pts_exponent=5):
        from bezier import _py_curve_helpers

        left, right = self._call_function_under_test(nodes)
        # Using the exponent means that ds = 1/2**exp, which
        # can be computed without roundoff.
        num_pts = 2 ** pts_exponent + 1
        left_half = np.linspace(0.0, 0.5, num_pts)
        right_half = np.linspace(0.5, 1.0, num_pts)
        unit_interval = np.linspace(0.0, 1.0, num_pts)
        pairs = [(left, left_half), (right, right_half)]
        for sub_curve, half in pairs:
            # Make sure sub_curve([0, 1]) == curve(half)
            self.assertEqual(
                _py_curve_helpers.evaluate_multi(nodes, half),
                _py_curve_helpers.evaluate_multi(sub_curve, unit_interval),
            )

    def test_line(self):
        nodes = np.asfortranarray([[0.0, 4.0], [1.0, 6.0]])
        expected_l = np.asfortranarray([[0.0, 2.0], [1.0, 3.5]])
        expected_r = np.asfortranarray([[2.0, 4.0], [3.5, 6.0]])
        self._helper(nodes, expected_l, expected_r)

    def test_line_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(shape=(2, 2), seed=88991, num_bits=8)
        self._points_check(nodes)

    def test_quadratic(self):
        nodes = np.asfortranarray([[0.0, 4.0, 7.0], [1.0, 6.0, 3.0]])
        expected_l = np.asfortranarray([[0.0, 2.0, 3.75], [1.0, 3.5, 4.0]])
        expected_r = np.asfortranarray([[3.75, 5.5, 7.0], [4.0, 4.5, 3.0]])
        self._helper(nodes, expected_l, expected_r)

    def test_quadratic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(shape=(2, 3), seed=10764, num_bits=8)
        self._points_check(nodes)

    def test_cubic(self):
        nodes = np.asfortranarray([[0.0, 4.0, 7.0, 6.0], [1.0, 6.0, 3.0, 5.0]])
        expected_l = np.asfortranarray(
            [[0.0, 2.0, 3.75, 4.875], [1.0, 3.5, 4.0, 4.125]]
        )
        expected_r = np.asfortranarray(
            [[4.875, 6.0, 6.5, 6.0], [4.125, 4.25, 4.0, 5.0]]
        )
        self._helper(nodes, expected_l, expected_r)

    def test_cubic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(shape=(2, 4), seed=990077, num_bits=8)
        self._points_check(nodes)

    def test_dynamic_subdivision_matrix(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(shape=(2, 5), seed=103, num_bits=8)
        self._points_check(nodes)


class Test_evaluate_multi_barycentric(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes, lambda1, lambda2):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.evaluate_multi_barycentric(
            nodes, lambda1, lambda2
        )

    def test_non_unity(self):
        nodes = np.asfortranarray(
            [[0.0, 0.5, 1.5, 2.0], [0.0, 3.0, 4.0, 8.0], [0.0, 1.0, 1.0, 1.0]]
        )
        lambda1 = np.asfortranarray([0.25, 0.5, 0.75])
        lambda2 = np.asfortranarray([0.25, 0.125, -0.75])
        result = self._call_function_under_test(nodes, lambda1, lambda2)
        expected = np.asfortranarray(
            [
                [0.125, 0.0859375, 0.421875],
                [0.453125, 0.390625, -2.109375],
                [0.109375, 0.119140625, -0.421875],
            ]
        )
        self.assertEqual(result, expected)


class Test_evaluate_multi(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes, s_vals):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.evaluate_multi(nodes, s_vals)

    def test_linear(self):
        num_vals = 129
        s_vals = np.linspace(0.0, 1.0, num_vals)
        # B(s) = [s + 1, 1 - 2 s, 3 s - 7]
        nodes = np.asfortranarray([[1.0, 2.0], [1.0, -1.0], [-7.0, -4.0]])
        result = self._call_function_under_test(nodes, s_vals)
        expected = np.empty((3, num_vals), order="F")
        expected[0, :] = 1.0 + s_vals
        expected[1, :] = 1.0 - 2.0 * s_vals
        expected[2, :] = -7.0 + 3.0 * s_vals
        self.assertEqual(result, expected)

    def test_quadratic(self):
        num_vals = 65
        s_vals = np.linspace(0.0, 1.0, num_vals)
        # B(s) = [s(4 - s), 2s(2s - 1)]
        nodes = np.asfortranarray([[0.0, 2.0, 3.0], [0.0, -1.0, 2.0]])
        result = self._call_function_under_test(nodes, s_vals)
        expected = np.empty((2, num_vals), order="F")
        expected[0, :] = s_vals * (4.0 - s_vals)
        expected[1, :] = 2.0 * s_vals * (2.0 * s_vals - 1.0)
        self.assertEqual(result, expected)

    def test_binomial_overflow_int32(self):
        s_vals = np.asfortranarray([0.5])
        degree = 30
        nodes = np.eye(degree + 1, order="F")

        expected = np.asfortranarray(
            [
                1.0,
                30.0,
                435.0,
                4060.0,
                27405.0,
                142506.0,
                593775.0,
                2035800.0,
                5852925.0,
                14307150.0,
                30045015.0,
                54627300.0,
                86493225.0,
                119759850.0,
                145422675.0,
                155117520.0,
                145422675.0,
                119759850.0,
                86493225.0,
                54627300.0,
                30045015.0,
                14307150.0,
                5852925.0,
                2035800.0,
                593775.0,
                142506.0,
                27405.0,
                4060.0,
                435.0,
                30.0,
                1.0,
            ]
        )
        evaluated = self._call_function_under_test(nodes, s_vals)
        binomial_coefficients = evaluated.flatten() * 2.0 ** degree
        self.assertEqual(expected, binomial_coefficients)

    @unittest.skipIf(
        base_utils.IS_LINUX and not base_utils.IS_64_BIT,
        "32-bit is skipped on Linux",
    )
    def test_binomial_roundoff(self):
        s_vals = np.asfortranarray([0.5])
        degree = 55
        nodes = np.eye(degree + 1, order="F")

        expected = np.asfortranarray(
            [
                1.0,
                55.0,
                1485.0,
                26235.0,
                341055.0,
                3478761.0,
                28989675.0,
                202927725.0,
                1217566350.0,
                6358402050.0,
                29248649430.0,
                119653565850.0,
                438729741450.0,
                1451182990950.0,
                4353548972850.0,
                11899700525790.0,
                29749251314475.0,
                68248282427325.0,
                144079707346575.0,
                280576272201225.0,
                505037289962205.0,
                841728816603675.0,
                1300853625660225.0,
                1866442158555975.0,
                2488589544741300.0,
                3085851035479212.0,
                3560597348629860.0 - 0.5,
                3824345300380220.0 - 0.5,
                3824345300380220.0 - 0.5,
                3560597348629860.0 - 0.5,
                3085851035479212.0 - 0.5,
                2488589544741300.0 - 0.5,
                1866442158555974.0 + 0.5,
                1300853625660225.0 - 0.5 ** 2,
                841728816603675.0 - 0.5 ** 3,
                505037289962205.0 - 0.5 ** 4,
                280576272201225.0 - 0.5 ** 4,
                144079707346575.0 - 0.5 ** 5,
                68248282427325.0 - 0.5 ** 6,
                29749251314475.0 - 0.5 ** 7,
                11899700525790.0 - 0.5 ** 8,
                4353548972850.0 - 3 * 0.5 ** 11,
                1451182990950.0 - 0.5 ** 11,
                438729741450.0 - 3 * 0.5 ** 14,
                119653565850.0 - 3 * 0.5 ** 16,
                29248649430.0 - 3 * 0.5 ** 18,
                6358402050.0 - 3 * 0.5 ** 20,
                1217566350.0 - 0.5 ** 21,
                202927725.0 - 3 * 0.5 ** 25,
                28989675.0 - 0.5 ** 26,
                3478761.0 - 0.5 ** 29,
                341055.0 - 3 * 0.5 ** 34,
                26235.0 - 0.5 ** 36,
                1485.0 - 0.5 ** 40,
                55.0 - 5 * 0.5 ** 47,
                1.0,
            ],
        )
        evaluated = self._call_function_under_test(nodes, s_vals)
        binomial_coefficients = evaluated.flatten() * 2.0 ** degree
        self.assertEqual(expected, binomial_coefficients)


class Test_vec_size(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes, s_val):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.vec_size(nodes, s_val)

    def test_linear(self):
        nodes = np.asfortranarray([[0.0, 3.0], [0.0, -4.0]])
        size = self._call_function_under_test(nodes, 0.25)
        self.assertEqual(size, 0.25 * 5.0)

    def test_quadratic(self):
        nodes = np.asfortranarray([[0.0, 2.0, 1.0], [0.0, 3.0, 6.0]])
        size = self._call_function_under_test(nodes, 0.5)
        self.assertEqual(size, 3.25)


class Test_compute_length(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.compute_length(nodes)

    def _scipy_skip(self):
        if SCIPY_INT is None:  # pragma: NO COVER
            self.skipTest("SciPy not installed")

    def test_invalid_size(self):
        nodes = np.empty((2, 0), order="F")
        with self.assertRaises(ValueError) as exc_info:
            self._call_function_under_test(nodes)

        expected = ("Curve should have at least one node.",)
        self.assertEqual(exc_info.exception.args, expected)

    def test_degree_zero(self):
        nodes = np.asfortranarray([[0.0], [0.0]])
        length = self._call_function_under_test(nodes)
        self.assertEqual(length, 0.0)

    def test_linear(self):
        nodes = np.asfortranarray([[0.0, 3.0], [0.0, 4.0]])
        length = self._call_function_under_test(nodes)
        self.assertEqual(length, 5.0)

    def test_quadratic(self):
        self._scipy_skip()
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 2.0, 0.0]])
        length = self._call_function_under_test(nodes)
        # pylint: disable=no-member,assignment-from-no-return
        # 2 INT_0^1 SQRT(16 s^2  - 16 s + 5) ds = SQRT(5) + sinh^{-1}(2)/2
        arcs2 = np.arcsinh(2.0)
        # pylint: enable=no-member,assignment-from-no-return
        expected = np.sqrt(5.0) + 0.5 * arcs2
        local_eps = abs(SPACING(expected))
        self.assertAlmostEqual(length, expected, delta=local_eps)

    def test_cubic(self):
        self._scipy_skip()
        nodes = np.asfortranarray([[0.0, 1.0, 2.0, 3.5], [0.0, 2.0, 0.0, 0.0]])
        length = self._call_function_under_test(nodes)
        # x(s) = s (s^2 + 6) / 2
        # y(s) = 6 s (s - 1)^2
        # x'(s)^2 + y'(s)^2 = (9/4)(145s^4 - 384s^3 + 356s^2 - 128s + 20)
        expected = float.fromhex("0x1.05dd184047a7bp+2")
        local_eps = abs(SPACING(expected))
        self.assertAlmostEqual(length, expected, delta=local_eps)


class Test_elevate_nodes(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.elevate_nodes(nodes)

    def test_linear(self):
        nodes = np.asfortranarray([[0.0, 2.0], [0.0, 4.0]])
        result = self._call_function_under_test(nodes)
        expected = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 2.0, 4.0]])
        self.assertEqual(result, expected)

    def test_quadratic(self):
        nodes = np.asfortranarray(
            [[0.0, 3.0, 6.0], [0.5, 0.5, 0.5], [0.75, 3.0, 2.25]]
        )
        result = self._call_function_under_test(nodes)
        expected = np.asfortranarray(
            [
                [0.0, 2.0, 4.0, 6.0],
                [0.5, 0.5, 0.5, 0.5],
                [0.75, 2.25, 2.75, 2.25],
            ]
        )
        self.assertEqual(result, expected)


class Test_de_casteljau_one_round(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes, lambda1, lambda2):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.de_casteljau_one_round(
            nodes, lambda1, lambda2
        )

    def test_it(self):
        nodes = np.asfortranarray([[0.0, 3.0], [1.0, 5.0]])
        result = self._call_function_under_test(nodes, 0.25, 0.75)
        self.assertEqual(result, np.asfortranarray([[2.25], [4.0]]))


class Test_specialize_curve(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes, start, end):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.specialize_curve(nodes, start, end)

    def test_linear(self):
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        result = self._call_function_under_test(nodes, 0.25, 0.75)
        expected = np.asfortranarray([[0.25, 0.75], [0.25, 0.75]])
        self.assertEqual(result, expected)

    def test_against_subdivision(self):
        import bezier

        nodes = np.asfortranarray([[0.0, 1.0, 3.0], [1.0, 6.0, 5.0]])
        curve = bezier.Curve(nodes, 2)
        left, right = curve.subdivide()
        left_nodes = self._call_function_under_test(nodes, 0.0, 0.5)
        self.assertEqual(left.nodes, left_nodes)
        right_nodes = self._call_function_under_test(nodes, 0.5, 1.0)
        self.assertEqual(right.nodes, right_nodes)

    def test_cubic(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0, 1.0, 3.0], [0.0, -1.0, -2.0, 2.0]]
        )
        result = self._call_function_under_test(nodes, 0.125, 0.625)
        expected = (
            np.asfortranarray(
                [[171, 375, 499, 735], [-187, -423, -579, -335]], dtype=FLOAT64
            )
            / 512.0
        )
        self.assertEqual(result, expected)

    def test_quartic(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0, 1.0, 3.0, 3.0], [5.0, 6.0, 7.0, 6.0, 7.0]]
        )
        result = self._call_function_under_test(nodes, 0.5, 0.75)
        expected = np.asfortranarray(
            [
                [1.5625, 1.78125, 2.015625, 2.2578125, 2.47265625],
                [6.375, 6.4375, 6.46875, 6.484375, 6.5234375],
            ]
        )
        self.assertEqual(result, expected)


class Test_evaluate_hodograph(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(s, nodes):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.evaluate_hodograph(s, nodes)

    def test_line(self):
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        first_deriv1 = self._call_function_under_test(0.25, nodes)
        expected = np.asfortranarray(nodes[:, [1]] - nodes[:, [0]])
        self.assertEqual(first_deriv1, expected)
        # Make sure it is the same elsewhere since
        # the derivative curve is degree 0.
        first_deriv2 = self._call_function_under_test(0.75, nodes)
        self.assertEqual(first_deriv1, first_deriv2)

    def test_quadratic(self):
        nodes = np.asfortranarray([[0.0, 0.5, 1.25], [0.0, 1.0, 0.25]])
        # This defines the curve
        #  B(s) = [s(s + 4)/4, s(8 - 7s)/4]
        # B'(s) = [(2 + s)/2, (4 - 7s)/2]
        for s_val in (0.0, 0.25, 0.5, 0.625, 0.875):
            first_deriv = self._call_function_under_test(s_val, nodes)
            self.assertEqual(first_deriv.shape, (2, 1))
            self.assertEqual(first_deriv[0, 0], (2.0 + s_val) / 2.0)
            self.assertEqual(first_deriv[1, 0], (4.0 - 7.0 * s_val) / 2.0)

    def test_cubic(self):
        nodes = np.asfortranarray(
            [[0.0, 0.25, 0.75, 1.25], [0.0, 1.0, 0.5, 1.0]]
        )
        # This defines the curve
        #  B(s) = [s(3 + 3s - s^2)/4, s(5s^2 - 9s + 6)/2]
        # B'(s) = [3(1 + 2s - s^2)/4, 3(5s^2 - 6s + 2)/2]
        for s_val in (0.125, 0.5, 0.75, 1.0, 1.125):
            first_deriv = self._call_function_under_test(s_val, nodes)
            self.assertEqual(first_deriv.shape, (2, 1))
            x_prime = 3.0 * (1.0 + 2.0 * s_val - s_val * s_val) / 4.0
            self.assertEqual(first_deriv[0, 0], x_prime)
            y_prime = 3.0 * (5.0 * s_val * s_val - 6.0 * s_val + 2.0) / 2.0
            self.assertEqual(first_deriv[1, 0], y_prime)


class Test_get_curvature(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes, tangent_vec, s):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.get_curvature(nodes, tangent_vec, s)

    @staticmethod
    def _get_tangent_vec(s, nodes):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.evaluate_hodograph(s, nodes)

    def test_line(self):
        s = 0.5
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        tangent_vec = self._get_tangent_vec(s, nodes)
        result = self._call_function_under_test(nodes, tangent_vec, s)
        self.assertEqual(result, 0.0)

    def test_elevated_line(self):
        s = 0.25
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 0.5, 1.0]])
        tangent_vec = self._get_tangent_vec(s, nodes)
        result = self._call_function_under_test(nodes, tangent_vec, s)
        self.assertEqual(result, 0.0)

    def test_quadratic(self):
        s = 0.5
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 1.0, 0.0]])
        tangent_vec = self._get_tangent_vec(s, nodes)
        result = self._call_function_under_test(nodes, tangent_vec, s)
        self.assertEqual(result, -4.0)


class Test_newton_refine(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes, point, s):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.newton_refine(nodes, point, s)

    def test_it(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0, 3.0, 2.0], [0.0, -1.0, 2.0, 2.0], [0.0, 1.0, 2.0, 4.0]]
        )
        # curve(1/2) = p
        point = np.asfortranarray([[1.75], [0.625], [1.625]])
        new_s = self._call_function_under_test(nodes, point, 0.25)
        self.assertEqual(110.0 * new_s, 57.0)


class Test_locate_point(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes, point):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.locate_point(nodes, point)

    def test_it(self):
        nodes = np.asfortranarray(
            [[0.0, 3.0, 1.0], [0.0, 0.0, 1.0], [0.0, -1.0, 3.0]]
        )
        # C(1/8) = p
        point = np.asfortranarray([[43.0], [1.0], [-11.0]]) / 64
        result = self._call_function_under_test(nodes, point)
        self.assertEqual(result, 0.125)

    def test_no_match(self):
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 1.0, 0.0]])
        point = np.asfortranarray([[0.5], [2.0]])
        self.assertIsNone(self._call_function_under_test(nodes, point))

    def test_failure_on_invalid(self):
        nodes = np.asfortranarray(
            [[0.0, -1.0, 1.0, -0.75], [2.0, 0.0, 1.0, 1.625]]
        )
        point = np.asfortranarray([[-0.25], [1.375]])
        with self.assertRaises(ValueError):
            self._call_function_under_test(nodes, point)

    def test_outside_left(self):
        # Newton's method pushes the value slightly to the left of ``0.0``.
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 0.0]])
        point = np.asfortranarray([[0.0], [0.0]])
        result = self._call_function_under_test(nodes, point)
        self.assertEqual(result, 0.0)

    def test_outside_right(self):
        # Newton's method pushes the value slightly to the right of ``1.0``.
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 0.0]])
        point = np.asfortranarray([[2.0], [0.0]])
        result = self._call_function_under_test(nodes, point)
        self.assertEqual(result, 1.0)


class Test_reduce_pseudo_inverse(utils.NumPyTestCase):
    EPS = 0.5 ** 52

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.reduce_pseudo_inverse(nodes)

    def test_to_constant(self):
        nodes = np.asfortranarray([[-2.0, -2.0], [1.0, 1.0]])
        result = self._call_function_under_test(nodes)
        expected = np.asfortranarray([[-2.0], [1.0]])
        self.assertEqual(result, expected)

    def test_to_linear(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 2.0, 4.0]])
        result = self._call_function_under_test(nodes)
        expected = np.asfortranarray([[0.0, 2.0], [0.0, 4.0]])
        self.assertEqual(result, expected)

    def _actually_inverse_helper(self, degree):
        from bezier import _py_curve_helpers
        from bezier import _py_helpers

        nodes = np.eye(degree + 2, order="F")
        reduction_mat = self._call_function_under_test(nodes)
        id_mat = np.eye(degree + 1, order="F")
        elevation_mat = _py_curve_helpers.elevate_nodes(id_mat)
        result = _py_helpers.matrix_product(elevation_mat, reduction_mat)
        return result, id_mat

    def test_to_linear_actually_inverse(self):
        result, id_mat = self._actually_inverse_helper(1)
        self.assertEqual(result, id_mat)

    def test_from_quadratic_not_elevated(self):
        from bezier import _py_curve_helpers

        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.5, 0.0]])
        result = self._call_function_under_test(nodes)
        expected = np.asfortranarray([[0.0, 2.0], [0.5, 0.5]])
        self.assertEqual(result, expected)
        re_elevated = _py_curve_helpers.elevate_nodes(result)
        self.assertTrue(np.any(nodes != re_elevated))

    def test_to_quadratic(self):
        nodes = np.asfortranarray(
            [
                [0.0, 2.0, 4.0, 6.0],
                [0.5, 0.5, 0.5, 0.5],
                [0.75, 2.25, 2.75, 2.25],
            ]
        )
        result = self._call_function_under_test(nodes)
        expected = np.asfortranarray(
            [[0.0, 3.0, 6.0], [0.5, 0.5, 0.5], [0.75, 3.0, 2.25]]
        )
        self.assertEqual(result, expected)

    def test_to_quadratic_actually_inverse(self):
        result, id_mat = self._actually_inverse_helper(2)
        max_err = np.abs(result - id_mat).max()
        self.assertLess(max_err, self.EPS)

    def test_to_cubic(self):
        nodes = np.asfortranarray([[0.0, 0.75, 2.0, 2.75, 2.0]])
        result = self._call_function_under_test(nodes)
        expected = np.asfortranarray([[0.0, 1.0, 3.0, 2.0]])
        self.assertEqual(result, expected)

    def test_to_cubic_actually_inverse(self):
        result, id_mat = self._actually_inverse_helper(3)
        max_err = np.abs(result - id_mat).max()
        self.assertLess(max_err, self.EPS)

    def test_unsupported_degree(self):
        from bezier import _py_helpers

        degree = 5
        nodes = utils.get_random_nodes(
            shape=(2, degree + 1), seed=3820, num_bits=8
        )
        with self.assertRaises(_py_helpers.UnsupportedDegree) as exc_info:
            self._call_function_under_test(nodes)
        self.assertEqual(exc_info.exception.degree, degree)
        self.assertEqual(exc_info.exception.supported, (1, 2, 3, 4))


class Test_projection_error(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes, projected):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.projection_error(nodes, projected)

    def test_it(self):
        nodes = np.asfortranarray([[0.0, 3.0], [4.0, 0.0]])
        result = self._call_function_under_test(nodes, nodes)
        self.assertEqual(result, 0.0)
        projected = np.asfortranarray([[0.5, 2.5], [4.5, 0.5]])
        result = self._call_function_under_test(nodes, projected)
        self.assertEqual(5.0 * result, 1.0)

    def test_nodes_zero(self):
        nodes = np.asfortranarray([[0.0], [0.0]])
        result = self._call_function_under_test(nodes, nodes)
        self.assertEqual(result, 0.0)


class Test_maybe_reduce(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.maybe_reduce(nodes)

    def _low_degree_helper(self, nodes):
        was_reduced, new_nodes = self._call_function_under_test(nodes)
        self.assertFalse(was_reduced)
        self.assertIs(new_nodes, nodes)

    def test_low_degree(self):
        nodes = np.asfortranarray([[1.0], [1.0]])
        self._low_degree_helper(nodes)
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        self._low_degree_helper(nodes)
        # NOTE: This **should** be reduced, but we don't bother reducing
        #       to a point (since it isn't a curve).
        nodes = np.asfortranarray([[2.0, 2.0], [2.0, 2.0]])
        was_reduced, new_nodes = self._call_function_under_test(nodes)
        self.assertTrue(was_reduced)
        expected = np.asfortranarray([[2.0], [2.0]])
        self.assertEqual(new_nodes, expected)

    def test_to_linear(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [3.0, 3.5, 4.0]])
        was_reduced, new_nodes = self._call_function_under_test(nodes)
        self.assertTrue(was_reduced)
        expected = np.asfortranarray([[0.0, 2.0], [3.0, 4.0]])
        self.assertEqual(expected, new_nodes)

    def test_to_quadratic(self):
        nodes = np.asfortranarray([[3.0, 2.0, 1.0, 0.0], [0.0, 2.0, 2.0, 0.0]])
        was_reduced, new_nodes = self._call_function_under_test(nodes)
        self.assertTrue(was_reduced)
        expected = np.asfortranarray([[3.0, 1.5, 0.0], [0.0, 3.0, 0.0]])
        self.assertEqual(expected, new_nodes)

    def test_from_cubic_not_elevated(self):
        nodes = np.asfortranarray(
            [[0.0, -1.0, 1.0, -0.75], [2.0, 0.0, 1.0, 1.625]]
        )
        was_reduced, new_nodes = self._call_function_under_test(nodes)
        self.assertFalse(was_reduced)
        self.assertIs(new_nodes, nodes)

    def test_to_cubic(self):
        nodes = np.asfortranarray(
            [[0.0, 0.75, 2.0, 3.5, 5.0], [0.0, 1.5, 2.5, 3.0, 3.0]]
        )
        was_reduced, new_nodes = self._call_function_under_test(nodes)
        self.assertTrue(was_reduced)
        expected = np.asfortranarray(
            [[0.0, 1.0, 3.0, 5.0], [0.0, 2.0, 3.0, 3.0]]
        )
        self.assertEqual(expected, new_nodes)

    def test_unsupported_degree(self):
        from bezier import _py_helpers

        degree = 5
        nodes = utils.get_random_nodes(
            shape=(2, degree + 1), seed=77618, num_bits=8
        )
        with self.assertRaises(_py_helpers.UnsupportedDegree) as exc_info:
            self._call_function_under_test(nodes)
        self.assertEqual(exc_info.exception.degree, degree)
        self.assertEqual(exc_info.exception.supported, (0, 1, 2, 3, 4))


class Test_full_reduce(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _py_curve_helpers

        return _py_curve_helpers.full_reduce(nodes)

    def test_linear(self):
        nodes = np.asfortranarray([[5.5, 5.5]])
        new_nodes = self._call_function_under_test(nodes)
        expected = np.asfortranarray([[5.5]])
        self.assertEqual(expected, new_nodes)

    def test_single(self):
        nodes = np.asfortranarray([[0.0, 2.0, 4.0, 6.0], [0.0, 4.0, 6.0, 6.0]])
        new_nodes = self._call_function_under_test(nodes)
        expected = np.asfortranarray([[0.0, 3.0, 6.0], [0.0, 6.0, 6.0]])
        self.assertEqual(expected, new_nodes)

    def test_multiple(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0, 3.0], [4.0, 4.5, 5.0, 5.5]])
        new_nodes = self._call_function_under_test(nodes)
        expected = np.asfortranarray([[0.0, 3.0], [4.0, 5.5]])
        self.assertEqual(expected, new_nodes)

    def test_no_reduce(self):
        nodes = np.asfortranarray(
            [[0.0, -1.0, 1.0, -0.75], [2.0, 0.0, 1.0, 1.625]]
        )
        new_nodes = self._call_function_under_test(nodes)
        self.assertIs(new_nodes, nodes)

    def test_unsupported_degree(self):
        from bezier import _py_helpers

        degree = 5
        nodes = utils.get_random_nodes(
            shape=(2, degree + 1), seed=360009, num_bits=8
        )
        with self.assertRaises(_py_helpers.UnsupportedDegree) as exc_info:
            self._call_function_under_test(nodes)
        self.assertEqual(exc_info.exception.degree, degree)
        self.assertEqual(exc_info.exception.supported, (0, 1, 2, 3, 4))
