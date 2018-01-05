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
import pytest

try:
    from bezier import _HAS_SURFACE_SPEEDUP as HAS_SURFACE_SPEEDUP
except ImportError:  # pragma: NO COVER
    HAS_SURFACE_SPEEDUP = False
from tests.unit import utils


UNIT_TRIANGLE = np.asfortranarray([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
])
FLOAT64 = np.float64  # pylint: disable=no-member
# pylint: disable=invalid-name,no-member
slow = pytest.mark.skipif(
    pytest.config.getoption('--ignore-slow') and not HAS_SURFACE_SPEEDUP,
    reason='--ignore-slow ignores the slow tests',
)
# pylint: enable=invalid-name,no-member


class Test_polynomial_sign(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(poly_surface, degree):
        from bezier import _surface_helpers

        return _surface_helpers.polynomial_sign(poly_surface, degree)

    def test_positive(self):
        bernstein = np.asfortranarray(
            [[1.0], [2.0], [3.0], [4.0], [5.0], [6.0]])
        sign = self._call_function_under_test(bernstein, 2)
        self.assertEqual(sign, 1)

    def test_negative(self):
        bernstein = np.asfortranarray([[-1.0], [-2.0], [-1.0]])
        sign = self._call_function_under_test(bernstein, 1)
        self.assertEqual(sign, -1)

    def test_zero(self):
        bernstein = np.zeros((10, 1), order='F')
        sign = self._call_function_under_test(bernstein, 3)
        self.assertEqual(sign, 0)

    def test_mixed(self):
        bernstein = np.asfortranarray([[-1.0], [1.0], [-1.0]])
        sign = self._call_function_under_test(bernstein, 1)
        self.assertEqual(sign, 0)

    def test_max_iterations(self):
        bernstein = np.asfortranarray([[1.0], [2.0], [3.0]])
        subs = 'bezier._surface_helpers._MAX_POLY_SUBDIVISIONS'
        with unittest.mock.patch(subs, new=1):
            sign = self._call_function_under_test(bernstein, 1)
            self.assertEqual(sign, 1)

    def test_no_conclusion(self):
        bernstein = np.asfortranarray([[-1.0], [1.0], [2.0]])
        subs = 'bezier._surface_helpers._MAX_POLY_SUBDIVISIONS'
        with unittest.mock.patch(subs, new=0):
            with self.assertRaises(ValueError):
                self._call_function_under_test(bernstein, 1)

    def test_conclusion_from_corner_node(self):
        # NOTE: This comes from the surface defined by
        #          [0.0 , 0.0  ]
        #          [0.5 , 0.5  ]
        #          [1.0 , 0.625]
        #          [0.0 , 0.5  ]
        #          [0.5 , 0.5  ]
        #          [0.25, 1.0  ]
        bernstein = np.asfortranarray([
            [1.0], [0.5], [0.0], [0.75], [0.4375], [1.0]])
        sign = self._call_function_under_test(bernstein, 2)
        self.assertEqual(sign, 0)


class Test_two_by_two_det(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(mat):
        from bezier import _surface_helpers

        return _surface_helpers.two_by_two_det(mat)

    def test_integers(self):
        mat = np.asfortranarray([
            [1.0, 2.0],
            [3.0, 4.0],
        ])
        self.assertEqual(self._call_function_under_test(mat), -2.0)

    def test_better_than_numpy(self):
        mat = np.asfortranarray([
            [-24.0, 3.0],
            [-27.0, 0.0],
        ]) / 16.0
        actual_det = self._call_function_under_test(mat)
        self.assertEqual(actual_det, 81.0 / 256.0)

        np_det = np.linalg.det(mat)
        self.assertNotEqual(actual_det, np_det)
        self.assertLess(abs(actual_det - np_det), 1e-16)


class Test_quadratic_jacobian_polynomial(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _surface_helpers

        return _surface_helpers.quadratic_jacobian_polynomial(nodes)

    def test_it(self):
        # B(L1, L2, L3) = [L1^2 + L2^2, L2^2 + L3^2]
        nodes = np.asfortranarray([
            [1.0, 0.0],
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
        ])
        bernstein = self._call_function_under_test(nodes)
        expected = np.asfortranarray([
            [0.0], [2.0], [0.0], [-2.0], [2.0], [0.0]])
        self.assertEqual(bernstein, expected)

    def test_against_det(self):
        from bezier import _surface_helpers

        # B(L1, L2, L3) = [s (t + 2), s^2 + 4 t]
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 1.0],
            [0.0, 2.0],
            [1.5, 2.0],
            [0.0, 4.0],
        ])

        st_vals = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.0, 0.5],
            [0.5, 0.5],
            [0.0, 1.0],
        ])
        as_det = _surface_helpers.jacobian_det(nodes, 2, st_vals)
        as_det = as_det.reshape((6, 1), order='F')
        # B_s = [t + 2, 2*s]
        # B_t = [s, 4]
        # det(DB) = -2 (s^2 - 2t - 4)
        bernstein = self._call_function_under_test(nodes)
        evaluated_bernstein = _surface_helpers.evaluate_cartesian_multi(
            bernstein, 2, st_vals, 1)
        self.assertEqual(evaluated_bernstein, as_det)


class Test_cubic_jacobian_polynomial(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _surface_helpers

        return _surface_helpers.cubic_jacobian_polynomial(nodes)

    def test_it(self):
        # B(L1, L2, L3) = [L1^3 + L2^3, L2^3 + L3^3]
        nodes = np.asfortranarray([
            [1.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
        ])
        bernstein = self._call_function_under_test(nodes)
        shape = (15, 1)
        self.assertEqual(bernstein.shape, shape)
        expected = np.zeros(shape, order='F')
        expected[2, 0] = 1.5
        expected[9, 0] = -1.5
        expected[11, 0] = 1.5
        self.assertEqual(bernstein, expected)


class Test__de_casteljau_one_round(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, lambda1, lambda2, lambda3):
        from bezier import _surface_helpers

        return _surface_helpers._de_casteljau_one_round(
            nodes, degree, lambda1, lambda2, lambda3)

    def test_linear(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        s_val, t_val = 0.5, 0.375
        expected = np.asfortranarray([
            [s_val, t_val],
        ])

        result = self._call_function_under_test(
            nodes, 1, 1.0 - s_val - t_val, s_val, t_val)
        self.assertEqual(result, expected)

    def test_quadratic(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(6, 2), seed=97764, num_bits=8)

        p200, p110, p020, p101, p011, p002 = nodes
        s_val = 0.25
        t_val = 0.125

        expected = np.asfortranarray([
            (1.0 - s_val - t_val) * p200 + s_val * p110 + t_val * p101,
            (1.0 - s_val - t_val) * p110 + s_val * p020 + t_val * p011,
            (1.0 - s_val - t_val) * p101 + s_val * p011 + t_val * p002,
        ])
        result = self._call_function_under_test(
            nodes, 2, 1.0 - s_val - t_val, s_val, t_val)
        self.assertEqual(result, expected)

    def test_cubic(self):
        from bezier import _helpers

        nodes = np.asfortranarray([
            [0.0, 0.0],
            [3.25, 1.5],
            [6.5, 1.5],
            [10.0, 0.0],
            [1.5, 3.25],
            [5.0, 5.0],
            [10.0, 5.25],
            [1.5, 6.5],
            [5.25, 10.0],
            [0.0, 10.0],
        ])

        s_val = 0.25
        t_val = 0.375
        lambda1 = 1.0 - s_val - t_val
        transform = np.asfortranarray([
            [lambda1, s_val, 0., 0., t_val, 0., 0., 0., 0., 0.],
            [0., lambda1, s_val, 0., 0., t_val, 0., 0., 0., 0.],
            [0., 0., lambda1, s_val, 0., 0., t_val, 0., 0., 0.],
            [0., 0., 0., 0., lambda1, s_val, 0., t_val, 0., 0.],
            [0., 0., 0., 0., 0., lambda1, s_val, 0., t_val, 0.],
            [0., 0., 0., 0., 0., 0., 0., lambda1, s_val, t_val],
        ])
        expected = _helpers.matrix_product(transform, nodes)
        result = self._call_function_under_test(
            nodes, 3, lambda1, s_val, t_val)
        self.assertEqual(result, expected)


@utils.needs_surface_speedup
class Test_speedup_de_casteljau_one_round(Test__de_casteljau_one_round):

    @staticmethod
    def _call_function_under_test(nodes, degree, lambda1, lambda2, lambda3):
        from bezier import _surface_speedup

        return _surface_speedup.de_casteljau_one_round(
            nodes, degree, lambda1, lambda2, lambda3)


class Test_make_transform(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(degree, weights_a, weights_b, weights_c):
        from bezier import _surface_helpers

        return _surface_helpers.make_transform(
            degree, weights_a, weights_b, weights_c)

    def _helper(self, degree, weights, expected0, expected1, expected2):
        result = self._call_function_under_test(
            degree, weights[0, :], weights[1, :], weights[2, :])

        self.assertIsInstance(result, dict)
        self.assertEqual(len(result), 3)
        self.assertEqual(result[0], expected0)
        self.assertEqual(result[1], expected1)
        self.assertEqual(result[2], expected2)

    def test_linear(self):
        weights = np.asfortranarray([
            [1.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
        ])
        expected0 = np.asfortranarray(weights[[0], :])
        expected1 = np.asfortranarray(weights[[1], :])
        expected2 = np.asfortranarray(weights[[2], :])
        self._helper(1, weights, expected0, expected1, expected2)

    def test_quadratic(self):
        weights = np.asfortranarray([
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ])
        expected0 = np.asfortranarray([
            [0.0, 0.5, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, 0.5, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5, 0.5],
        ])
        expected1 = np.asfortranarray([
            [0.5, 0.0, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.0, 0.5],
        ])
        expected2 = np.asfortranarray([
            [0.5, 0.5, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.5, 0.0],
        ])
        self._helper(2, weights, expected0, expected1, expected2)


class Test_reduced_to_matrix(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(shape, degree, vals_by_weight):
        from bezier import _surface_helpers

        return _surface_helpers.reduced_to_matrix(
            shape, degree, vals_by_weight)

    def test_it(self):
        shape = (6, 2)
        degree = 2

        expected = np.asfortranarray([
            [1.0, 0.0],
            [-1.0, 1.0],
            [0.0, 1.0],
            [0.0, -1.0],
            [-1.0, -1.0],
            [2.0, 0.0],
        ])
        vals_by_weight = {
            (0, 0): expected[[0], :],
            (0, 1): expected[[1], :],
            (1, 1): expected[[2], :],
            (0, 2): expected[[3], :],
            (1, 2): expected[[4], :],
            (2, 2): expected[[5], :],
        }

        result = self._call_function_under_test(shape, degree, vals_by_weight)
        self.assertEqual(result, expected)


class Test__specialize_surface(utils.NumPyTestCase):

    WEIGHTS0 = np.asfortranarray([1.0, 0.0, 0.0])
    WEIGHTS1 = np.asfortranarray([0.5, 0.5, 0.0])
    WEIGHTS2 = np.asfortranarray([0.0, 1.0, 0.0])
    WEIGHTS3 = np.asfortranarray([0.5, 0.0, 0.5])
    WEIGHTS4 = np.asfortranarray([0.0, 0.5, 0.5])
    WEIGHTS5 = np.asfortranarray([0.0, 0.0, 1.0])

    @staticmethod
    def _call_function_under_test(
            nodes, degree, weights_a, weights_b, weights_c):
        from bezier import _surface_helpers

        return _surface_helpers._specialize_surface(
            nodes, degree, weights_a, weights_b, weights_c)

    def _helper(self, degree, expected_a, expected_b, expected_c, expected_d):
        from bezier import _helpers

        num_nodes = ((degree + 1) * (degree + 2)) // 2
        id_mat = _helpers.eye(num_nodes)

        computed_a = self._call_function_under_test(
            id_mat, degree,
            self.WEIGHTS0, self.WEIGHTS1, self.WEIGHTS3)
        computed_b = self._call_function_under_test(
            id_mat, degree,
            self.WEIGHTS4, self.WEIGHTS3, self.WEIGHTS1)
        computed_c = self._call_function_under_test(
            id_mat, degree,
            self.WEIGHTS1, self.WEIGHTS2, self.WEIGHTS4)
        computed_d = self._call_function_under_test(
            id_mat, degree,
            self.WEIGHTS3, self.WEIGHTS4, self.WEIGHTS5)

        self.assertEqual(computed_a, expected_a)
        self.assertEqual(computed_b, expected_b)
        self.assertEqual(computed_c, expected_c)
        self.assertEqual(computed_d, expected_d)

    def test_known_linear(self):
        from bezier import _surface_helpers

        self._helper(
            1, _surface_helpers.LINEAR_SUBDIVIDE_A,
            _surface_helpers.LINEAR_SUBDIVIDE_B,
            _surface_helpers.LINEAR_SUBDIVIDE_C,
            _surface_helpers.LINEAR_SUBDIVIDE_D)

    def test_known_quadratic(self):
        from bezier import _surface_helpers

        self._helper(
            2, _surface_helpers.QUADRATIC_SUBDIVIDE_A,
            _surface_helpers.QUADRATIC_SUBDIVIDE_B,
            _surface_helpers.QUADRATIC_SUBDIVIDE_C,
            _surface_helpers.QUADRATIC_SUBDIVIDE_D)

    def test_known_cubic(self):
        from bezier import _surface_helpers

        self._helper(
            3, _surface_helpers.CUBIC_SUBDIVIDE_A,
            _surface_helpers.CUBIC_SUBDIVIDE_B,
            _surface_helpers.CUBIC_SUBDIVIDE_C,
            _surface_helpers.CUBIC_SUBDIVIDE_D)

    def test_known_quartic(self):
        from bezier import _surface_helpers

        self._helper(
            4, _surface_helpers.QUARTIC_SUBDIVIDE_A,
            _surface_helpers.QUARTIC_SUBDIVIDE_B,
            _surface_helpers.QUARTIC_SUBDIVIDE_C,
            _surface_helpers.QUARTIC_SUBDIVIDE_D)


@utils.needs_surface_speedup
class Test_speedup_specialize_surface(Test__specialize_surface):

    @staticmethod
    def _call_function_under_test(
            nodes, degree, weights_a, weights_b, weights_c):
        from bezier import _surface_speedup

        return _surface_speedup.specialize_surface(
            nodes, degree, weights_a, weights_b, weights_c)


class Test__subdivide_nodes(utils.NumPyTestCase):

    REF_TRIANGLE = utils.ref_triangle_uniform_nodes(5)

    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _surface_helpers

        return _surface_helpers._subdivide_nodes(nodes, degree)

    def _helper(self, nodes, degree, expected_a,
                expected_b, expected_c, expected_d):
        nodes_a, nodes_b, nodes_c, nodes_d = self._call_function_under_test(
            nodes, degree)

        self.assertEqual(nodes_a, expected_a)
        self.assertEqual(nodes_b, expected_b)
        self.assertEqual(nodes_c, expected_c)
        self.assertEqual(nodes_d, expected_d)

    def _points_check(self, nodes, degree):
        from bezier import _surface_helpers

        _, dimension = nodes.shape
        sub_surfaces = self._call_function_under_test(
            nodes, degree)

        ref_triangle = self.REF_TRIANGLE
        quarter_a = 0.5 * ref_triangle
        quarters = [
            quarter_a,
            np.asfortranarray([0.5, 0.5]) - quarter_a,  # B
            quarter_a + np.asfortranarray([0.5, 0.0]),  # C
            quarter_a + np.asfortranarray([0.0, 0.5]),  # D
        ]

        for sub_surface, quarter in zip(sub_surfaces, quarters):
            # Make sure sub_surface(ref_triangle) == surface(quarter)
            main_vals = _surface_helpers.evaluate_cartesian_multi(
                nodes, degree, quarter, dimension)
            sub_vals = _surface_helpers.evaluate_cartesian_multi(
                sub_surface, degree, ref_triangle, dimension)
            self.assertEqual(main_vals, sub_vals)

    def test_linear(self):
        expected_a = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5],
        ])
        expected_b = np.asfortranarray([
            [0.5, 0.5],
            [0.0, 0.5],
            [0.5, 0.0],
        ])
        expected_c = np.asfortranarray([
            [0.5, 0.0],
            [1.0, 0.0],
            [0.5, 0.5],
        ])
        expected_d = np.asfortranarray([
            [0.0, 0.5],
            [0.5, 0.5],
            [0.0, 1.0],
        ])
        self._helper(
            UNIT_TRIANGLE, 1,
            expected_a, expected_b, expected_c, expected_d)

    @slow
    def test_line_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(3, 2), seed=123987, num_bits=8)
        self._points_check(nodes, 1)

    def test_quadratic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.25],
            [1.0, 0.0],
            [0.5, 0.75],
            [0.0, 1.0],
            [0.0, 0.5],
        ])
        expected_a = np.asfortranarray([
            [0.0, 0.0],
            [0.25, 0.125],
            [0.5, 0.125],
            [0.25, 0.375],
            [0.25, 0.5],
            [0.25, 0.5],
        ])
        expected_b = np.asfortranarray([
            [0.25, 0.625],
            [0.25, 0.625],
            [0.25, 0.5],
            [0.5, 0.5],
            [0.25, 0.5],
            [0.5, 0.125],
        ])
        expected_c = np.asfortranarray([
            [0.5, 0.125],
            [0.75, 0.125],
            [1.0, 0.0],
            [0.5, 0.5],
            [0.5, 0.5],
            [0.25, 0.625],
        ])
        expected_d = np.asfortranarray([
            [0.25, 0.5],
            [0.25, 0.625],
            [0.25, 0.625],
            [0.25, 0.625],
            [0.0, 0.75],
            [0.0, 0.5],
        ])
        self._helper(
            nodes, 2, expected_a, expected_b, expected_c, expected_d)

    @slow
    def test_quadratic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(6, 2), seed=45001, num_bits=8)
        self._points_check(nodes, 2)

    def test_cubic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [3.25, 1.5],
            [6.5, 1.5],
            [10.0, 0.0],
            [1.5, 3.25],
            [5.0, 5.0],
            [10.0, 5.25],
            [1.5, 6.5],
            [5.25, 10.0],
            [0.0, 10.0],
        ])
        expected_a = np.asfortranarray([
            [0.0, 0.0],
            [1.625, 0.75],
            [3.25, 1.125],
            [4.90625, 1.125],
            [0.75, 1.625],
            [2.4375, 2.4375],
            [4.3125, 2.875],
            [1.125, 3.25],
            [2.875, 4.3125],
            [1.125, 4.90625],
        ])
        expected_b = np.asfortranarray([
            [6.96875, 6.96875],
            [4.8125, 6.65625],
            [2.875, 5.96875],
            [1.125, 4.90625],
            [6.65625, 4.8125],
            [4.75, 4.75],
            [2.875, 4.3125],
            [5.96875, 2.875],
            [4.3125, 2.875],
            [4.90625, 1.125],
        ])
        expected_c = np.asfortranarray([
            [4.90625, 1.125],
            [6.5625, 1.125],
            [8.25, 0.75],
            [10.0, 0.0],
            [5.96875, 2.875],
            [7.875, 2.9375],
            [10.0, 2.625],
            [6.65625, 4.8125],
            [8.8125, 5.125],
            [6.96875, 6.96875],
        ])
        expected_d = np.asfortranarray([
            [1.125, 4.90625],
            [2.875, 5.96875],
            [4.8125, 6.65625],
            [6.96875, 6.96875],
            [1.125, 6.5625],
            [2.9375, 7.875],
            [5.125, 8.8125],
            [0.75, 8.25],
            [2.625, 10.0],
            [0.0, 10.0],
        ])
        self._helper(
            nodes, 3, expected_a, expected_b, expected_c, expected_d)

    @slow
    def test_cubic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(10, 2), seed=346323, num_bits=8)
        self._points_check(nodes, 3)

    @slow
    def test_quartic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(15, 2), seed=741002, num_bits=8)
        self._points_check(nodes, 4)

    @slow
    def test_on_the_fly(self):
        # Test for a degree where the subdivision is done on the fly
        # rather than via a stored matrix.
        nodes = utils.get_random_nodes(
            shape=(21, 2), seed=446, num_bits=8)
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        self._points_check(nodes, 5)


@utils.needs_surface_speedup
class Test_speedup_subdivide_nodes(Test__subdivide_nodes):

    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _surface_speedup

        return _surface_speedup.subdivide_nodes(nodes, degree)


class Test_jacobian_s(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, dimension):
        from bezier import _surface_helpers

        return _surface_helpers.jacobian_s(nodes, degree, dimension)

    def test_linear(self):
        nodes = np.asfortranarray([[0.0], [1.0], [np.nan]])
        result = self._call_function_under_test(nodes, 1, 1)
        expected = np.asfortranarray([[1.0]])
        self.assertEqual(result, expected)

    def test_quadratic(self):
        nodes = np.asfortranarray([
            [0.0, 1.0],
            [1.0, 11.0],
            [5.0, 7.0],
            [4.0, -2.0],
            [-1.0, 6.0],
            [np.nan, np.nan],
        ])
        result = self._call_function_under_test(nodes, 2, 2)
        expected = 2.0 * np.asfortranarray([
            [1.0, 10.0],
            [4.0, -4.0],
            [-5.0, 8.0],
        ])
        self.assertEqual(result, expected)

    def test_cubic(self):
        nodes = np.arange(10, dtype=FLOAT64)[:, np.newaxis]**2
        result = self._call_function_under_test(nodes, 3, 1)
        expected = 3 * np.asfortranarray([
            [1], [3], [5], [9], [11], [15]], dtype=FLOAT64)
        self.assertEqual(result, expected)

    def test_quartic(self):
        nodes = np.arange(15, dtype=FLOAT64)[:, np.newaxis]**2
        result = self._call_function_under_test(nodes, 4, 1)
        expected = 4 * np.asfortranarray([
            [1], [3], [5], [7], [11], [13],
            [15], [19], [21], [25]], dtype=FLOAT64)
        self.assertEqual(result, expected)


class Test_jacobian_t(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, dimension):
        from bezier import _surface_helpers

        return _surface_helpers.jacobian_t(nodes, degree, dimension)

    def test_linear(self):
        nodes = np.asfortranarray([[0.0], [np.nan], [1.0]])
        result = self._call_function_under_test(nodes, 1, 1)
        expected = np.asfortranarray([[1.0]])
        self.assertEqual(result, expected)

    def test_quadratic(self):
        nodes = np.asfortranarray([
            [4.0, -2.0],
            [0.0, 1.0],
            [np.nan, np.nan],
            [5.0, 7.0],
            [-1.0, 6.0],
            [1.0, 12.0],
        ])
        result = self._call_function_under_test(nodes, 2, 2)
        expected = 2.0 * np.asfortranarray([
            [1.0, 9.0],
            [-1.0, 5.0],
            [-4.0, 5.0],
        ])
        self.assertEqual(result, expected)

    def test_cubic(self):
        nodes = np.arange(10, dtype=FLOAT64)[:, np.newaxis]**2
        result = self._call_function_under_test(nodes, 3, 1)
        expected = 3 * np.asfortranarray(
            [[16], [24], [32], [33], [39], [32]], dtype=FLOAT64)
        self.assertEqual(result, expected)

    def test_quartic(self):
        nodes = np.arange(15, dtype=FLOAT64)[:, np.newaxis]**2
        result = self._call_function_under_test(nodes, 4, 1)
        expected = 4 * np.asfortranarray([
            [25], [35], [45], [55], [56],
            [64], [72], [63], [69], [52],
        ], dtype=FLOAT64)
        self.assertEqual(result, expected)


class Test__jacobian_both(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, dimension):
        from bezier import _surface_helpers

        return _surface_helpers._jacobian_both(nodes, degree, dimension)

    def test_linear(self):
        # B(s, t) = -2s + 2t + 3
        nodes = np.asfortranarray([
            [3.0],
            [1.0],
            [5.0],
        ])
        result = self._call_function_under_test(nodes, 1, 1)

        # B_s = -2
        # B_t = 2
        expected = np.asfortranarray([
            [-2.0, 2.0],
        ])
        self.assertEqual(result, expected)

    def test_quadratic(self):
        # B(s, t) = [
        #     4 s t - 2 s + 5 t^2 - 6 t + 3,
        #     -s (s - 2 t),
        #     8 s^2 - 10 s t - 2 s - 13 t^2 + 12 t + 1,
        # ]
        #
        nodes = np.asfortranarray([
            [3.0, 0.0, 1.0],
            [2.0, 0.0, 0.0],
            [1.0, -1.0, 7.0],
            [0.0, 0.0, 7.0],
            [1.0, 1.0, 1.0],
            [2.0, 0.0, 0.0],
        ])
        result = self._call_function_under_test(nodes, 2, 3)

        # B_s = [
        #     4 t - 2,
        #     -2 s + 2 t,
        #     16 s - 10 t - 2,
        # ]
        # B_t = [
        #     4 s + 10 t - 6,
        #     2 s,
        #    -10 s - 26 t + 12,
        # ]
        expected = np.asfortranarray([
            [-2.0, 0.0, -2.0, -6.0, 0.0, 12.0],
            [-2.0, -2.0, 14.0, -2.0, 2.0, 2.0],
            [2.0, 2.0, -12.0, 4.0, 0.0, -14.0],
        ])
        self.assertEqual(result, expected)

    def test_cubic(self):
        # B(s, t) = [
        #     -2s^3 + 9s^2t + 12st^2 - 12st + 3s - 2t^3 + 6t,
        #     (-10s^3 - 30s^2t + 30s^2 - 36st^2 + 42st -
        #          18s - 15t^3 + 30t^2 - 18t + 7),
        # ]
        nodes = np.asfortranarray([
            [0.0, 7.0],
            [1.0, 1.0],
            [2.0, 5.0],
            [1.0, 9.0],
            [2.0, 1.0],
            [1.0, 2.0],
            [3.0, 3.0],
            [4.0, 5.0],
            [5.0, 1.0],
            [4.0, 4.0],
        ])
        result = self._call_function_under_test(nodes, 3, 2)

        # B_s = [
        #     -6s^2 + 18st + 12t^2 - 12t + 3,
        #     -30s^2 - 60st + 60s - 36t^2 + 42t - 18,
        # ]
        # B_t = [
        #     9s^2 + 24st - 12s - 6t^2 + 6,
        #     -30s^2 - 72st + 42s - 45t^2 + 60t - 18,
        # ]
        expected = np.asfortranarray([
            [3.0, -18.0, 6.0, -18.0],
            [3.0, 12.0, 0.0, 3.0],
            [-3.0, 12.0, 3.0, -6.0],
            [-3.0, 3.0, 6.0, 12.0],
            [6.0, 3.0, 12.0, -3.0],
            [3.0, -12.0, 0.0, -3.0],
        ])
        self.assertEqual(result, expected)


@utils.needs_surface_speedup
class Test_speedup_jacobian_both(Test__jacobian_both):

    @staticmethod
    def _call_function_under_test(nodes, degree, dimension):
        from bezier import _surface_speedup

        return _surface_speedup.jacobian_both(nodes, degree, dimension)


class Test__jacobian_det(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, st_vals):
        from bezier import _surface_helpers

        return _surface_helpers._jacobian_det(nodes, degree, st_vals)

    def test_linear(self):
        import bezier

        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 2.0],
        ])
        degree = 1
        surface = bezier.Surface(nodes, degree=degree, _copy=False)
        self.assertTrue(surface.is_valid)
        st_vals = np.asfortranarray(np.random.random((13, 2)))
        result = self._call_function_under_test(nodes, degree, st_vals)
        expected = 2.0 * np.ones(13, order='F')
        self.assertEqual(result, expected)

    def test_nonlinear(self):
        import bezier

        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.0, 0.5],
            [1.0, 1.0],
            [0.0, 1.0],
        ])
        degree = 2
        surface = bezier.Surface(nodes, degree=degree, _copy=False)
        self.assertTrue(surface.is_valid)
        # B(s, t) = [s(t + 1), t(s + 1)]
        st_vals = np.asfortranarray([
            [0.125, 0.125],
            [0.5, 0.375],
            [0.25, 0.75],
            [1.0, 0.0],
        ])
        result = self._call_function_under_test(nodes, degree, st_vals)
        # det(DB) = s + t + 1
        expected = np.asfortranarray([1.25, 1.875, 2.0, 2.0])
        self.assertEqual(result, expected)


@utils.needs_surface_speedup
class Test_speedup_jacobian_det(Test__jacobian_det):

    @staticmethod
    def _call_function_under_test(nodes, degree, st_vals):
        from bezier import _surface_speedup

        return _surface_speedup.jacobian_det(nodes, degree, st_vals)


class Test_classify_intersection(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection, edge_nodes1, edge_nodes2):
        from bezier import _surface_helpers

        return _surface_helpers.classify_intersection(
            intersection, edge_nodes1, edge_nodes2)

    def test_simple(self):
        first = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        edge_nodes1 = (None, None, first)
        second = np.asfortranarray([
            [0.25, 0.0],
            [0.75, 1.0],
        ])
        edge_nodes2 = (second, None, None)
        intersection = make_intersect(2, 0.5, 0, 0.5)
        result = self._call_function_under_test(
            intersection, edge_nodes1, edge_nodes2)
        self.assertIs(result, get_enum('SECOND'))

        # Swap and classify.
        intersection = make_intersect(0, 0.5, 2, 0.5)
        result = self._call_function_under_test(
            intersection, edge_nodes2, edge_nodes1)
        self.assertIs(result, get_enum('FIRST'))

    def test_corner_end(self):
        intersection = make_intersect(2, 1.0, 2, 0.5)
        with self.assertRaises(ValueError):
            self._call_function_under_test(intersection, (), ())

    def test_corner_start(self):
        import bezier

        surface = bezier.Surface.from_nodes(np.asfortranarray([
            [1.0, 1.0],
            [0.0, 0.0],
            [2.0, 0.0],
        ]))
        edge_nodes1 = tuple(edge._nodes for edge in surface.edges)
        second = np.asfortranarray([
            [1.0, 0.0],
            [1.0, 2.0],
        ])
        edge_nodes2 = (None, None, second)
        intersection = make_intersect(0, 0.0, 2, 0.5)
        result = self._call_function_under_test(
            intersection, edge_nodes1, edge_nodes2)
        self.assertIs(result, get_enum('FIRST'))

    def test_tangent(self):
        first = np.asfortranarray([
            [0.0, 0.0],
            [1.5, 1.0],
            [3.0, 0.0],
        ])
        edge_nodes1 = (first, None, None)
        second = np.asfortranarray([
            [1.0, 0.0],
            [1.5, 1.0],
            [2.0, 0.0],
        ])
        edge_nodes2 = (None, second, None)
        intersection = make_intersect(0, 0.5, 1, 0.5)
        result = self._call_function_under_test(
            intersection, edge_nodes1, edge_nodes2)
        self.assertIs(result, get_enum('TANGENT_FIRST'))

    def test_ignored_corner(self):
        import bezier

        surface1 = bezier.Surface(UNIT_TRIANGLE, 1)
        edge_nodes1 = tuple(edge._nodes for edge in surface1.edges)
        surface2 = bezier.Surface.from_nodes(np.asfortranarray([
            [0.0, 0.0],
            [-1.0, 0.0],
            [0.0, -1.0],
        ]))
        edge_nodes2 = tuple(edge._nodes for edge in surface2.edges)

        intersection = make_intersect(0, 0.0, 0, 0.0)
        result = self._call_function_under_test(
            intersection, edge_nodes1, edge_nodes2)
        self.assertIs(result, get_enum('IGNORED_CORNER'))


class Test_classify_tangent_intersection(unittest.TestCase):

    QUADRATIC1 = np.asfortranarray([
        [1.0, 0.0],
        [1.5, 1.0],
        [2.0, 0.0],
    ])
    QUADRATIC2 = np.asfortranarray([
        [0.0, 0.0],
        [1.5, 1.0],
        [3.0, 0.0],
    ])
    QUADRATIC3 = np.asfortranarray([
        [1.0, 1.0],
        [1.5, 0.0],
        [2.0, 1.0],
    ])

    @staticmethod
    def _call_function_under_test(
            intersection, nodes1, tangent1, nodes2, tangent2):
        from bezier import _surface_helpers

        return _surface_helpers.classify_tangent_intersection(
            intersection, nodes1, tangent1, nodes2, tangent2)

    def _call_helper(self, intersection, first, second):
        from bezier import _curve_helpers

        nodes1 = first._nodes
        tangent1 = _curve_helpers.evaluate_hodograph(intersection.s, nodes1)
        nodes2 = second._nodes
        tangent2 = _curve_helpers.evaluate_hodograph(intersection.t, nodes2)

        return self._call_function_under_test(
            intersection, nodes1, tangent1, nodes2, tangent2)

    def test_first_curvature(self):
        import bezier

        first = bezier.Curve(self.QUADRATIC1[::-1, :], 2)
        second = bezier.Curve(self.QUADRATIC2[::-1, :], 2)
        intersection = make_intersect(1, 0.5, 1, 0.5)

        result = self._call_helper(intersection, first, second)
        self.assertIs(result, get_enum('TANGENT_FIRST'))

    def test_second_curvature(self):
        import bezier

        first = bezier.Curve(self.QUADRATIC1, 2)
        second = bezier.Curve(self.QUADRATIC2, 2)
        intersection = make_intersect(1, 0.5, 0, 0.5)

        result = self._call_helper(intersection, first, second)
        self.assertIs(result, get_enum('TANGENT_SECOND'))

    def test_same_direction_same_curvature(self):
        import bezier

        first = bezier.Curve.from_nodes(np.asfortranarray([
            [1.0, 0.25],
            [-0.5, -0.25],
            [0.0, 0.25],
        ]))
        second = bezier.Curve.from_nodes(np.asfortranarray([
            [0.75, 0.25],
            [-0.25, -0.25],
            [-0.25, 0.25],
        ]))
        intersection = make_intersect(0, 0.5, 0, 0.5)
        with self.assertRaises(NotImplementedError):
            self._call_helper(intersection, first, second)

    def test_opposed_same_curvature(self):
        import bezier

        first = bezier.Curve.from_nodes(np.asfortranarray([
            [0.0, 0.25],
            [-0.5, -0.25],
            [1.0, 0.25],
        ]))
        second = bezier.Curve.from_nodes(np.asfortranarray([
            [0.75, 0.25],
            [-0.25, -0.25],
            [-0.25, 0.25],
        ]))
        intersection = make_intersect(1, 0.5, 2, 0.5)
        with self.assertRaises(NotImplementedError):
            self._call_helper(intersection, first, second)

    def test_opposed_same_sign_curvature_no_overlap(self):
        import bezier

        first = bezier.Curve(self.QUADRATIC1[::-1, :], 2)
        second = bezier.Curve(self.QUADRATIC3, 2)
        intersection = make_intersect(2, 0.5, 1, 0.5)

        result = self._call_helper(intersection, first, second)
        self.assertIs(result, get_enum('OPPOSED'))

    def test_opposed_same_sign_curvature_with_overlap(self):
        import bezier

        first = bezier.Curve(self.QUADRATIC1, 2)
        second = bezier.Curve(self.QUADRATIC3[::-1, :], 2)
        intersection = make_intersect(1, 0.5, 1, 0.5)

        with self.assertRaises(NotImplementedError):
            self._call_helper(intersection, first, second)

    def test_opposed_opp_sign_curvature_no_overlap(self):
        import bezier

        first = bezier.Curve(self.QUADRATIC1[::-1, :], 2)
        second = bezier.Curve(self.QUADRATIC2, 2)
        intersection = make_intersect(1, 0.5, 2, 0.5)

        result = self._call_helper(intersection, first, second)
        self.assertIs(result, get_enum('OPPOSED'))

    def test_opposed_opp_sign_curvature_with_overlap(self):
        import bezier

        first = bezier.Curve(self.QUADRATIC1, 2)
        second = bezier.Curve(self.QUADRATIC2[::-1, :], 2)
        intersection = make_intersect(1, 0.5, 0, 0.5)

        with self.assertRaises(NotImplementedError):
            self._call_helper(intersection, first, second)


class Test_ignored_edge_corner(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(
            edge_tangent, corner_tangent, corner_previous_edge):
        from bezier import _surface_helpers

        return _surface_helpers.ignored_edge_corner(
            edge_tangent, corner_tangent, corner_previous_edge)

    def test_first_across(self):
        edge_tangent = np.asfortranarray([[1.0, 0.0]])
        corner_tangent = np.asfortranarray([[0.0, 1.0]])
        self.assertFalse(
            self._call_function_under_test(edge_tangent, corner_tangent, None))

    def test_outside(self):
        edge_tangent = np.asfortranarray([[-1.0, 1.0]])
        corner_tangent = np.asfortranarray([[0.5, 0.5]])
        corner_previous_edge = np.asfortranarray([
            [0.5, 2.0],
            [0.5, 0.5],
        ])
        result = self._call_function_under_test(
            edge_tangent, corner_tangent, corner_previous_edge)
        self.assertTrue(result)

    def test_straddle(self):
        edge_tangent = np.asfortranarray([[1.0, 0.0]])
        corner_tangent = np.asfortranarray([[1.0, -1.0]])
        corner_previous_edge = np.asfortranarray([
            [1.0, 1.0],
            [0.5, 0.0],
        ])
        result = self._call_function_under_test(
            edge_tangent, corner_tangent, corner_previous_edge)
        self.assertFalse(result)


class Test_ignored_double_corner(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(
            intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2):
        from bezier import _surface_helpers

        return _surface_helpers.ignored_double_corner(
            intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2)

    def test_ignored(self):
        import bezier

        surface1 = bezier.Surface.from_nodes(np.asfortranarray([
            [1.0, 0.0],
            [1.5, 0.25],
            [0.5, 1.0],
        ]))
        edge_nodes1 = tuple(edge._nodes for edge in surface1.edges)
        surface2 = bezier.Surface(UNIT_TRIANGLE, 1)
        edge_nodes2 = tuple(edge._nodes for edge in surface2.edges)
        intersection = make_intersect(0, 0.0, 1, 0.0)
        tangent_s = np.asfortranarray([[0.5, 0.25]])
        tangent_t = np.asfortranarray([[-1.0, 1.0]])

        result = self._call_function_under_test(
            intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2)
        self.assertTrue(result)

    def test_overlap_first(self):
        import bezier

        surface1 = bezier.Surface(UNIT_TRIANGLE, 1)
        edge_nodes1 = tuple(edge._nodes for edge in surface1.edges)
        surface2 = bezier.Surface.from_nodes(np.asfortranarray([
            [1.0, 0.0],
            [1.0, 1.0],
            [0.5, 0.25],
        ]))
        edge_nodes2 = tuple(edge._nodes for edge in surface2.edges)
        intersection = make_intersect(1, 0.0, 0, 0.0)
        tangent_s = np.asfortranarray([[-1.0, 1.0]])
        tangent_t = np.asfortranarray([[0.0, 1.0]])

        result = self._call_function_under_test(
            intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2)
        self.assertFalse(result)

    def test_overlap_second(self):
        import bezier

        surface1 = bezier.Surface.from_nodes(np.asfortranarray([
            [1.0, 0.0],
            [1.0, 1.0],
            [0.5, 0.25],
        ]))
        edge_nodes1 = tuple(edge._nodes for edge in surface1.edges)
        surface2 = bezier.Surface(UNIT_TRIANGLE, 1)
        edge_nodes2 = tuple(edge._nodes for edge in surface2.edges)
        intersection = make_intersect(0, 0.0, 1, 0.0)
        tangent_s = np.asfortranarray([[0.0, 1.0]])
        tangent_t = np.asfortranarray([[-1.0, 1.0]])

        result = self._call_function_under_test(
            intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2)
        self.assertFalse(result)

    def test_segment_contained(self):
        import bezier

        surface1 = bezier.Surface.from_nodes(np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.5],
            [0.5, 1.0],
        ]))
        edge_nodes1 = tuple(edge._nodes for edge in surface1.edges)
        surface2 = bezier.Surface(UNIT_TRIANGLE, 1)
        edge_nodes2 = tuple(edge._nodes for edge in surface2.edges)
        intersection = make_intersect(0, 0.0, 0, 0.0)
        tangent_s = np.asfortranarray([[1.0, 0.5]])
        tangent_t = np.asfortranarray([[1.0, 0.0]])

        result = self._call_function_under_test(
            intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2)
        self.assertFalse(result)


class Test_ignored_corner(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(
            intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2):
        from bezier import _surface_helpers

        return _surface_helpers.ignored_corner(
            intersection, tangent_s, tangent_t, edge_nodes1, edge_nodes2)

    def test_not_corner(self):
        intersection = make_intersect(0, 0.5, 2, 0.5)
        result = self._call_function_under_test(
            intersection, None, None, (), ())
        self.assertFalse(result)

    def test_s_corner(self):
        import bezier

        surface = bezier.Surface(UNIT_TRIANGLE, degree=1, _copy=False)
        edge_nodes1 = tuple(edge._nodes for edge in surface.edges)
        edge_nodes2 = ()
        intersection = make_intersect(2, 0.0, None, 0.5)

        patch = unittest.mock.patch(
            'bezier._surface_helpers.ignored_edge_corner',
            return_value=unittest.mock.sentinel.edge_result)
        with patch as mocked:
            result = self._call_function_under_test(
                intersection, unittest.mock.sentinel.tangent_s,
                unittest.mock.sentinel.tangent_t, edge_nodes1, edge_nodes2)

        self.assertIs(result, unittest.mock.sentinel.edge_result)
        self.assertEqual(mocked.call_count, 1)
        call = mocked.mock_calls[0]
        _, positional, keyword = call
        self.assertEqual(keyword, {})
        self.assertEqual(len(positional), 3)
        self.assertIs(positional[0], unittest.mock.sentinel.tangent_t)
        self.assertIs(positional[1], unittest.mock.sentinel.tangent_s)
        _, previous_edge, _ = edge_nodes1
        self.assertEqual(positional[2], previous_edge)

    def test_t_corner(self):
        import bezier

        surface = bezier.Surface(UNIT_TRIANGLE, degree=1, _copy=False)
        edge_nodes1 = ()
        edge_nodes2 = tuple(edge._nodes for edge in surface.edges)
        intersection = make_intersect(None, 0.5, 1, 0.0)

        patch = unittest.mock.patch(
            'bezier._surface_helpers.ignored_edge_corner',
            return_value=unittest.mock.sentinel.edge_result)
        with patch as mocked:
            result = self._call_function_under_test(
                intersection, unittest.mock.sentinel.tangent_s,
                unittest.mock.sentinel.tangent_t, edge_nodes1, edge_nodes2)

        self.assertIs(result, unittest.mock.sentinel.edge_result)
        self.assertEqual(mocked.call_count, 1)
        call = mocked.mock_calls[0]
        _, positional, keyword = call
        self.assertEqual(keyword, {})
        self.assertEqual(len(positional), 3)
        self.assertIs(positional[0], unittest.mock.sentinel.tangent_s)
        self.assertIs(positional[1], unittest.mock.sentinel.tangent_t)
        previous_edge, _, _ = edge_nodes2
        self.assertEqual(positional[2], previous_edge)

    def test_double_corner(self):
        intersection = make_intersect(0, 0.0, 2, 0.0)

        patch = unittest.mock.patch(
            'bezier._surface_helpers.ignored_double_corner',
            return_value=unittest.mock.sentinel.double_result)
        with patch as mocked:
            result = self._call_function_under_test(
                intersection, unittest.mock.sentinel.tangent_s,
                unittest.mock.sentinel.tangent_t, (), ())

        self.assertIs(result, unittest.mock.sentinel.double_result)
        mocked.assert_called_once_with(
            intersection, unittest.mock.sentinel.tangent_s,
            unittest.mock.sentinel.tangent_t, (), ())


class Test_handle_ends(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(index1, s, index2, t):
        from bezier import _surface_helpers

        return _surface_helpers.handle_ends(index1, s, index2, t)

    def test_neither(self):
        edge_end, intersection_args = self._call_function_under_test(
            0, 0.5, 1, 0.5)
        self.assertFalse(edge_end)
        self.assertEqual(intersection_args, (0, 0.5, 1, 0.5))

    def test_s(self):
        edge_end, intersection_args = self._call_function_under_test(
            2, 1.0, 1, 0.25)
        self.assertTrue(edge_end)
        self.assertEqual(intersection_args, (0, 0.0, 1, 0.25))

    def test_t(self):
        edge_end, intersection_args = self._call_function_under_test(
            2, 0.75, 0, 1.0)
        self.assertTrue(edge_end)
        self.assertEqual(intersection_args, (2, 0.75, 1, 0.0))

    def test_both(self):
        edge_end, intersection_args = self._call_function_under_test(
            1, 1.0, 1, 1.0)
        self.assertTrue(edge_end)
        self.assertEqual(intersection_args, (2, 0.0, 2, 0.0))


class Test_to_front(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection, intersections, unused):
        from bezier import _surface_helpers

        return _surface_helpers.to_front(intersection, intersections, unused)

    def test_no_change(self):
        intersection = make_intersect(
            1, 0.5, 2, 0.5, interior_curve=get_enum('FIRST'))
        result = self._call_function_under_test(intersection, [], [])
        self.assertIs(result, intersection)

    def test_move_s(self):
        from bezier import _intersection_helpers

        enum_val = get_enum('FIRST')
        intersection = make_intersect(
            2, 1.0, None, None, interior_curve=enum_val)

        result = self._call_function_under_test(intersection, [], [])
        self.assertIsNot(result, intersection)
        self.assertIsInstance(result, _intersection_helpers.Intersection)
        self.assertEqual(result.s, 0.0)
        self.assertEqual(result.index_first, 0)
        self.assertIsNone(result.t)
        self.assertIsNone(result.index_second)
        self.assertEqual(result.interior_curve, enum_val)

    def test_move_s_to_existing(self):
        intersection = make_intersect(
            1, 1.0, None, None, interior_curve=get_enum('FIRST'))
        existing_int = make_intersect(
            2, 0.0, 0, 0.5, interior_curve=get_enum('FIRST'))

        result = self._call_function_under_test(
            intersection, [existing_int], [])
        self.assertIs(result, existing_int)

    def test_move_s_to_existing_and_remove(self):
        intersection = make_intersect(
            1, 1.0, None, None, interior_curve=get_enum('FIRST'))
        existing_int = make_intersect(
            2, 0.0, 0, 0.5, interior_curve=get_enum('FIRST'))
        unused = [existing_int]

        result = self._call_function_under_test(
            intersection, [existing_int], unused)
        self.assertIs(result, existing_int)
        self.assertEqual(unused, [])

    def test_move_t(self):
        from bezier import _intersection_helpers

        enum_val = get_enum('SECOND')
        intersection = make_intersect(
            None, None, 1, 1.0, interior_curve=enum_val)

        result = self._call_function_under_test(intersection, [], [])
        self.assertIsNot(result, intersection)
        self.assertIsInstance(result, _intersection_helpers.Intersection)
        self.assertIsNone(result.s)
        self.assertIsNone(result.index_first)
        self.assertEqual(result.t, 0.0)
        self.assertEqual(result.index_second, 2)
        self.assertEqual(result.interior_curve, enum_val)

    def test_move_t_to_existing(self):
        intersection = make_intersect(
            None, None, 0, 1.0, interior_curve=get_enum('SECOND'))
        existing_int = make_intersect(
            2, 0.5, 1, 0.0, interior_curve=get_enum('SECOND'))

        result = self._call_function_under_test(
            intersection, [existing_int], [])
        self.assertIs(result, existing_int)

    def test_move_t_to_existing_and_remove(self):
        intersection = make_intersect(
            None, None, 0, 1.0, interior_curve=get_enum('SECOND'))
        existing_int = make_intersect(
            2, 0.5, 1, 0.0, interior_curve=get_enum('SECOND'))
        unused = [existing_int]

        result = self._call_function_under_test(
            intersection, [existing_int], unused)
        self.assertIs(result, existing_int)
        self.assertEqual(unused, [])


class Test_get_next_first(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection, intersections):
        from bezier import _surface_helpers

        return _surface_helpers.get_next_first(intersection, intersections)

    def test_move_to_corner(self):
        from bezier import _intersection_helpers

        intersection = make_intersect(1, 0.25, None, None)
        result = self._call_function_under_test(intersection, [])
        self.assertIsInstance(result, _intersection_helpers.Intersection)
        self.assertEqual(result.s, 1.0)
        self.assertEqual(result.index_first, 1)
        self.assertIsNone(result.t)
        self.assertIsNone(result.index_second)
        self.assertIs(result.interior_curve, get_enum('FIRST'))

    def test_move_to_existing(self):
        intersection = make_intersect(2, 0.25, None, None)
        intersections = [
            # An "acceptable" intersection that will be overtaken by the
            # next since 0.25 < 0.5 < 0.875.
            make_intersect(
                2, 0.875, None, None, interior_curve=get_enum('SECOND')),
            make_intersect(
                2, 0.5, None, None, interior_curve=get_enum('FIRST')),
            # On a different curve.
            make_intersect(1, None, None, None),
            # Same curve, but before.
            make_intersect(
                2, 0.125, None, None, interior_curve=get_enum('FIRST')),
            # Past the already accepted intersection.
            make_intersect(
                2, 0.625, None, None, interior_curve=get_enum('FIRST')),
        ]
        result = self._call_function_under_test(intersection, intersections)
        self.assertIs(result, intersections[1])


class Test_get_next_second(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection, intersections):
        from bezier import _surface_helpers

        return _surface_helpers.get_next_second(intersection, intersections)

    def test_move_to_corner(self):
        from bezier import _intersection_helpers

        intersection = make_intersect(None, None, 1, 0.625)
        result = self._call_function_under_test(intersection, [])
        self.assertIsInstance(result, _intersection_helpers.Intersection)
        self.assertIsNone(result.s)
        self.assertIsNone(result.index_first)
        self.assertEqual(result.t, 1.0)
        self.assertEqual(result.index_second, 1)
        self.assertIs(result.interior_curve, get_enum('SECOND'))

    def test_move_to_existing(self):
        intersection = make_intersect(None, None, 1, 0.125)
        intersections = [
            # An "acceptable" intersection that will be overtaken by the
            # next since 0.125 < 0.625 < 0.75.
            make_intersect(
                None, None, 1, 0.75, interior_curve=get_enum('FIRST')),
            make_intersect(
                None, None, 1, 0.625, interior_curve=get_enum('SECOND')),
            # On a different curve.
            make_intersect(None, None, 2, None),
            # Same curve, but before.
            make_intersect(
                None, None, 1, 0.0625, interior_curve=get_enum('FIRST')),
            # Past the already accepted intersection.
            make_intersect(
                None, None, 1, 0.6875, interior_curve=get_enum('SECOND')),
        ]
        result = self._call_function_under_test(intersection, intersections)
        self.assertIs(result, intersections[1])


class Test_get_next(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection, intersections, unused):
        from bezier import _surface_helpers

        return _surface_helpers.get_next(
            intersection, intersections, unused)

    def test_remove_from_unused(self):
        # Also tests branch through "first".
        return_value = make_intersect(
            1, None, 0, None, interior_curve=get_enum('SECOND'))
        unused = [return_value]
        intersection = make_intersect(
            2, None, 2, None, interior_curve=get_enum('FIRST'))

        patch = unittest.mock.patch(
            'bezier._surface_helpers.get_next_first',
            return_value=return_value)
        with patch as mocked:
            result = self._call_function_under_test(
                intersection, unittest.mock.sentinel.intersections, unused)

        self.assertIs(result, return_value)
        self.assertEqual(unused, [])
        mocked.assert_called_once_with(
            intersection, unittest.mock.sentinel.intersections)

    def test_second(self):
        return_value = make_intersect(
            0, None, 2, None, interior_curve=get_enum('FIRST'))
        intersection = make_intersect(
            0, None, 2, None, interior_curve=get_enum('SECOND'))

        patch = unittest.mock.patch(
            'bezier._surface_helpers.get_next_second',
            return_value=return_value)
        with patch as mocked:
            result = self._call_function_under_test(
                intersection, unittest.mock.sentinel.intersections, [])

        self.assertIs(result, return_value)
        mocked.assert_called_once_with(
            intersection, unittest.mock.sentinel.intersections)

    def test_invalid_classification(self):
        intersection = make_intersect(
            1, None, 1, None, interior_curve=get_enum('OPPOSED'))
        with self.assertRaises(ValueError):
            self._call_function_under_test(intersection, [], [])


class Test_ends_to_curve(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(start_node, end_node):
        from bezier import _surface_helpers

        return _surface_helpers.ends_to_curve(start_node, end_node)

    def test_bad_classification(self):
        start_node = make_intersect(0, 0.5, 0, 0.5)
        end_node = make_intersect(0, 0.75, 0, 0.75)
        with self.assertRaises(ValueError):
            self._call_function_under_test(start_node, end_node)

    def _on_different_curves(self, interior_curve):
        start_node = make_intersect(
            0, 0.5, 2, 0.5, interior_curve=interior_curve)
        end_node = make_intersect(1, 0.5, 1, 0.5)
        with self.assertRaises(ValueError):
            self._call_function_under_test(start_node, end_node)

    def test_first_on_different_curves(self):
        self._on_different_curves(get_enum('FIRST'))

    def test_second_on_different_curves(self):
        self._on_different_curves(get_enum('SECOND'))

    def test_first(self):
        start_node = make_intersect(
            0, 0.5, None, None, interior_curve=get_enum('FIRST'))
        end_node = make_intersect(0, 0.75, None, None)
        result = self._call_function_under_test(start_node, end_node)
        self.assertEqual(result, (0, 0.5, 0.75))

    def test_second(self):
        start_node = make_intersect(
            None, None, 2, 0.125, interior_curve=get_enum('SECOND'))
        end_node = make_intersect(None, None, 2, 0.25)
        result = self._call_function_under_test(start_node, end_node)
        self.assertEqual(result, (5, 0.125, 0.25))


class Test_no_intersections(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes1, degree1, nodes2, degree2):
        from bezier import _surface_helpers

        return _surface_helpers.no_intersections(
            nodes1, degree1, nodes2, degree2)

    def test_disjoint(self):
        nodes2 = UNIT_TRIANGLE + np.asfortranarray([[5.0, 0.0]])
        edges_info, contained = self._call_function_under_test(
            UNIT_TRIANGLE, 1, nodes2, 1)
        self.assertEqual(edges_info, [])
        self.assertIsNone(contained)

    def test_first_contained(self):
        nodes2 = 4.0 * UNIT_TRIANGLE - np.asfortranarray([[1.0, 1.0]])
        edges_info, contained = self._call_function_under_test(
            UNIT_TRIANGLE, 1, nodes2, 1)
        self.assertIsNone(edges_info)
        self.assertTrue(contained)

    def test_second_contained(self):
        nodes1 = 4.0 * UNIT_TRIANGLE - np.asfortranarray([[1.0, 1.0]])
        edges_info, contained = self._call_function_under_test(
            nodes1, 1, UNIT_TRIANGLE, 1)
        self.assertIsNone(edges_info)
        self.assertIsNotNone(contained)
        self.assertFalse(contained)


class Test_tangent_only_intersections(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(all_types):
        from bezier import _surface_helpers

        return _surface_helpers.tangent_only_intersections(all_types)

    def test_too_few_types(self):
        with self.assertRaises(ValueError):
            self._call_function_under_test(set())

    def test_too_many_types(self):
        all_types = set([get_enum('FIRST'), get_enum('SECOND')])
        with self.assertRaises(ValueError):
            self._call_function_under_test(all_types)

    def test_bad_types(self):
        all_types = set([get_enum('FIRST')])
        with self.assertRaises(ValueError):
            self._call_function_under_test(all_types)

    def test_ignored_types(self):
        all_types = set([get_enum('OPPOSED')])
        edges_info, contained = self._call_function_under_test(all_types)
        self.assertEqual(edges_info, [])
        self.assertIsNone(contained)

        all_types = set([get_enum('IGNORED_CORNER')])
        edges_info, contained = self._call_function_under_test(all_types)
        self.assertEqual(edges_info, [])
        self.assertIsNone(contained)

    def test_first(self):
        all_types = set([get_enum('TANGENT_FIRST')])
        edges_info, contained = self._call_function_under_test(all_types)
        self.assertIsNone(edges_info)
        self.assertTrue(contained)

    def test_second(self):
        all_types = set([get_enum('TANGENT_SECOND')])
        edges_info, contained = self._call_function_under_test(all_types)
        self.assertIsNone(edges_info)
        self.assertIsNotNone(contained)
        self.assertFalse(contained)


class Test_basic_interior_combine(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(intersections, **kwargs):
        from bezier import _surface_helpers

        return _surface_helpers.basic_interior_combine(
            intersections, **kwargs)

    def test_it(self):
        intersection1 = make_intersect(
            1, 0.25, 0, 0.4375, interior_curve=get_enum('FIRST'))
        intersection2 = make_intersect(
            1, 0.5, 2, 0.75, interior_curve=get_enum('SECOND'))

        edges_info, contained = self._call_function_under_test(
            [intersection1, intersection2])
        expected = [
            (
                (5, 0.75, 1.0),
                (3, 0.0, 0.4375),
                (1, 0.25, 0.5),
            ),
        ]
        self.assertEqual(edges_info, expected)
        self.assertIsNone(contained)

    def test_two_curved_polygons(self):
        intersections = [
            make_intersect(
                0, 0.625, 0, 0.25, interior_curve=get_enum('FIRST')),
            make_intersect(
                0, 0.375, 0, 0.75, interior_curve=get_enum('SECOND')),
            make_intersect(
                0, 0.3125, 1, 0.25, interior_curve=get_enum('FIRST')),
            make_intersect(
                0, 0.6875, 2, 0.75, interior_curve=get_enum('SECOND')),
        ]

        edges_info, contained = self._call_function_under_test(intersections)
        expected = [
            (
                (5, 0.75, 1.0),
                (3, 0.0, 0.25),
                (0, 0.625, 0.6875),
            ), (
                (0, 0.3125, 0.375),
                (3, 0.75, 1.0),
                (4, 0.0, 0.25),
            ),
        ]
        self.assertEqual(edges_info, expected)
        self.assertIsNone(contained)

    def test_first_contained(self):
        intersections = [
            make_intersect(
                2, 0.0, 1, 0.5, interior_curve=get_enum('FIRST')),
            make_intersect(
                0, 0.0, 2, 0.625, interior_curve=get_enum('FIRST')),
            make_intersect(
                1, 0.0, 0, 0.25, interior_curve=get_enum('FIRST')),
        ]
        edges_info, contained = self._call_function_under_test(intersections)
        self.assertIsNone(edges_info)
        self.assertTrue(contained)

    def test_second_contained(self):
        intersections = [
            make_intersect(
                0, 0.875, 1, 0.0, interior_curve=get_enum('SECOND')),
            make_intersect(
                1, 0.75, 2, 0.0, interior_curve=get_enum('SECOND')),
            make_intersect(
                2, 0.875, 0, 0.0, interior_curve=get_enum('SECOND')),
        ]
        edges_info, contained = self._call_function_under_test(intersections)
        self.assertIsNone(edges_info)
        self.assertIsNotNone(contained)
        self.assertFalse(contained)

    def _too_many_edges_helper(self, to_front, get_next, **kwargs):
        start = make_intersect(
            1, 0.0, 1, 0.0, interior_curve=get_enum('SECOND'))
        with self.assertRaises(RuntimeError):
            self._call_function_under_test([start], **kwargs)

        max_edges = kwargs.pop('max_edges', 10)
        self.assertEqual(kwargs, {})

        self.assertEqual(to_front.call_count, max_edges)
        self.assertEqual(get_next.call_count, max_edges + 1)

    @unittest.mock.patch(
        'bezier._surface_helpers.get_next',
        return_value=unittest.mock.sentinel.next_)
    @unittest.mock.patch(
        'bezier._surface_helpers.to_front',
        return_value=unittest.mock.sentinel.front)
    def test_too_many_edges(self, to_front, get_next):
        self._too_many_edges_helper(to_front, get_next)

    @unittest.mock.patch(
        'bezier._surface_helpers.get_next',
        return_value=unittest.mock.sentinel.next_)
    @unittest.mock.patch(
        'bezier._surface_helpers.to_front',
        return_value=unittest.mock.sentinel.front)
    def test_too_many_edges_explicit_max(self, to_front, get_next):
        self._too_many_edges_helper(to_front, get_next, max_edges=3)


class Test_combine_intersections(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(
            intersections, nodes1, degree1, nodes2, degree2, all_types):
        from bezier import _surface_helpers

        return _surface_helpers.combine_intersections(
            intersections, nodes1, degree1, nodes2, degree2, all_types)

    def test_empty(self):
        nodes2 = np.asfortranarray([
            [-1.0, 0.0],
            [-1.0, 1.0],
            [-2.0, 1.0],
        ])
        edges_info, contained = self._call_function_under_test(
            [], UNIT_TRIANGLE, 1, nodes2, 1, None)
        self.assertEqual(edges_info, [])
        self.assertIsNone(contained)

    def test_basic(self):
        nodes2 = np.asfortranarray([
            [0.75, -0.25],
            [0.75, 0.75],
            [-0.25, 0.75],
        ])
        intersections = [
            make_intersect(
                0, 0.75, 0, 0.25, interior_curve=get_enum('SECOND')),
            make_intersect(
                0, 0.5, 2, 0.75, interior_curve=get_enum('FIRST')),
            make_intersect(
                1, 0.25, 0, 0.5, interior_curve=get_enum('FIRST')),
            make_intersect(
                1, 0.75, 1, 0.5, interior_curve=get_enum('SECOND')),
            make_intersect(
                2, 0.25, 1, 0.75, interior_curve=get_enum('FIRST')),
            make_intersect(
                2, 0.5, 2, 0.25, interior_curve=get_enum('SECOND')),
        ]
        edges_info, contained = self._call_function_under_test(
            intersections, UNIT_TRIANGLE, 1, nodes2, 1, None)
        expected = [
            (
                (5, 0.25, 0.75),
                (0, 0.5, 0.75),
                (3, 0.25, 0.5),
                (1, 0.25, 0.75),
                (4, 0.5, 0.75),
                (2, 0.25, 0.5),
            ),
        ]
        self.assertEqual(edges_info, expected)
        self.assertIsNone(contained)

    def test_tangent(self):
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [0.5, -0.5],
            [1.0, 0.0],
            [0.25, 0.5],
            [0.75, 0.5],
            [0.5, 1.0],
        ])
        nodes2 = np.asfortranarray([
            [-1.0, -0.25],
            [2.0, -0.25],
            [0.5, 1.5],
        ])

        # Surface1-Edge0(0.5) = Surface2-Edge2(0.5)
        enum_val = get_enum('TANGENT_FIRST')
        all_types = set([enum_val])
        edges_info, contained = self._call_function_under_test(
            [], nodes1, 2, nodes2, 1, all_types)
        self.assertIsNone(edges_info)
        self.assertTrue(contained)


class Test__evaluate_barycentric(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, lambda1, lambda2, lambda3):
        from bezier import _surface_helpers

        return _surface_helpers._evaluate_barycentric(
            nodes, degree, lambda1, lambda2, lambda3)

    def test_linear(self):
        lambda_vals = (0.25, 0.5, 0.25)
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.5],
            [0.0, 1.25],
        ])

        expected = np.asfortranarray([[0.5, 0.5625]])
        result = self._call_function_under_test(nodes, 1, *lambda_vals)
        self.assertEqual(result, expected)

    def test_quadratic(self):
        lambda_vals = (0.0, 0.25, 0.75)
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.5],
            [0.5, 1.25],
            [0.0, 1.25],
            [0.0, 0.5],
        ])

        expected = np.asfortranarray([[0.0625, 0.78125]])
        result = self._call_function_under_test(nodes, 2, *lambda_vals)
        self.assertEqual(result, expected)

    def test_quadratic_dimension3(self):
        lambda_vals = (0.125, 0.375, 0.5)
        nodes = np.asfortranarray([
            [0.0, 0.0, 1.0],
            [0.5, 0.0, 0.25],
            [1.0, 0.5, 0.0],
            [0.5, 1.25, 1.25],
            [0.0, 1.25, 0.5],
            [0.0, 0.5, -1.0],
        ])

        expected = np.asfortranarray([[0.25, 0.8203125, 0.1328125]])
        result = self._call_function_under_test(nodes, 2, *lambda_vals)
        self.assertEqual(result, expected)

    def test_cubic(self):
        lambda_vals = (0.125, 0.5, 0.375)
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.25, 0.0],
            [0.75, 0.25],
            [1.0, 0.0],
            [0.0, 0.25],
            [0.375, 0.25],
            [0.5, 0.25],
            [0.0, 0.5],
            [0.25, 0.75],
            [0.0, 1.0],
        ])

        expected = np.asfortranarray([[0.447265625, 0.37060546875]])
        result = self._call_function_under_test(nodes, 3, *lambda_vals)
        self.assertEqual(result, expected)

    def test_quartic(self):
        import math

        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(15, 2), seed=11112222, num_bits=8)

        lambda_vals = (0.125, 0.375, 0.5)
        index = 0
        expected = np.asfortranarray([[0.0, 0.0]])
        for k in range(4 + 1):
            for j in range(4 + 1 - k):
                i = 4 - j - k
                denom = (math.factorial(i) * math.factorial(j) *
                         math.factorial(k))
                coeff = 24 / denom
                expected += (
                    coeff * lambda_vals[0]**i * lambda_vals[1]**j *
                    lambda_vals[2]**k * nodes[[index], :])
                index += 1

        result = self._call_function_under_test(nodes, 4, *lambda_vals)
        self.assertEqual(result, expected)


@utils.needs_surface_speedup
class Test_speedup_evaluate_barycentric(Test__evaluate_barycentric):

    @staticmethod
    def _call_function_under_test(nodes, degree, lambda1, lambda2, lambda3):
        from bezier import _surface_speedup

        return _surface_speedup.evaluate_barycentric(
            nodes, degree, lambda1, lambda2, lambda3)


class Test__evaluate_barycentric_multi(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, param_vals, dimension):
        from bezier import _surface_helpers

        return _surface_helpers._evaluate_barycentric_multi(
            nodes, degree, param_vals, dimension)

    def test_basic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [2.0, 1.0],
            [-3.0, 2.0],
        ])
        expected = np.asfortranarray([
            [0.0, 0.0],
            [2.0, 1.0],
            [-0.5, 1.5],
        ])

        param_vals = np.asfortranarray([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.5, 0.5],
        ])
        result = self._call_function_under_test(nodes, 1, param_vals, 2)
        self.assertEqual(result, expected)

    def test_outside_domain(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [3.0, -1.0],
            [1.0, 0.0],
        ])
        expected = np.asfortranarray([
            [1.0, -0.25],
            [0.0, 1.0],
            [2.375, -0.75],
        ])

        param_vals = np.asfortranarray([
            [0.25, 0.25, 0.25],
            [-1.0, -1.0, 3.0],
            [0.125, 0.75, 0.125],
        ])
        result = self._call_function_under_test(nodes, 1, param_vals, 2)
        self.assertEqual(result, expected)


@utils.needs_surface_speedup
class Test_speedup_evaluate_barycentric_multi(
        Test__evaluate_barycentric_multi):

    @staticmethod
    def _call_function_under_test(nodes, degree, param_vals, dimension):
        from bezier import _surface_speedup

        return _surface_speedup.evaluate_barycentric_multi(
            nodes, degree, param_vals, dimension)


class Test__evaluate_cartesian_multi(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, param_vals, dimension):
        from bezier import _surface_helpers

        return _surface_helpers._evaluate_cartesian_multi(
            nodes, degree, param_vals, dimension)

    def test_basic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.75],
            [2.0, 1.0],
            [-1.5, 1.0],
            [-0.5, 1.5],
            [-3.0, 2.0],
        ])
        expected = np.asfortranarray([
            [-1.75, 1.75],
            [0.0, 0.0],
            [0.25, 1.0625],
            [-0.625, 1.046875],
        ])

        param_vals = np.asfortranarray([
            [0.25, 0.75],
            [0.0, 0.0],
            [0.5, 0.25],
            [0.25, 0.375],
        ])
        result = self._call_function_under_test(nodes, 2, param_vals, 2)
        self.assertEqual(result, expected)

    def test_outside_domain(self):
        nodes = np.asfortranarray([
            [0.0, 2.0],
            [1.0, 1.0],
            [1.0, 0.0],
        ])
        expected = np.asfortranarray([
            [0.5, 1.25],
            [2.0, -3.0],
            [0.875, 1.0],
        ])

        param_vals = np.asfortranarray([
            [0.25, 0.25],
            [-1.0, 3.0],
            [0.75, 0.125],
        ])
        result = self._call_function_under_test(nodes, 1, param_vals, 2)
        self.assertEqual(result, expected)


@utils.needs_surface_speedup
class Test_speedup_evaluate_cartesian_multi(Test__evaluate_cartesian_multi):

    @staticmethod
    def _call_function_under_test(nodes, degree, param_vals, dimension):
        from bezier import _surface_speedup

        return _surface_speedup.evaluate_cartesian_multi(
            nodes, degree, param_vals, dimension)


class Test__compute_edge_nodes(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _surface_helpers

        return _surface_helpers._compute_edge_nodes(nodes, degree)

    def _check(self, nodes, degree, expected1, expected2, expected3):
        nodes1, nodes2, nodes3 = self._call_function_under_test(nodes, degree)
        self.assertEqual(nodes1, expected1)
        self.assertEqual(nodes2, expected2)
        self.assertEqual(nodes3, expected3)

    def test_linear(self):
        nodes = np.asfortranarray([
            [1.0, 2.0],
            [4.0, 2.5],
            [0.0, 4.0],
        ])
        p100, p010, p001 = nodes

        self._check(
            nodes, 1,
            np.asfortranarray(np.vstack([p100, p010])),
            np.asfortranarray(np.vstack([p010, p001])),
            np.asfortranarray(np.vstack([p001, p100])))

    def test_quadratic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.25, 0.5],
            [2.0, 1.0],
            [-1.5, 0.75],
            [0.0, 2.0],
            [-3.0, 3.0],
        ])
        p200, p110, p020, p101, p011, p002 = nodes

        self._check(
            nodes, 2,
            np.asfortranarray(np.vstack([p200, p110, p020])),
            np.asfortranarray(np.vstack([p020, p011, p002])),
            np.asfortranarray(np.vstack([p002, p101, p200])))

    def test_cubic(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.328125, 0.1484375],
            [0.65625, 0.1484375],
            [1.0, 0.0],
            [0.1484375, 0.328125],
            [0.5, 0.5],
            [1.0, 0.53125],
            [0.1484375, 0.65625],
            [0.53125, 1.0],
            [0.0, 1.0],
        ])
        (p300, p210, p120, p030, p201,
         unused_p111, p021, p102, p012, p003) = nodes

        self._check(
            nodes, 3,
            np.asfortranarray(np.vstack([p300, p210, p120, p030])),
            np.asfortranarray(np.vstack([p030, p021, p012, p003])),
            np.asfortranarray(np.vstack([p003, p102, p201, p300])))


@utils.needs_surface_speedup
class Test_speedup_compute_edge_nodes(Test__compute_edge_nodes):

    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _surface_speedup

        return _surface_speedup.compute_edge_nodes(nodes, degree)


def make_intersect(*args, **kwargs):
    from bezier import _intersection_helpers

    return _intersection_helpers.Intersection(*args, **kwargs)


def get_enum(str_val):
    from bezier import _surface_helpers

    return _surface_helpers.IntersectionClassification[str_val]
