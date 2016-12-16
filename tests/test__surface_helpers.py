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

from tests import utils


UNIT_TRIANGLE = np.array([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
])


def make_intersect(left, s, right, t):
    from bezier import _intersection_helpers

    return _intersection_helpers.Intersection(left, s, right, t)


def get_enum(str_val):
    from bezier import _surface_helpers

    return _surface_helpers.IntersectionClassification[str_val]


class Test_polynomial_sign(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(poly_surface):
        from bezier import _surface_helpers

        return _surface_helpers.polynomial_sign(poly_surface)

    def _helper(self, bernstein, expected):
        import bezier

        poly_surface = bezier.Surface(bernstein)
        result = self._call_function_under_test(poly_surface)
        self.assertEqual(result, expected)

    def test_positive(self):
        # pylint: disable=no-member
        bernstein = np.array([[1.0, 2.0, 3.0, 4.0, 5.0, 6.0]]).T
        # pylint: enable=no-member
        self._helper(bernstein, 1)

    def test_negative(self):
        # pylint: disable=no-member
        bernstein = np.array([[-1.0, -2.0, -1.0]]).T
        # pylint: enable=no-member
        self._helper(bernstein, -1)

    def test_zero(self):
        bernstein = np.zeros((10, 1))
        self._helper(bernstein, 0)

    def test_mixed(self):
        # pylint: disable=no-member
        bernstein = np.array([[-1.0, 1.0, -1.0]]).T
        # pylint: enable=no-member
        self._helper(bernstein, 0)

    def test_max_iterations(self):
        # pylint: disable=no-member
        bernstein = np.array([[1.0, 2.0, 3.0]]).T
        # pylint: enable=no-member
        subs = 'bezier._surface_helpers._MAX_POLY_SUBDIVISIONS'
        with mock.patch(subs, new=1):
            self._helper(bernstein, 1)

    def test_no_conclusion(self):
        # pylint: disable=no-member
        bernstein = np.array([[-1.0, 1.0, 2.0]]).T
        # pylint: enable=no-member
        subs = 'bezier._surface_helpers._MAX_POLY_SUBDIVISIONS'
        with mock.patch(subs, new=0):
            with self.assertRaises(ValueError):
                self._helper(bernstein, None)


class Test__2x2_det(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(mat):
        from bezier import _surface_helpers

        return _surface_helpers._2x2_det(mat)

    def test_integers(self):
        mat = np.array([
            [1.0, 2.0],
            [3.0, 4.0],
        ])
        self.assertEqual(self._call_function_under_test(mat), -2.0)

    def test_better_than_numpy(self):
        mat = np.array([
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
        nodes = np.array([
            [1.0, 0.0],
            [0.0, 0.0],
            [1.0, 1.0],
            [0.0, 0.0],
            [0.0, 0.0],
            [0.0, 1.0],
        ])
        bernstein = self._call_function_under_test(nodes)
        self.assertEqual(bernstein.shape, (6, 1))
        expected = np.array([[0.0, 2.0, 0.0, -2.0, 2.0, 0.0]])
        expected = expected.T  # pylint: disable=no-member
        self.assertEqual(bernstein, expected)


class Test_cubic_jacobian_polynomial(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _surface_helpers

        return _surface_helpers.cubic_jacobian_polynomial(nodes)

    def test_it(self):
        # B(L1, L2, L3) = [L1^3 + L2^3, L2^3 + L3^3]
        nodes = np.array([
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
        expected = np.zeros(shape)
        expected[2, 0] = 1.5
        expected[9, 0] = -1.5
        expected[11, 0] = 1.5
        self.assertEqual(bernstein, expected)


class Test_de_casteljau_one_round(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, lambda1, lambda2, lambda3):
        from bezier import _surface_helpers

        return _surface_helpers.de_casteljau_one_round(
            nodes, degree, lambda1, lambda2, lambda3)

    def test_linear(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        s_val, t_val = 0.5, 0.375
        expected = np.array([
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

        q100 = (1.0 - s_val - t_val) * p200 + s_val * p110 + t_val * p101
        q010 = (1.0 - s_val - t_val) * p110 + s_val * p020 + t_val * p011
        q001 = (1.0 - s_val - t_val) * p101 + s_val * p011 + t_val * p002

        expected = np.vstack([q100, q010, q001])
        result = self._call_function_under_test(
            nodes, 2, 1.0 - s_val - t_val, s_val, t_val)
        self.assertEqual(result, expected)

    def test_cubic(self):
        nodes = np.array([
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
        transform = np.array([
            [lambda1, s_val, 0., 0., t_val, 0., 0., 0., 0., 0.],
            [0., lambda1, s_val, 0., 0., t_val, 0., 0., 0., 0.],
            [0., 0., lambda1, s_val, 0., 0., t_val, 0., 0., 0.],
            [0., 0., 0., 0., lambda1, s_val, 0., t_val, 0., 0.],
            [0., 0., 0., 0., 0., lambda1, s_val, 0., t_val, 0.],
            [0., 0., 0., 0., 0., 0., 0., lambda1, s_val, t_val],
        ])
        # pylint: disable=no-member
        expected = transform.dot(nodes)
        # pylint: enable=no-member
        result = self._call_function_under_test(
            nodes, 3, lambda1, s_val, t_val)
        self.assertEqual(result, expected)


class Test__make_transform(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(degree, weights_a, weights_b, weights_c):
        from bezier import _surface_helpers

        return _surface_helpers._make_transform(
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
        weights = np.array([
            [1.0, 0.0, 0.0],
            [0.5, 0.5, 0.0],
            [0.5, 0.0, 0.5],
        ])
        expected0 = weights[[0], :]
        expected1 = weights[[1], :]
        expected2 = weights[[2], :]
        self._helper(1, weights, expected0, expected1, expected2)

    def test_quadratic(self):
        weights = np.array([
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
        ])
        expected0 = np.array([
            [0.0, 0.5, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.0, 0.5, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.0, 0.5, 0.5],
        ])
        expected1 = np.array([
            [0.5, 0.0, 0.0, 0.5, 0.0, 0.0],
            [0.0, 0.5, 0.0, 0.0, 0.5, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.0, 0.5],
        ])
        expected2 = np.array([
            [0.5, 0.5, 0.0, 0.0, 0.0, 0.0],
            [0.0, 0.5, 0.5, 0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0, 0.5, 0.5, 0.0],
        ])
        self._helper(2, weights, expected0, expected1, expected2)


class Test__reduced_to_matrix(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(shape, degree, vals_by_weight):
        from bezier import _surface_helpers

        return _surface_helpers._reduced_to_matrix(
            shape, degree, vals_by_weight)

    def test_it(self):
        shape = (6, 2)
        degree = 2

        expected = np.array([
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


class Test_specialize_surface(utils.NumPyTestCase):

    WEIGHTS0 = (1.0, 0.0, 0.0)
    WEIGHTS1 = (0.5, 0.5, 0.0)
    WEIGHTS2 = (0.0, 1.0, 0.0)
    WEIGHTS3 = (0.5, 0.0, 0.5)
    WEIGHTS4 = (0.0, 0.5, 0.5)
    WEIGHTS5 = (0.0, 0.0, 1.0)

    @staticmethod
    def _call_function_under_test(nodes, degree,
                                  weights_a, weights_b, weights_c):
        from bezier import _surface_helpers

        return _surface_helpers.specialize_surface(
            nodes, degree, weights_a, weights_b, weights_c)

    def _helpers(self, degree, all_nodes, inds_a, inds_b, inds_c, inds_d):
        num_nodes = len(inds_a)
        id_mat = np.eye(num_nodes)

        expected_a = self._call_function_under_test(
            id_mat, degree,
            self.WEIGHTS0, self.WEIGHTS1, self.WEIGHTS3)
        expected_b = self._call_function_under_test(
            id_mat, degree,
            self.WEIGHTS4, self.WEIGHTS3, self.WEIGHTS1)
        expected_c = self._call_function_under_test(
            id_mat, degree,
            self.WEIGHTS1, self.WEIGHTS2, self.WEIGHTS4)
        expected_d = self._call_function_under_test(
            id_mat, degree,
            self.WEIGHTS3, self.WEIGHTS4, self.WEIGHTS5)

        self.assertEqual(all_nodes[inds_a, :], expected_a)
        self.assertEqual(all_nodes[inds_b, :], expected_b)
        self.assertEqual(all_nodes[inds_c, :], expected_c)
        self.assertEqual(all_nodes[inds_d, :], expected_d)

    def test_known_linear(self):
        from bezier import _surface_helpers

        all_nodes = _surface_helpers.LINEAR_SUBDIVIDE
        self._helpers(1, all_nodes,
                      (0, 1, 3), (4, 3, 1),
                      (1, 2, 4), (3, 4, 5))

    def test_known_quadratic(self):
        from bezier import _surface_helpers

        all_nodes = _surface_helpers.QUADRATIC_SUBDIVIDE
        self._helpers(2, all_nodes,
                      (0, 1, 2, 5, 6, 9),
                      (11, 10, 9, 7, 6, 2),
                      (2, 3, 4, 7, 8, 11),
                      (9, 10, 11, 12, 13, 14))

    def test_known_cubic(self):
        from bezier import _surface_helpers

        all_nodes = _surface_helpers.CUBIC_SUBDIVIDE
        self._helpers(3, all_nodes,
                      (0, 1, 2, 3, 7, 8, 9, 13, 14, 18),
                      (21, 20, 19, 18, 16, 15, 14, 10, 9, 3),
                      (3, 4, 5, 6, 10, 11, 12, 16, 17, 21),
                      (18, 19, 20, 21, 22, 23, 24, 25, 26, 27))


class Test__mean_centroid(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(candidates):
        from bezier import _surface_helpers

        return _surface_helpers._mean_centroid(candidates)

    def test_it(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        surface1 = bezier.Surface(nodes)
        surface2 = bezier.Surface(
            nodes, base_x=0.5, base_y=0.25, width=0.75)
        surface3 = bezier.Surface(
            nodes, base_x=0.25, base_y=1.25, width=0.5)

        centroid_x, centroid_y = self._call_function_under_test(
            [surface1, surface2, surface3])
        self.assertEqual(centroid_x, 0.5)
        self.assertEqual(centroid_y, 0.75)


class Test_jacobian_s(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, dimension):
        from bezier import _surface_helpers

        return _surface_helpers.jacobian_s(nodes, degree, dimension)

    def test_linear(self):
        nodes = np.array([[0.0, 1.0, np.nan]])
        nodes = nodes.T  # pylint: disable=no-member
        result = self._call_function_under_test(nodes, 1, 1)
        expected = np.array([[1.0]])
        self.assertEqual(result, expected)

    def test_quadratic(self):
        nodes = np.array([
            [0.0, 1.0],
            [1.0, 11.0],
            [5.0, 7.0],
            [4.0, -2.0],
            [-1.0, 6.0],
            [np.nan, np.nan],
        ])
        result = self._call_function_under_test(nodes, 2, 2)
        expected = 2.0 * np.array([
            [1.0, 10.0],
            [4.0, -4.0],
            [-5.0, 8.0],
        ])
        self.assertEqual(result, expected)

    def test_cubic(self):
        nodes = np.arange(10, dtype=float)[:, np.newaxis]**2
        result = self._call_function_under_test(nodes, 3, 1)
        expected = 3 * np.array([[1, 3, 5, 9, 11, 15]], dtype=float)
        expected = expected.T  # pylint: disable=no-member
        self.assertEqual(result, expected)

    def test_quartic(self):
        nodes = np.arange(15, dtype=float)[:, np.newaxis]**2
        result = self._call_function_under_test(nodes, 4, 1)
        expected = 4 * np.array([
            [1, 3, 5, 7, 11, 13, 15, 19, 21, 25]], dtype=float)
        expected = expected.T  # pylint: disable=no-member
        self.assertEqual(result, expected)


class Test_jacobian_t(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, dimension):
        from bezier import _surface_helpers

        return _surface_helpers.jacobian_t(nodes, degree, dimension)

    def test_linear(self):
        nodes = np.array([[0.0, np.nan, 1.0]])
        nodes = nodes.T  # pylint: disable=no-member
        result = self._call_function_under_test(nodes, 1, 1)
        expected = np.array([[1.0]])
        self.assertEqual(result, expected)

    def test_quadratic(self):
        nodes = np.array([
            [4.0, -2.0],
            [0.0, 1.0],
            [np.nan, np.nan],
            [5.0, 7.0],
            [-1.0, 6.0],
            [1.0, 12.0],
        ])
        result = self._call_function_under_test(nodes, 2, 2)
        expected = 2.0 * np.array([
            [1.0, 9.0],
            [-1.0, 5.0],
            [-4.0, 5.0],
        ])
        self.assertEqual(result, expected)

    def test_cubic(self):
        nodes = np.arange(10, dtype=float)[:, np.newaxis]**2
        result = self._call_function_under_test(nodes, 3, 1)
        expected = 3 * np.array([[16, 24, 32, 33, 39, 32]], dtype=float)
        expected = expected.T  # pylint: disable=no-member
        self.assertEqual(result, expected)

    def test_quartic(self):
        nodes = np.arange(15, dtype=float)[:, np.newaxis]**2
        result = self._call_function_under_test(nodes, 4, 1)
        expected = 4 * np.array([
            [25, 35, 45, 55, 56, 64, 72, 63, 69, 52]], dtype=float)
        expected = expected.T  # pylint: disable=no-member
        self.assertEqual(result, expected)


class Test_newton_refine(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(surface, x_val, y_val, s, t):
        from bezier import _surface_helpers

        return _surface_helpers.newton_refine(surface, x_val, y_val, s, t)

    def test_it(self):
        import bezier

        surface = bezier.Surface(np.array([
            [0.0, 0.0],
            [0.5, -0.25],
            [1.0, 0.0],
            [0.0, 0.5],
            [0.5, 0.5],
            [-0.25, 0.875],
        ]))
        # This surface is given by
        #     [(4 s - t^2) / 4, (4 s^2 + 4 s t - t^2 - 4 s + 8 t) / 8]
        s, t = 0.25, 0.5
        # At our points, the Jacobian is
        #     [1, -1/4]
        #     [0,  1  ]
        # hence there will be no round-off when applying the inverse.
        x_val, y_val = surface.evaluate_cartesian(0.5, 0.25)
        new_s, new_t = self._call_function_under_test(
            surface, x_val, y_val, s, t)
        self.assertEqual(new_s, 247.0 / 512.0)
        self.assertEqual(new_t, 31.0 / 128.0)


class Test_locate_point(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(surface, x_val, y_val):
        from bezier import _surface_helpers

        return _surface_helpers.locate_point(surface, x_val, y_val)

    def test_it(self):
        import bezier

        surface = bezier.Surface(UNIT_TRIANGLE)
        x_val = 0.25
        y_val = 0.625
        s, t = self._call_function_under_test(surface, x_val, y_val)
        self.assertEqual(s, x_val)
        self.assertEqual(t, y_val)

    def test_extra_newton_step(self):
        import bezier

        surface = bezier.Surface(UNIT_TRIANGLE)
        x_val = 1.0 / 3.0
        y_val = 1.0 / 3.0
        with mock.patch('bezier._surface_helpers._LOCATE_EPS', new=-1.0):
            s, t = self._call_function_under_test(surface, x_val, y_val)

        self.assertEqual(s, x_val)
        self.assertEqual(t, y_val)

    def test_no_match(self):
        import bezier

        surface = bezier.Surface(UNIT_TRIANGLE)
        x_val = -0.125
        y_val = 0.25
        self.assertIsNone(
            self._call_function_under_test(surface, x_val, y_val))


class Test_classify_intersection(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection):
        from bezier import _surface_helpers

        return _surface_helpers.classify_intersection(intersection)

    def test_simple(self):
        import bezier

        left = bezier.Curve(np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ]))
        right = bezier.Curve(np.array([
            [0.25, 0.0],
            [0.75, 1.0],
        ]))
        intersection = make_intersect(left, 0.5, right, 0.5)
        result = self._call_function_under_test(intersection)
        self.assertIs(result, get_enum('second'))

        # Swap and classify.
        intersection = make_intersect(right, 0.5, left, 0.5)
        result = self._call_function_under_test(intersection)
        self.assertIs(result, get_enum('first'))

    def test_corner_end(self):
        import bezier

        left = bezier.Curve(np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ]))
        right = bezier.Curve(np.array([
            [1.0, 0.0],
            [1.0, 2.0],
        ]))
        intersection = make_intersect(left, 1.0, right, 0.5)
        with self.assertRaises(ValueError):
            self._call_function_under_test(intersection)

    def test_corner_start(self):
        import bezier

        left = bezier.Curve(np.array([
            [1.0, 1.0],
            [0.0, 0.0],
        ]))
        right = bezier.Curve(np.array([
            [1.0, 0.0],
            [1.0, 2.0],
        ]))
        intersection = make_intersect(left, 0.0, right, 0.5)
        result = self._call_function_under_test(intersection)
        self.assertIs(result, get_enum('first'))

    def test_tangent(self):
        import bezier

        left = bezier.Curve(np.array([
            [0.0, 0.0],
            [1.5, 1.0],
            [3.0, 0.0],
        ]))
        right = bezier.Curve(np.array([
            [1.0, 0.0],
            [1.5, 1.0],
            [2.0, 0.0],
        ]))
        intersection = make_intersect(left, 0.5, right, 0.5)
        result = self._call_function_under_test(intersection)
        self.assertIs(result, get_enum('tangent_first'))


class Test__classify_tangent_intersection(unittest.TestCase):

    QUADRATIC1 = np.array([
        [1.0, 0.0],
        [1.5, 1.0],
        [2.0, 0.0],
    ])
    QUADRATIC2 = np.array([
        [0.0, 0.0],
        [1.5, 1.0],
        [3.0, 0.0],
    ])
    QUADRATIC3 = np.array([
        [1.0, 1.0],
        [1.5, 0.0],
        [2.0, 1.0],
    ])

    @staticmethod
    def _call_function_under_test(intersection, tangent1, tangent2):
        from bezier import _surface_helpers

        return _surface_helpers._classify_tangent_intersection(
            intersection, tangent1, tangent2)

    def _call_helper(self, intersection):
        from bezier import _curve_helpers

        tangent1 = _curve_helpers.evaluate_hodograph(
            intersection.left._nodes,
            intersection.left.degree, intersection.s)
        tangent2 = _curve_helpers.evaluate_hodograph(
            intersection.right._nodes,
            intersection.right.degree, intersection.t)

        return self._call_function_under_test(
            intersection, tangent1, tangent2)

    def test_left_curvature(self):
        import bezier

        left = bezier.Curve(self.QUADRATIC1[::-1, :])
        right = bezier.Curve(self.QUADRATIC2[::-1, :])
        intersection = make_intersect(left, 0.5, right, 0.5)

        result = self._call_helper(intersection)
        self.assertIs(result, get_enum('tangent_first'))

    def test_right_curvature(self):
        import bezier

        left = bezier.Curve(self.QUADRATIC1)
        right = bezier.Curve(self.QUADRATIC2)
        intersection = make_intersect(left, 0.5, right, 0.5)

        result = self._call_helper(intersection)
        self.assertIs(result, get_enum('tangent_second'))

    def test_same_direction_same_curvature(self):
        import bezier

        left = bezier.Curve(np.array([
            [1.0, 0.25],
            [-0.5, -0.25],
            [0.0, 0.25],
        ]))
        right = bezier.Curve(np.array([
            [0.75, 0.25],
            [-0.25, -0.25],
            [-0.25, 0.25],
        ]))
        intersection = make_intersect(left, 0.5, right, 0.5)
        with self.assertRaises(NotImplementedError):
            self._call_helper(intersection)

    def test_opposed_same_curvature(self):
        import bezier

        left = bezier.Curve(np.array([
            [0.0, 0.25],
            [-0.5, -0.25],
            [1.0, 0.25],
        ]))
        right = bezier.Curve(np.array([
            [0.75, 0.25],
            [-0.25, -0.25],
            [-0.25, 0.25],
        ]))
        intersection = make_intersect(left, 0.5, right, 0.5)
        with self.assertRaises(NotImplementedError):
            self._call_helper(intersection)

    def test_opposed_same_sign_curvature_no_overlap(self):
        import bezier

        left = bezier.Curve(self.QUADRATIC1[::-1, :])
        right = bezier.Curve(self.QUADRATIC3)
        intersection = make_intersect(left, 0.5, right, 0.5)

        result = self._call_helper(intersection)
        self.assertIs(result, get_enum('opposed'))

    def test_opposed_same_sign_curvature_with_overlap(self):
        import bezier

        left = bezier.Curve(self.QUADRATIC1)
        right = bezier.Curve(self.QUADRATIC3[::-1, :])
        intersection = make_intersect(left, 0.5, right, 0.5)

        with self.assertRaises(NotImplementedError):
            self._call_helper(intersection)

    def test_opposed_opp_sign_curvature_no_overlap(self):
        import bezier

        left = bezier.Curve(self.QUADRATIC1[::-1, :])
        right = bezier.Curve(self.QUADRATIC2)
        intersection = make_intersect(left, 0.5, right, 0.5)

        result = self._call_helper(intersection)
        self.assertIs(result, get_enum('opposed'))

    def test_opposed_opp_sign_curvature_with_overlap(self):
        import bezier

        left = bezier.Curve(self.QUADRATIC1)
        right = bezier.Curve(self.QUADRATIC2[::-1, :])
        intersection = make_intersect(left, 0.5, right, 0.5)

        with self.assertRaises(NotImplementedError):
            self._call_helper(intersection)


class Test_edge_cycle(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(edge1, edge2, edge3):
        from bezier import _surface_helpers

        return _surface_helpers.edge_cycle(edge1, edge2, edge3)

    def test_it(self):
        edge1 = mock.Mock()
        edge2 = mock.Mock()
        edge3 = mock.Mock()
        self.assertIsNone(
            self._call_function_under_test(edge1, edge2, edge3))

        self.assertEqual(edge1._edge_index, 0)
        self.assertIs(edge1._next_edge, edge2)
        self.assertIs(edge1._previous_edge, edge3)

        self.assertEqual(edge2._edge_index, 1)
        self.assertIs(edge2._next_edge, edge3)
        self.assertIs(edge2._previous_edge, edge1)

        self.assertEqual(edge3._edge_index, 2)
        self.assertIs(edge3._next_edge, edge1)
        self.assertIs(edge3._previous_edge, edge2)


class Test_handle_corners(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection):
        from bezier import _surface_helpers

        return _surface_helpers.handle_corners(intersection)

    def test_neither(self):
        left = mock.Mock()
        right = mock.Mock()
        intersection = make_intersect(left, 0.5, right, 0.5)

        self.assertFalse(self._call_function_under_test(intersection))
        self.assertEqual(intersection.s, 0.5)
        self.assertIs(intersection.left, left)
        self.assertEqual(intersection.t, 0.5)
        self.assertIs(intersection.right, right)

    def test_s(self):
        left = mock.Mock(next_edge=mock.sentinel.next_left)
        right = mock.Mock()
        intersection = make_intersect(left, 1.0, right, 0.25)

        self.assertTrue(self._call_function_under_test(intersection))
        self.assertEqual(intersection.s, 0.0)
        self.assertIs(intersection.left, mock.sentinel.next_left)
        self.assertEqual(intersection.t, 0.25)
        self.assertIs(intersection.right, right)

    def test_t(self):
        left = mock.Mock()
        right = mock.Mock(next_edge=mock.sentinel.next_right)
        intersection = make_intersect(left, 0.75, right, 1.0)

        self.assertTrue(self._call_function_under_test(intersection))
        self.assertEqual(intersection.s, 0.75)
        self.assertIs(intersection.left, left)
        self.assertEqual(intersection.t, 0.0)
        self.assertIs(intersection.right, mock.sentinel.next_right)


class Test_verify_duplicates(unittest.TestCase):

    LEFT = mock.Mock(edge_index=1)
    RIGHT = mock.Mock(edge_index=2)

    @staticmethod
    def _call_function_under_test(duplicates, uniques):
        from bezier import _surface_helpers

        return _surface_helpers.verify_duplicates(duplicates, uniques)

    def test_empty(self):
        self.assertIsNone(self._call_function_under_test([], []))

    def test_success(self):
        uniq = make_intersect(self.LEFT, 0.0, self.RIGHT, 0.25)
        self.assertIsNone(
            self._call_function_under_test([uniq], [uniq]))

    def test_success_triple(self):
        uniq = make_intersect(self.LEFT, 0.0, self.RIGHT, 0.0)
        self.assertIsNone(
            self._call_function_under_test([uniq, uniq, uniq], [uniq]))

    def test_failed_uniqueness(self):
        uniq = make_intersect(self.LEFT, 0.375, self.RIGHT, 0.75)
        with self.assertRaises(ValueError):
            self._call_function_under_test([], [uniq, uniq])

    def test_bad_duplicate(self):
        dupe = make_intersect(self.LEFT, 0.75, self.RIGHT, 0.25)
        uniq = make_intersect(self.LEFT, 0.25, self.RIGHT, 0.75)
        with self.assertRaises(ValueError):
            self._call_function_under_test([dupe], [uniq])

    def test_bad_single_corner(self):
        uniq = make_intersect(self.LEFT, 0.125, self.RIGHT, 0.125)
        with self.assertRaises(ValueError):
            self._call_function_under_test([uniq], [uniq])

    def test_bad_double_corner(self):
        uniq = make_intersect(self.LEFT, 0.0, self.RIGHT, 1.0)
        with self.assertRaises(ValueError):
            self._call_function_under_test([uniq, uniq, uniq], [uniq])

    def test_bad_count(self):
        uniq = make_intersect(self.LEFT, 0.375, self.RIGHT, 0.75)
        with self.assertRaises(ValueError):
            self._call_function_under_test([uniq, uniq], [uniq])


class Test_combine_intersections(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersections, surface1, surface2):
        from bezier import _surface_helpers

        return _surface_helpers.combine_intersections(
            intersections, surface1, surface2)

    def test_empty(self):
        import bezier

        surface1 = bezier.Surface(UNIT_TRIANGLE)
        surface2 = bezier.Surface(np.array([
            [-1.0, 0.0],
            [-1.0, 1.0],
            [-2.0, 1.0],
        ]))
        result = self._call_function_under_test([], surface1, surface2)
        self.assertEqual(result, [])
