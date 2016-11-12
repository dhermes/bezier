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

from tests import utils


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
        import mock

        # pylint: disable=no-member
        bernstein = np.array([[1.0, 2.0, 3.0]]).T
        # pylint: enable=no-member
        with mock.patch('bezier._surface_helpers.MAX_SUBDIVISIONS', new=1):
            self._helper(bernstein, 1)

    def test_no_conclusion(self):
        import mock

        # pylint: disable=no-member
        bernstein = np.array([[-1.0, 1.0, 2.0]]).T
        # pylint: enable=no-member
        with mock.patch('bezier._surface_helpers.MAX_SUBDIVISIONS', new=0):
            with self.assertRaises(ValueError):
                self._helper(bernstein, None)


class Test_quadratic_jacobian_polynomial(unittest.TestCase):

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
        # pylint: disable=no-member
        expected = np.array([[0.0, 2.0, 0.0, -2.0, 2.0, 0.0]]).T
        # pylint: enable=no-member
        self.assertTrue(np.all(bernstein == expected))


class Test_de_casteljau_one_round(unittest.TestCase):

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
        self.assertTrue(np.all(result == expected))

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
        self.assertTrue(np.all(result == expected))

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
        self.assertTrue(np.all(result == expected))
