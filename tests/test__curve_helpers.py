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


class Test_make_subdivision_matrix(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(degree):
        from bezier import _curve_helpers

        return _curve_helpers.make_subdivision_matrix(degree)

    def _helper(self, degree, expected):
        result = self._call_function_under_test(degree)
        self.assertTrue(np.all(result == expected))

    def test_linear(self):
        from bezier import curve

        self._helper(1, curve._LINEAR_SUBDIVIDE)

    def test_quadratic(self):
        from bezier import curve

        self._helper(2, curve._QUADRATIC_SUBDIVIDE)

    def test_cubic(self):
        from bezier import curve

        self._helper(3, curve._CUBIC_SUBDIVIDE)


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


class Test__vec_size(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, s_val):
        from bezier import _curve_helpers

        return _curve_helpers._vec_size(nodes, degree, s_val)

    def test_linear(self):
        nodes = np.array([
            [0.0, 0.0],
            [3.0, -4.0],
        ])
        size = self._call_function_under_test(nodes, 1, 0.25)
        self.assertEqual(size, 0.25 * 5.0)

    def test_quadratic(self):
        nodes = np.array([
            [0.0, 0.0],
            [2.0, 3.0],
            [1.0, 6.0],
        ])
        size = self._call_function_under_test(nodes, 2, 0.5)
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


class Test_elevate_nodes(unittest.TestCase):

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
        self.assertTrue(np.all(result == expected))

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
        self.assertTrue(np.all(result == expected))


class Test_de_casteljau_one_round(unittest.TestCase):

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
        self.assertTrue(np.all(result == np.array([[2.25, 4.0]])))


class Test_specialize_curve(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, start, end):
        from bezier import _curve_helpers

        return _curve_helpers.specialize_curve(
            nodes, degree, start, end)

    def test_it(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        result = self._call_function_under_test(nodes, 1, 0.25, 0.75)
        expected = np.array([
            [0.25, 0.25],
            [0.75, 0.75],
        ])
        self.assertTrue(np.all(result == expected))

    def test_againt_subdivision(self):
        import bezier

        nodes = np.array([
            [0.0, 1.0],
            [1.0, 6.0],
            [3.0, 5.0],
        ])
        curve = bezier.Curve(nodes)
        left, right = curve.subdivide()

        left_nodes = self._call_function_under_test(nodes, 2, 0.0, 0.5)
        self.assertTrue(np.all(left.nodes == left_nodes))

        right_nodes = self._call_function_under_test(nodes, 2, 0.5, 1.0)
        self.assertTrue(np.all(right.nodes == right_nodes))
