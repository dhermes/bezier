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


class Test_vector_close(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(vec1, vec2, **kwargs):
        from bezier import _helpers

        return _helpers.vector_close(vec1, vec2, **kwargs)

    def test_identical(self):
        vec1 = np.array([0.5, 4.0])
        self.assertTrue(self._call_function_under_test(vec1, vec1))

    def test_far_apart(self):
        vec1 = np.array([0.0, 6.0])
        vec2 = np.array([1.0, -4.0])
        self.assertFalse(self._call_function_under_test(vec1, vec2))

    def test_close_but_different(self):
        vec1 = np.array([2.25, -3.5])
        vec2 = vec1 + np.array([-5.0, 12.0]) / 2.0**43
        self.assertTrue(self._call_function_under_test(vec1, vec2))

    def test_custom_epsilon(self):
        vec1 = np.array([3.0, 4.0])
        vec2 = np.array([2.0, 5.0])
        self.assertTrue(self._call_function_under_test(vec1, vec2, eps=0.5))
        self.assertFalse(self._call_function_under_test(vec1, vec2))

    def test_near_zero(self):
        vec1 = np.array([0.0, 0.0])
        vec2 = np.array([3.0, 4.0]) / 2.0**45
        self.assertTrue(self._call_function_under_test(vec1, vec2))

    def test_near_zero_fail(self):
        vec1 = np.array([1.0, 0.0]) / 2.0**20
        vec2 = np.array([0.0, 0.0])
        self.assertFalse(self._call_function_under_test(vec1, vec2))


class Test_in_interval(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(value, start, end):
        from bezier import _helpers

        return _helpers.in_interval(value, start, end)

    def test_interior(self):
        self.assertTrue(self._call_function_under_test(
            1.5, 1.0, 2.0))

    def test_barely_inside(self):
        local_epsilon = np.spacing(1.0)  # pylint: disable=no-member
        self.assertTrue(self._call_function_under_test(
            1.0 + local_epsilon, 1.0, 2.0))

    def test_barely_outside(self):
        local_epsilon = np.spacing(1.0)  # pylint: disable=no-member
        self.assertFalse(self._call_function_under_test(
            1.0 - local_epsilon / 2.0, 1.0, 2.0))

    def test_outside(self):
        self.assertFalse(self._call_function_under_test(
            -1.0, 1.0, 2.0))


class Test_bbox(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _helpers

        return _helpers.bbox(nodes)

    def test_it(self):
        nodes = np.array([
            [0.0, 5.0],
            [1.0, 3.0],
        ])
        left, right, bottom, top = self._call_function_under_test(nodes)
        self.assertEqual(left, 0.0)
        self.assertEqual(right, 1.0)
        self.assertEqual(bottom, 3.0)
        self.assertEqual(top, 5.0)


class Test_contains(unittest.TestCase):

    UNIT_SQUARE = np.array([
        [0.0, 0.0],
        [1.0, 1.0],
    ])

    @staticmethod
    def _call_function_under_test(nodes, x_val, y_val):
        from bezier import _helpers

        return _helpers.contains(nodes, x_val, y_val)

    def test_x_outside(self):
        self.assertFalse(
            self._call_function_under_test(self.UNIT_SQUARE, -1.0, 0.5))

    def test_y_outside(self):
        self.assertFalse(
            self._call_function_under_test(self.UNIT_SQUARE, 0.625, 4.25))

    def test_inside(self):
        self.assertTrue(
            self._call_function_under_test(self.UNIT_SQUARE, 0.25, 0.75))


class Test_contains_nd(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, point):
        from bezier import _helpers

        return _helpers.contains_nd(nodes, point)

    def test_below(self):
        nodes = np.array([
            [0.0, 1.0],
            [0.5, 0.0],
            [1.0, 2.0],
        ])
        point = np.array([[-0.5, 1.0]])
        self.assertFalse(self._call_function_under_test(nodes, point))

    def test_above(self):
        nodes = np.array([
            [0.0, -4.0, 2.0],
            [-1.0, 1.0, 3.0],
        ])
        point = np.array([[-0.5, 2.0, 2.5]])
        self.assertFalse(self._call_function_under_test(nodes, point))

    def test_inside(self):
        nodes = np.array([
            [0.0, 1.0, 2.0, 3.0],
            [1.0, -2.0, -4.0, 1.0],
        ])
        point = np.array([[0.5, 0.0, 0.0, 2.0]])
        self.assertTrue(self._call_function_under_test(nodes, point))


class Test_cross_product(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(vec0, vec1):
        from bezier import _helpers

        return _helpers.cross_product(vec0, vec1)

    def test_it(self):
        vec0 = np.array([[1.0, 7.0]]) / 8.0
        vec1 = np.array([[-11.0, 24.0]]) / 32.0
        result = self._call_function_under_test(vec0, vec1)

        vec0_as_3d = np.hstack([vec0, [[0.0]]])
        vec1_as_3d = np.hstack([vec1, [[0.0]]])

        actual_cross = np.cross(vec0_as_3d, vec1_as_3d)
        expected = np.array([[0.0, 0.0, result]])
        self.assertEqual(actual_cross, expected)


class Test_n_bits_away(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(value1, value2, **kwargs):
        from bezier import _helpers

        return _helpers.n_bits_away(value1, value2, **kwargs)

    def test_first_zero(self):
        self.assertTrue(
            self._call_function_under_test(0.0, 0.0))

    def test_second_zero(self):
        self.assertFalse(
            self._call_function_under_test(1.0, 0.0))

    def test_default_single_bit(self):
        value1 = 1.0
        value2 = 1.0 + 2.0**(-52)
        value3 = 1.0 - 2.0**(-53)

        self.assertTrue(
            self._call_function_under_test(value1, value1))
        self.assertTrue(
            self._call_function_under_test(value1, value2))
        self.assertTrue(
            self._call_function_under_test(value2, value1))
        self.assertTrue(
            self._call_function_under_test(value1, value3))
        self.assertTrue(
            self._call_function_under_test(value3, value1))

    def test_non_default_n(self):
        value1 = 1.5
        value2 = 1.5 + 2.0**(-43)

        self.assertFalse(
            self._call_function_under_test(value1, value2))
        self.assertTrue(
            self._call_function_under_test(value1, value2, num_bits=1000))

    def test_very_close(self):
        value1 = float.fromhex('0x1.fffffffffffd3p-3')
        value2 = float.fromhex('0x1.fffffffffffcep-3')

        self.assertFalse(
            self._call_function_under_test(value1, value2))
        self.assertFalse(
            self._call_function_under_test(value1, value2, num_bits=4))
        self.assertTrue(
            self._call_function_under_test(value1, value2, num_bits=5))
