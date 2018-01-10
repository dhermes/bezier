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

import numpy as np

from tests.unit import utils


SPACING = np.spacing  # pylint: disable=no-member


class Test__vector_close(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(vec1, vec2, **kwargs):
        from bezier import _helpers

        return _helpers._vector_close(vec1, vec2, **kwargs)

    def test_identical(self):
        vec1 = np.asfortranarray([[0.5, 4.0]])
        self.assertTrue(self._call_function_under_test(vec1, vec1))

    def test_far_apart(self):
        vec1 = np.asfortranarray([[0.0, 6.0]])
        vec2 = np.asfortranarray([[1.0, -4.0]])
        self.assertFalse(self._call_function_under_test(vec1, vec2))

    def test_close_but_different(self):
        vec1 = np.asfortranarray([[2.25, -3.5]])
        vec2 = vec1 + np.asfortranarray([[-5.0, 12.0]]) / 2.0**43
        self.assertTrue(self._call_function_under_test(vec1, vec2))

    def test_custom_epsilon(self):
        vec1 = np.asfortranarray([[3.0, 4.0]])
        vec2 = np.asfortranarray([[2.0, 5.0]])
        self.assertTrue(self._call_function_under_test(vec1, vec2, eps=0.5))
        self.assertFalse(self._call_function_under_test(vec1, vec2))

    def test_near_zero(self):
        vec1 = np.asfortranarray([[0.0, 0.0]])
        vec2 = np.asfortranarray([[3.0, 4.0]]) / 2.0**45
        self.assertTrue(self._call_function_under_test(vec1, vec2))

    def test_near_zero_fail(self):
        vec1 = np.asfortranarray([[1.0, 0.0]]) / 2.0**20
        vec2 = np.asfortranarray([[0.0, 0.0]])
        self.assertFalse(self._call_function_under_test(vec1, vec2))


@utils.needs_helpers_speedup
class Test_speedup_vector_close(Test__vector_close):

    @staticmethod
    def _call_function_under_test(vec1, vec2, **kwargs):
        from bezier import _helpers_speedup

        return _helpers_speedup.vector_close(vec1, vec2, **kwargs)


class Test__in_interval(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(value, start, end):
        from bezier import _helpers

        return _helpers._in_interval(value, start, end)

    def test_interior(self):
        self.assertTrue(self._call_function_under_test(
            1.5, 1.0, 2.0))

    def test_barely_inside(self):
        local_epsilon = SPACING(1.0)
        self.assertTrue(self._call_function_under_test(
            1.0 + local_epsilon, 1.0, 2.0))

    def test_barely_outside(self):
        local_epsilon = SPACING(1.0)
        self.assertFalse(self._call_function_under_test(
            1.0 - local_epsilon / 2.0, 1.0, 2.0))

    def test_outside(self):
        self.assertFalse(self._call_function_under_test(
            -1.0, 1.0, 2.0))


@utils.needs_helpers_speedup
class Test_speedup_in_interval(Test__in_interval):

    @staticmethod
    def _call_function_under_test(value, start, end):
        from bezier import _helpers_speedup

        return _helpers_speedup.in_interval(value, start, end)


class Test__bbox(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _helpers

        return _helpers._bbox(nodes)

    def test_it(self):
        nodes = np.asfortranarray([
            [0.0, 5.0],
            [1.0, 3.0],
        ])
        left, right, bottom, top = self._call_function_under_test(nodes)
        self.assertEqual(left, 0.0)
        self.assertEqual(right, 1.0)
        self.assertEqual(bottom, 3.0)
        self.assertEqual(top, 5.0)

    def test_lots_of_values(self):
        nodes = np.asfortranarray([
            [1.0, 0.0],
            [2.0, 1.0],
            [-1.0, 2.0],
            [5.0, -3.0],
            [4.0, 4.0],
            [0.0, 0.0],
        ])
        left, right, bottom, top = self._call_function_under_test(nodes)
        self.assertEqual(left, -1.0)
        self.assertEqual(right, 5.0)
        self.assertEqual(bottom, -3.0)
        self.assertEqual(top, 4.0)


@utils.needs_helpers_speedup
class Test_speedup_bbox(Test__bbox):

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _helpers_speedup

        return _helpers_speedup.bbox(nodes)


class Test__contains_nd(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, point):
        from bezier import _helpers

        return _helpers._contains_nd(nodes, point)

    def test_below(self):
        nodes = np.asfortranarray([
            [0.0, 1.0],
            [0.5, 0.0],
            [1.0, 2.0],
        ])
        point = np.asfortranarray([[-0.5, 1.0]])
        self.assertFalse(self._call_function_under_test(nodes, point))

    def test_above(self):
        nodes = np.asfortranarray([
            [0.0, -4.0, 2.0],
            [-1.0, 1.0, 3.0],
        ])
        point = np.asfortranarray([[-0.5, 2.0, 2.5]])
        self.assertFalse(self._call_function_under_test(nodes, point))

    def test_inside(self):
        nodes = np.asfortranarray([
            [0.0, 1.0, 2.0, 3.0],
            [1.0, -2.0, -4.0, 1.0],
        ])
        point = np.asfortranarray([[0.5, 0.0, 0.0, 2.0]])
        self.assertTrue(self._call_function_under_test(nodes, point))

    def test_shape_mismatch(self):
        nodes = np.asfortranarray([
            [0.0, 1.0],
            [1.0, 3.0],
            [2.0, 6.0],
        ])
        point = np.asfortranarray([[0.0, 1.5, 1.0]])

        with self.assertRaises(ValueError):
            self._call_function_under_test(nodes, point)


@utils.needs_helpers_speedup
class Test_speedup_contains_nd(Test__contains_nd):

    @staticmethod
    def _call_function_under_test(nodes, point):
        from bezier import _helpers_speedup

        return _helpers_speedup.contains_nd(nodes, point)


class Test__cross_product(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(vec0, vec1):
        from bezier import _helpers

        return _helpers._cross_product(vec0, vec1)

    def test_it(self):
        vec0 = np.asfortranarray([[1.0, 7.0]]) / 8.0
        vec1 = np.asfortranarray([[-11.0, 24.0]]) / 32.0
        result = self._call_function_under_test(vec0, vec1)

        vec0_as_3d = np.zeros((1, 3), order='F')
        vec0_as_3d[0, :2] = vec0
        vec1_as_3d = np.zeros((1, 3), order='F')
        vec1_as_3d[0, :2] = vec1

        actual_cross = np.asfortranarray(np.cross(vec0_as_3d, vec1_as_3d))
        expected = np.asfortranarray([[0.0, 0.0, result]])
        self.assertEqual(actual_cross, expected)


@utils.needs_helpers_speedup
class Test_speedup_cross_product(Test__cross_product):

    @staticmethod
    def _call_function_under_test(vec0, vec1):
        from bezier import _helpers_speedup

        return _helpers_speedup.cross_product(vec0, vec1)


class Test__ulps_away(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(value1, value2, **kwargs):
        from bezier import _helpers

        return _helpers._ulps_away(value1, value2, **kwargs)

    def test_first_zero(self):
        self.assertTrue(
            self._call_function_under_test(0.0, 0.0))

    def test_second_zero(self):
        self.assertFalse(
            self._call_function_under_test(1.0, 0.0))

    def test_first_both_signs(self):
        delta = 0.5**52
        self.assertTrue(
            self._call_function_under_test(1.0, 1.0 + delta))
        self.assertTrue(
            self._call_function_under_test(-1.0, -1.0 - delta))

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


@utils.needs_helpers_speedup
class Test_speedup_ulps_away(Test__ulps_away):

    @staticmethod
    def _call_function_under_test(value1, value2, **kwargs):
        from bezier import _helpers_speedup

        return _helpers_speedup.ulps_away(value1, value2, **kwargs)


class Test_matrix_product(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(mat1, mat2):
        from bezier import _helpers

        return _helpers.matrix_product(mat1, mat2)

    def test_it(self):
        mat1 = np.asfortranarray([
            [1.0, 2.0],
            [3.0, 4.0],
            [5.0, 6.0],
        ])
        mat2 = np.asfortranarray([
            [7.0, 8.0, 9.0],
            [10.0, 11.0, 12.0],
        ])
        result = self._call_function_under_test(mat1, mat2)

        expected = np.asfortranarray([
            [27.0, 30.0, 33.0],
            [61.0, 68.0, 75.0],
            [95.0, 106.0, 117.0],
        ])
        self.assertEqual(result, expected)
        # Make sure our data is F-contiguous.
        self.assertTrue(result.flags.f_contiguous)
        self.assertFalse(result.flags.c_contiguous)
        # matrix_product() has the side-effect of returning a "view"
        # since it returns the transpose of a product of transposes.
        self.assertFalse(result.flags.owndata)


class Test__wiggle_interval(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(value, **kwargs):
        from bezier import _helpers

        return _helpers._wiggle_interval(value, **kwargs)

    def test_at_endpoint(self):
        # Really just making sure the function doesn't raise.
        result, success = self._call_function_under_test(0.0)
        self.assertTrue(success)
        self.assertEqual(result, 0.0)

        result, success = self._call_function_under_test(1.0)
        self.assertTrue(success)
        self.assertEqual(result, 1.0)

    def test_near_endpoint(self):
        _, success = self._call_function_under_test(1.0 + 2.0**(-20))
        self.assertFalse(success)

    def test_outside_below(self):
        _, success = self._call_function_under_test(-0.25)
        self.assertFalse(success)

    def test_outside_above(self):
        _, success = self._call_function_under_test(1.5)
        self.assertFalse(success)

    def test_valid(self):
        # Really just making sure the function doesn't raise.
        result, success = self._call_function_under_test(0.25)
        self.assertTrue(success)
        self.assertEqual(result, 0.25)

    def test_wiggle_below(self):
        value = -2.0**(-60)
        result, success = self._call_function_under_test(value)
        self.assertTrue(success)
        self.assertEqual(result, 0.0)

    def test_wiggle_above(self):
        value = 1 + 2.0**(-52)
        result, success = self._call_function_under_test(value)
        self.assertTrue(success)
        self.assertEqual(result, 1.0)

    def test_outer_boundary(self):
        # Values near / at the left-hand boundary.
        value = float.fromhex('-0x1.ffffffffffffep-46')
        self.assertEqual(
            self._call_function_under_test(value), (0.0, True))
        value = float.fromhex('-0x1.fffffffffffffp-46')
        self.assertEqual(
            self._call_function_under_test(value), (0.0, True))

        value = float.fromhex('-0x1.0000000000000p-45')
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)

        value = float.fromhex('-0x1.0000000000001p-45')
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)

        # Values near / at the right-hand boundary.
        value = float.fromhex('0x1.000000000007ep+0')
        self.assertEqual(
            self._call_function_under_test(value), (1.0, True))
        value = float.fromhex('0x1.000000000007fp+0')
        self.assertEqual(
            self._call_function_under_test(value), (1.0, True))

        value = float.fromhex('0x1.0000000000080p+0')
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)

        value = float.fromhex('0x1.0000000000081p+0')
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)

    def test_inner_boundary(self):
        # Values near / at the left-hand boundary.
        value = float.fromhex('0x1.ffffffffffffep-46')
        self.assertEqual(
            self._call_function_under_test(value), (0.0, True))
        value = float.fromhex('0x1.fffffffffffffp-46')
        self.assertEqual(
            self._call_function_under_test(value), (0.0, True))
        value = float.fromhex('0x1.0000000000000p-45')
        self.assertEqual(
            self._call_function_under_test(value), (value, True))
        value = float.fromhex('0x1.0000000000001p-45')
        self.assertEqual(
            self._call_function_under_test(value), (value, True))

        # Values near / at the right-hand boundary.
        value = float.fromhex('0x1.ffffffffffefep-1')
        self.assertEqual(
            self._call_function_under_test(value), (value, True))
        value = float.fromhex('0x1.ffffffffffeffp-1')
        self.assertEqual(
            self._call_function_under_test(value), (value, True))
        value = float.fromhex('0x1.fffffffffff00p-1')
        self.assertEqual(
            self._call_function_under_test(value), (value, True))
        value = float.fromhex('0x1.fffffffffff01p-1')
        self.assertEqual(
            self._call_function_under_test(value), (1.0, True))
        value = float.fromhex('0x1.fffffffffff02p-1')
        self.assertEqual(
            self._call_function_under_test(value), (1.0, True))

    def test_custom_wiggle(self):
        value = 1.25
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)

        result, success = self._call_function_under_test(value, wiggle=0.5)
        self.assertTrue(success)
        self.assertEqual(result, 1.0)

        value = 0.875
        self.assertEqual(
            self._call_function_under_test(value), (value, True))
        self.assertEqual(
            self._call_function_under_test(value, wiggle=0.25), (1.0, True))


@utils.needs_helpers_speedup
class Test_speedup_wiggle_interval(Test__wiggle_interval):

    def _call_function_under_test(self, value, **kwargs):
        from bezier import _helpers_speedup

        self.assertEqual(kwargs, {})
        return _helpers_speedup.wiggle_interval(value, **kwargs)

    def test_custom_wiggle(self):
        # Fortran implementation doesn't support optional wiggle. This
        # isn't because Fortran **can't** (just use "optional"), it's just
        # to allow the compiler to pre-compute 1 + wiggle / 1 - wiggle
        # rather than having to deal with it at run-time.
        pass
