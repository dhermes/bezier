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
import six

from tests.unit import utils


SPACING = np.spacing  # pylint: disable=no-member


class Test_vector_close(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(vec1, vec2, **kwargs):
        from bezier import _py_helpers

        return _py_helpers.vector_close(vec1, vec2, **kwargs)

    def test_identical(self):
        vec1 = np.asfortranarray([0.5, 4.0])
        self.assertTrue(self._call_function_under_test(vec1, vec1))

    def test_far_apart(self):
        vec1 = np.asfortranarray([0.0, 6.0])
        vec2 = np.asfortranarray([1.0, -4.0])
        self.assertFalse(self._call_function_under_test(vec1, vec2))

    def test_close_but_different(self):
        vec1 = np.asfortranarray([2.25, -3.5])
        vec2 = vec1 + np.asfortranarray([-5.0, 12.0]) / 2.0 ** 43
        self.assertTrue(self._call_function_under_test(vec1, vec2))

    def test_custom_epsilon(self):
        vec1 = np.asfortranarray([3.0, 4.0])
        vec2 = np.asfortranarray([2.0, 5.0])
        self.assertTrue(self._call_function_under_test(vec1, vec2, eps=0.5))
        self.assertFalse(self._call_function_under_test(vec1, vec2))

    def test_near_zero(self):
        vec1 = np.asfortranarray([0.0, 0.0])
        vec2 = np.asfortranarray([3.0, 4.0]) / 2.0 ** 45
        self.assertTrue(self._call_function_under_test(vec1, vec2))

    def test_near_zero_fail(self):
        vec1 = np.asfortranarray([1.0, 0.0]) / 2.0 ** 20
        vec2 = np.asfortranarray([0.0, 0.0])
        self.assertFalse(self._call_function_under_test(vec1, vec2))


class Test_in_interval(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(value, start, end):
        from bezier import _py_helpers

        return _py_helpers.in_interval(value, start, end)

    def test_interior(self):
        self.assertTrue(self._call_function_under_test(1.5, 1.0, 2.0))

    def test_barely_inside(self):
        # pylint: disable=assignment-from-no-return
        local_epsilon = SPACING(1.0)
        # pylint: enable=assignment-from-no-return
        self.assertTrue(
            self._call_function_under_test(1.0 + local_epsilon, 1.0, 2.0)
        )

    def test_barely_outside(self):
        # pylint: disable=assignment-from-no-return
        local_epsilon = SPACING(1.0)
        # pylint: enable=assignment-from-no-return
        self.assertFalse(
            self._call_function_under_test(1.0 - local_epsilon / 2.0, 1.0, 2.0)
        )

    def test_outside(self):
        self.assertFalse(self._call_function_under_test(-1.0, 1.0, 2.0))


class Test_bbox(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _py_helpers

        return _py_helpers.bbox(nodes)

    def test_it(self):
        nodes = np.asfortranarray([[0.0, 1.0], [5.0, 3.0]])
        left, right, bottom, top = self._call_function_under_test(nodes)
        self.assertEqual(left, 0.0)
        self.assertEqual(right, 1.0)
        self.assertEqual(bottom, 3.0)
        self.assertEqual(top, 5.0)

    def test_lots_of_values(self):
        nodes = np.asfortranarray(
            [[1.0, 2.0, -1.0, 5.0, 4.0, 0.0], [0.0, 1.0, 2.0, -3.0, 4.0, 0.0]]
        )
        left, right, bottom, top = self._call_function_under_test(nodes)
        self.assertEqual(left, -1.0)
        self.assertEqual(right, 5.0)
        self.assertEqual(bottom, -3.0)
        self.assertEqual(top, 4.0)


class Test_contains_nd(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes, point):
        from bezier import _py_helpers

        return _py_helpers.contains_nd(nodes, point)

    def test_below(self):
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [1.0, 0.0, 2.0]])
        point = np.asfortranarray([-0.5, 1.0])
        self.assertFalse(self._call_function_under_test(nodes, point))

    def test_above(self):
        nodes = np.asfortranarray([[0.0, -1.0], [-4.0, 1.0], [2.0, 3.0]])
        point = np.asfortranarray([-0.5, 2.0, 2.5])
        self.assertFalse(self._call_function_under_test(nodes, point))

    def test_inside(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0], [1.0, -2.0], [2.0, -4.0], [3.0, 1.0]]
        )
        point = np.asfortranarray([0.5, 0.0, 0.0, 2.0])
        self.assertTrue(self._call_function_under_test(nodes, point))

    def test_shape_mismatch(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [1.0, 3.0, 6.0]])
        point = np.asfortranarray([0.0, 1.5, 1.0])
        with self.assertRaises(ValueError):
            self._call_function_under_test(nodes, point)


class Test_cross_product(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(vec0, vec1):
        from bezier import _py_helpers

        return _py_helpers.cross_product(vec0, vec1)

    def test_it(self):
        vec0 = np.asfortranarray([1.0, 7.0]) / 8.0
        vec1 = np.asfortranarray([-11.0, 24.0]) / 32.0
        result = self._call_function_under_test(vec0, vec1)
        vec0_as_3d = np.zeros((3,), order="F")
        vec0_as_3d[:2] = vec0
        vec1_as_3d = np.zeros((3,), order="F")
        vec1_as_3d[:2] = vec1
        actual_cross = np.cross(vec0_as_3d, vec1_as_3d)
        expected = np.asfortranarray([0.0, 0.0, result])
        self.assertEqual(actual_cross, expected)


class Test_matrix_product(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(mat1, mat2):
        from bezier import _py_helpers

        return _py_helpers.matrix_product(mat1, mat2)

    def test_it(self):
        mat1 = np.asfortranarray([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
        mat2 = np.asfortranarray([[7.0, 8.0, 9.0], [10.0, 11.0, 12.0]])
        result = self._call_function_under_test(mat1, mat2)
        expected = np.asfortranarray(
            [[27.0, 30.0, 33.0], [61.0, 68.0, 75.0], [95.0, 106.0, 117.0]]
        )
        self.assertEqual(result, expected)
        # Make sure our data is F-contiguous.
        self.assertTrue(result.flags.f_contiguous)
        self.assertFalse(result.flags.c_contiguous)
        # matrix_product() has the side-effect of returning a "view"
        # since it returns the transpose of a product of transposes.
        self.assertFalse(result.flags.owndata)


class Test_wiggle_interval(unittest.TestCase):
    WIGGLE = 0.5 ** 44
    MACHINE_EPS = 0.5 ** 52

    @staticmethod
    def _call_function_under_test(value, **kwargs):
        from bezier import _py_helpers

        return _py_helpers.wiggle_interval(value, **kwargs)

    def test_at_endpoint(self):
        # Really just making sure the function doesn't raise.
        result, success = self._call_function_under_test(0.0)
        self.assertTrue(success)
        self.assertEqual(result, 0.0)
        result, success = self._call_function_under_test(1.0)
        self.assertTrue(success)
        self.assertEqual(result, 1.0)

    def test_near_endpoint(self):
        _, success = self._call_function_under_test(1.0 + 0.5 ** 20)
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
        value = -(0.5 ** 60)
        result, success = self._call_function_under_test(value)
        self.assertTrue(success)
        self.assertEqual(result, 0.0)

    def test_wiggle_above(self):
        value = 1.0 + self.MACHINE_EPS
        result, success = self._call_function_under_test(value)
        self.assertTrue(success)
        self.assertEqual(result, 1.0)

    def test_outer_boundary(self):
        # Values near / at the left-hand boundary.
        value = -self.WIGGLE + self.MACHINE_EPS * self.WIGGLE
        self.assertEqual(self._call_function_under_test(value), (0.0, True))
        value = -self.WIGGLE + self.MACHINE_EPS * self.WIGGLE / 2
        self.assertEqual(self._call_function_under_test(value), (0.0, True))
        value = -self.WIGGLE
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)
        value = -self.WIGGLE - self.MACHINE_EPS * self.WIGGLE
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)
        # Values near / at the right-hand boundary.
        value = 1.0 + self.WIGGLE - 2 * self.MACHINE_EPS
        self.assertEqual(self._call_function_under_test(value), (1.0, True))
        value = 1.0 + self.WIGGLE - self.MACHINE_EPS
        self.assertEqual(self._call_function_under_test(value), (1.0, True))
        value = 1.0 + self.WIGGLE
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)
        value = 1.0 + self.WIGGLE + self.MACHINE_EPS
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)

    def test_inner_boundary(self):
        # Values near / at the left-hand boundary.
        value = self.WIGGLE - self.WIGGLE * self.MACHINE_EPS
        self.assertEqual(self._call_function_under_test(value), (0.0, True))
        value = self.WIGGLE - self.WIGGLE * self.MACHINE_EPS / 2
        self.assertEqual(self._call_function_under_test(value), (0.0, True))
        value = self.WIGGLE
        self.assertEqual(self._call_function_under_test(value), (value, True))
        value = self.WIGGLE + self.WIGGLE * self.MACHINE_EPS
        self.assertEqual(self._call_function_under_test(value), (value, True))
        # Values near / at the right-hand boundary.
        value = 1.0 - self.WIGGLE - self.MACHINE_EPS
        self.assertEqual(self._call_function_under_test(value), (value, True))
        value = 1.0 - self.WIGGLE - self.MACHINE_EPS / 2
        self.assertEqual(self._call_function_under_test(value), (value, True))
        value = 1.0 - self.WIGGLE
        self.assertEqual(self._call_function_under_test(value), (value, True))
        value = 1.0 - self.WIGGLE + self.MACHINE_EPS / 2
        self.assertEqual(self._call_function_under_test(value), (1.0, True))
        value = 1.0 - self.WIGGLE + self.MACHINE_EPS
        self.assertEqual(self._call_function_under_test(value), (1.0, True))

    def test_custom_wiggle(self):
        value = 1.25
        _, success = self._call_function_under_test(value)
        self.assertFalse(success)
        result, success = self._call_function_under_test(value, wiggle=0.5)
        self.assertTrue(success)
        self.assertEqual(result, 1.0)
        value = 0.875
        self.assertEqual(self._call_function_under_test(value), (value, True))
        self.assertEqual(
            self._call_function_under_test(value, wiggle=0.25), (1.0, True)
        )


class Test_cross_product_compare(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(start, candidate1, candidate2):
        from bezier import _py_helpers

        return _py_helpers.cross_product_compare(start, candidate1, candidate2)

    def test_it(self):
        start = np.asfortranarray([0.0, 0.0])
        candidate1 = np.asfortranarray([1.0, 0.0])
        candidate2 = np.asfortranarray([1.0, 1.0])
        result = self._call_function_under_test(start, candidate1, candidate2)
        self.assertEqual(result, 1.0)


class Test_in_sorted(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(values, value):
        from bezier import _py_helpers

        return _py_helpers.in_sorted(values, value)

    def test_inside(self):
        values = [0, 5, 8, 12, 17]
        self.assertFalse(self._call_function_under_test(values, 4))
        self.assertTrue(self._call_function_under_test(values, 5))
        self.assertFalse(self._call_function_under_test(values, 6))
        self.assertFalse(self._call_function_under_test(values, 7))
        self.assertTrue(self._call_function_under_test(values, 8))
        self.assertFalse(self._call_function_under_test(values, 9))
        self.assertFalse(self._call_function_under_test(values, 10))
        self.assertFalse(self._call_function_under_test(values, 11))
        self.assertTrue(self._call_function_under_test(values, 12))
        self.assertFalse(self._call_function_under_test(values, 13))

    def test_left_endpoint(self):
        values = [9, 11, 220]
        self.assertFalse(self._call_function_under_test(values, 7))
        self.assertFalse(self._call_function_under_test(values, 8))
        self.assertTrue(self._call_function_under_test(values, 9))
        self.assertFalse(self._call_function_under_test(values, 10))

    def test_right_endpoint(self):
        values = [16, 18, 20]
        self.assertFalse(self._call_function_under_test(values, 19))
        self.assertTrue(self._call_function_under_test(values, 20))
        self.assertFalse(self._call_function_under_test(values, 21))
        self.assertFalse(self._call_function_under_test(values, 22))


class Test_simple_convex_hull(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(points):
        from bezier import _py_helpers

        return _py_helpers.simple_convex_hull(points)

    def test_triangle_centroid(self):
        points = np.asfortranarray(
            [[0.0, 0.0, 1.0, 3.0, 0.0], [0.0, 3.0, 1.0, 0.0, 3.0]]
        )
        polygon = self._call_function_under_test(points)
        expected = np.asfortranarray([[0.0, 3.0, 0.0], [0.0, 0.0, 3.0]])
        self.assertEqual(expected, polygon)

    def test_two_points(self):
        points = np.asfortranarray([[0.0, 1.0], [0.0, 0.0]])
        polygon = self._call_function_under_test(points)
        expected = np.asfortranarray([[0.0, 1.0], [0.0, 0.0]])
        self.assertEqual(expected, polygon)
        # Switch the order of the points.
        points = np.asfortranarray([[1.0, 0.0], [0.0, 0.0]])
        polygon = self._call_function_under_test(points)
        self.assertEqual(expected, polygon)

    def test_one_point(self):
        points = np.asfortranarray([[2.0], [3.0]])
        polygon = self._call_function_under_test(points)
        expected = points
        self.assertEqual(expected, polygon)

    def test_zero_points(self):
        points = np.empty((2, 0), order="F")
        polygon = self._call_function_under_test(points)
        self.assertEqual(polygon.shape, (2, 0))

    def test_10x10_grid(self):
        points = np.empty((2, 100), order="F")
        index = 0
        for i in six.moves.xrange(10):
            for j in six.moves.xrange(10):
                points[:, index] = i, j
                index += 1
        polygon = self._call_function_under_test(points)
        expected = np.asfortranarray(
            [[0.0, 9.0, 9.0, 0.0], [0.0, 0.0, 9.0, 9.0]]
        )
        self.assertEqual(expected, polygon)

    def test_almost_linear(self):
        from bezier import _py_helpers

        # In a previous implementation, this case broke the algorithm
        # because the middle point of the line was placed in both the
        # upper and lower hull (which used 4 points for the hull when
        # only 3 were allocated).
        points = np.asfortranarray(
            [
                [
                    -0.12878911375710406,
                    -0.08626630936431968,
                    -0.043743504971535306,
                ],
                [
                    -0.05306646729159134,
                    -0.0032018988543520074,
                    0.04666266958288733,
                ],
            ]
        )
        polygon = self._call_function_under_test(points)
        expected = points
        self.assertEqual(expected, polygon)
        # Also verify why the case failed previously.
        point0 = points[:, 0]
        point1 = points[:, 1]
        point2 = points[:, 2]
        compare_lower = _py_helpers.cross_product_compare(
            point0, point1, point2
        )
        self.assertGreater(compare_lower, 0.0)
        compare_upper = _py_helpers.cross_product_compare(
            point2, point1, point0
        )
        self.assertGreater(compare_upper, 0.0)


class Test_is_separating(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(direction, polygon1, polygon2):
        from bezier import _py_helpers

        return _py_helpers.is_separating(direction, polygon1, polygon2)

    def test_it(self):
        direction = np.asfortranarray([0.0, 1.0])
        polygon1 = np.asfortranarray([[1.0, 3.0, -2.0], [1.0, 4.0, 3.0]])
        polygon2 = np.asfortranarray(
            [[3.5, 6.5, 6.5, 3.5], [3.0, 3.0, 7.0, 7.0]]
        )
        self.assertTrue(
            self._call_function_under_test(direction, polygon1, polygon2)
        )


class Test_polygon_collide(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(polygon1, polygon2):
        from bezier import _py_helpers

        return _py_helpers.polygon_collide(polygon1, polygon2)

    def test_first_edge(self):
        polygon1 = np.asfortranarray([[1.0, 3.0, -2.0], [1.0, 4.0, 3.0]])
        polygon2 = np.asfortranarray(
            [[3.5, 6.5, 6.5, 3.5], [3.0, 3.0, 7.0, 7.0]]
        )
        self.assertFalse(self._call_function_under_test(polygon1, polygon2))
        # Check with arguments swapped.
        self.assertFalse(self._call_function_under_test(polygon2, polygon1))

    def test_colliding(self):
        polygon1 = np.asfortranarray([[0.0, 3.0, 0.0], [0.0, 0.0, 3.0]])
        polygon2 = np.asfortranarray([[1.0, 4.0, 1.0], [1.0, 1.0, 4.0]])
        self.assertTrue(self._call_function_under_test(polygon1, polygon2))

    def test_non_first_edge_polygon1(self):
        polygon1 = np.asfortranarray(
            [[1.0, 1.0, 0.0, 0.0], [0.0, 1.0, 1.0, 0.0]]
        )
        polygon2 = np.asfortranarray(
            [[2.0, 3.0, 3.0, 2.0], [0.0, 0.0, 1.0, 1.0]]
        )
        self.assertFalse(self._call_function_under_test(polygon1, polygon2))

    def test_non_first_edge_polygon2(self):
        polygon1 = np.asfortranarray(
            [[0.0, 2.0, 2.0, 0.0], [0.0, 0.0, 2.0, 2.0]]
        )
        polygon2 = np.asfortranarray([[1.0, 4.0, 4.0], [4.0, 1.0, 4.0]])
        self.assertFalse(self._call_function_under_test(polygon1, polygon2))


class Test_solve2x2(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(lhs, rhs):
        from bezier import _py_helpers

        return _py_helpers.solve2x2(lhs, rhs)

    def test_solve_without_row_swap(self):
        lhs = np.asfortranarray([[2.0, 3.0], [1.0, 2.0]])
        rhs = np.asfortranarray([31.0, 19.0])
        singular, x_val, y_val = self._call_function_under_test(lhs, rhs)
        self.assertFalse(singular)
        self.assertEqual(x_val, 5.0)
        self.assertEqual(y_val, 7.0)

    def test_solve_with_row_swap(self):
        lhs = np.asfortranarray([[1.0, 0.0], [4.0, 1.0]])
        rhs = np.asfortranarray([3.0, 13.0])
        singular, x_val, y_val = self._call_function_under_test(lhs, rhs)
        self.assertFalse(singular)
        self.assertEqual(x_val, 3.0)
        self.assertEqual(y_val, 1.0)

    def test_zero_column(self):
        lhs = np.zeros((2, 2), order="F")
        singular, x_val, y_val = self._call_function_under_test(lhs, None)
        self.assertTrue(singular)
        self.assertIsNone(x_val)
        self.assertIsNone(y_val)

    def test_singular_without_row_swap(self):
        lhs = np.asfortranarray([[2.0, 4.0], [1.0, 2.0]])
        singular, x_val, y_val = self._call_function_under_test(lhs, None)
        self.assertTrue(singular)
        self.assertIsNone(x_val)
        self.assertIsNone(y_val)

    def test_singular_with_row_swap(self):
        lhs = np.asfortranarray([[3.0, 1.0], [12.0, 4.0]])
        singular, x_val, y_val = self._call_function_under_test(lhs, None)
        self.assertTrue(singular)
        self.assertIsNone(x_val)
        self.assertIsNone(y_val)


class TestUnsupportedDegree(unittest.TestCase):
    @staticmethod
    def _get_target_class():
        from bezier import _py_helpers

        return _py_helpers.UnsupportedDegree

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_inherit_not_implemented(self):
        klass = self._get_target_class()
        self.assertTrue(issubclass(klass, NotImplementedError))

    def test_constructor_defaults(self):
        exc = self._make_one(3)
        self.assertEqual(exc.degree, 3)
        self.assertEqual(exc.supported, ())

    def test_constructor_explicit(self):
        exc = self._make_one(4, supported=(1, 2))
        self.assertEqual(exc.degree, 4)
        self.assertEqual(exc.supported, (1, 2))

    def test___str__zero_supported(self):
        exc = self._make_one(1)
        as_str = str(exc)
        expected = "degree=1"
        self.assertEqual(as_str, expected)

    def test___str__one_supported(self):
        exc = self._make_one(2, supported=(1,))
        as_str = str(exc)
        expected = "The only degree supported at this time is 1 (degree=2)"
        self.assertEqual(as_str, expected)

    def test___str__multiple_supported(self):
        exc = self._make_one(3, supported=(1, 2))
        as_str = str(exc)
        expected = (
            "The only degrees supported at this time are 1 and 2 (degree=3)"
        )
        self.assertEqual(as_str, expected)
        exc = self._make_one(4, supported=(1, 3, 2))
        as_str = str(exc)
        expected = (
            "The only degrees supported at this "
            "time are 1, 3 and 2 (degree=4)"
        )
        self.assertEqual(as_str, expected)
