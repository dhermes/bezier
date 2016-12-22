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


# pylint: disable=too-many-arguments
def check_intersection(test_case, intersection, expected,
                       curve1, curve2, s_val, t_val):
    from bezier import _intersection_helpers

    test_case.assertIsInstance(
        intersection, _intersection_helpers.Intersection)
    test_case.assertEqual(intersection.point, expected)
    test_case.assertIs(intersection.left, curve1)
    test_case.assertEqual(intersection.s, s_val)
    test_case.assertIs(intersection.right, curve2)
    test_case.assertEqual(intersection.t, t_val)
# pylint: enable=too-many-arguments


class Test__check_close(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(s, curve1, t, curve2):
        from bezier import _intersection_helpers

        return _intersection_helpers._check_close(
            s, curve1, t, curve2)

    def test_success(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ])
        curve = bezier.Curve(nodes)
        s_val = 0.5
        wiggle = 2.0**(-50)
        result = self._call_function_under_test(
            s_val, curve, s_val + wiggle, curve)

        expected = np.array([0.5, 0.5])
        self.assertEqual(result, expected)

    def test_failure(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve = bezier.Curve(nodes)
        # The nodes of thise curve are far away.
        with self.assertRaises(ValueError):
            self._call_function_under_test(0.0, curve, 1.0, curve)


class Test__check_parameters(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(s, t):
        from bezier import _intersection_helpers

        return _intersection_helpers._check_parameters(s, t)

    def test_at_endpoint(self):
        # Really just making sure the function doesn't raise.
        result = self._call_function_under_test(0.0, 0.5)
        self.assertEqual(result, (0.0, 0.5))
        result = self._call_function_under_test(0.5, 1.0)
        self.assertEqual(result, (0.5, 1.0))

    def test_near_endpoint(self):
        with self.assertRaises(ValueError):
            self._call_function_under_test(0.75, 1.0 + 2.0**(-20))

    def test_s_outside(self):
        with self.assertRaises(ValueError):
            self._call_function_under_test(-0.25, 0.5)

    def test_t_outside(self):
        with self.assertRaises(ValueError):
            self._call_function_under_test(0.25, 1.5)

    def test_valid(self):
        # Really just making sure the function doesn't raise.
        result = self._call_function_under_test(0.25, 0.5)
        self.assertEqual(result, (0.25, 0.5))

    def test_wiggle_below(self):
        value = -2.0**(-60)
        result = self._call_function_under_test(value, 0.25)
        self.assertEqual(result, (0.0, 0.25))

    def test_wiggle_above(self):
        value = 1 + 2.0**(-52)
        result = self._call_function_under_test(0.875, value)
        self.assertEqual(result, (0.875, 1.0))


class Test_bbox_intersect(unittest.TestCase):

    UNIT_SQUARE = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ])

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _intersection_helpers

        return _intersection_helpers.bbox_intersect(nodes1, nodes2)

    def test_intersect(self):
        from bezier import _intersection_helpers

        nodes = self.UNIT_SQUARE + np.array([[0.5, 0.5]])
        result = self._call_function_under_test(self.UNIT_SQUARE, nodes)
        self.assertIs(
            result, _intersection_helpers.BoxIntersectionType.intersection)

    def test_far_apart(self):
        from bezier import _intersection_helpers

        nodes = self.UNIT_SQUARE + np.array([[100.0, 100.0]])
        result = self._call_function_under_test(self.UNIT_SQUARE, nodes)
        self.assertIs(
            result, _intersection_helpers.BoxIntersectionType.disjoint)

    def test_disjoint_but_aligned(self):
        from bezier import _intersection_helpers

        nodes = self.UNIT_SQUARE + np.array([[1.0, 2.0]])
        result = self._call_function_under_test(self.UNIT_SQUARE, nodes)
        self.assertIs(
            result, _intersection_helpers.BoxIntersectionType.disjoint)

    def test_tangent(self):
        from bezier import _intersection_helpers

        nodes = self.UNIT_SQUARE + np.array([[1.0, 0.0]])
        result = self._call_function_under_test(self.UNIT_SQUARE, nodes)
        self.assertIs(
            result, _intersection_helpers.BoxIntersectionType.tangent)

    def test_almost_tangent(self):
        from bezier import _intersection_helpers

        x_val = 1.0 + np.spacing(1.0)  # pylint: disable=no-member
        nodes = self.UNIT_SQUARE + np.array([[x_val, 0.0]])
        result = self._call_function_under_test(self.UNIT_SQUARE, nodes)
        self.assertIs(
            result, _intersection_helpers.BoxIntersectionType.disjoint)


class Test_linearization_error(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(curve):
        from bezier import _intersection_helpers

        return _intersection_helpers.linearization_error(curve)

    def test_linear(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        self.assertEqual(error_val, 0.0)

    def test_degree_elevated_linear(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 2.0],
        ])
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        self.assertEqual(error_val, 0.0)

        nodes = np.array([
            [0.0, 0.0],
            [0.25, 0.5],
            [0.5, 1.0],
            [0.75, 1.5],
            [1.0, 2.0],
        ])
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        self.assertEqual(error_val, 0.0)

    def test_hidden_linear(self):
        import bezier

        # NOTE: This is the line 3 y = 4 x, but with the parameterization
        #       x(s) = 3 s (4 - 3 s).
        nodes = np.array([
            [0.0, 0.0],
            [6.0, 8.0],
            [3.0, 4.0],
        ])
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        # D^2 v = [-9, -12]
        expected = 0.125 * 2 * 1 * 15.0
        self.assertEqual(error_val, expected)

    def test_quadratic(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
            [5.0, 6.0],
        ])
        # NOTE: This is hand picked so that
        #             d Nodes = [1, 1], [4, 5]
        #           d^2 Nodes = [3, 4]
        #       so that sqrt(3^2 + 4^2) = 5.0
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        expected = 0.125 * 2 * 1 * 5.0
        self.assertEqual(error_val, expected)

        # For a degree two curve, the 2nd derivative is constant
        # so by subdividing, our error should drop by a factor
        # of (1/2)^2 = 4.
        left, right = curve.subdivide()
        error_left = self._call_function_under_test(left)
        error_right = self._call_function_under_test(right)
        self.assertEqual(error_left, 0.25 * expected)
        self.assertEqual(error_right, 0.25 * expected)

    def test_hidden_quadratic(self):
        import bezier

        # NOTE: This is the line y = 1 + x^2 / 4, but with the
        #       parameterization x(s) = (3 s - 1)^2.
        nodes = np.array([
            [1.0, 1.25],
            [-0.5, 0.5],
            [-0.5, 2.0],
            [1.0, -1.0],
            [4.0, 5.0],
        ])
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        # D^2 v = [1.5, 2.25], [1.5, -4.5], [1.5, 9]
        expected = 0.125 * 4 * 3 * np.sqrt(1.5**2 + 9.0**2)
        self.assertEqual(error_val, expected)

    def test_cubic(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
            [5.0, 6.0],
            [6.0, 7.0],
        ])
        # NOTE: This is hand picked so that
        #             d Nodes = [1, 1], [4, 5], [1, 1]
        #           d^2 Nodes = [3, 4], [-3, -4]
        #       so that sqrt(3^2 + 4^2) = 5.0
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        expected = 0.125 * 3 * 2 * 5.0
        self.assertEqual(error_val, expected)

    def test_quartic(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
            [5.0, 6.0],
            [6.0, 7.0],
            [4.0, 7.0],
        ])
        # NOTE: This is hand picked so that
        #             d Nodes = [1, 1], [4, 5], [1, 1], [-2, 0]
        #           d^2 Nodes = [3, 4], [-3, -4], [-3, -1]
        #       so that sqrt(3^2 + 4^2) = 5.0
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        expected = 0.125 * 4 * 3 * 5.0
        self.assertEqual(error_val, expected)

    def test_degree_weights_on_the_fly(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
            [7.0, 3.0],
            [11.0, 8.0],
            [15.0, 1.0],
            [16.0, -3.0],
        ])
        # NOTE: This is hand picked so that
        #             d Nodes = [1, 1], [6, 2], [4, 5], [4, -7], [1, -4]
        #           d^2 Nodes = [5, 1], [-2, 3], [0, -12], [-3, 3]
        #       so that sqrt(5^2 + 12^2) = 13.0
        curve = bezier.Curve(nodes)
        error_val = self._call_function_under_test(curve)
        expected = 0.125 * 5 * 4 * 13.0
        self.assertEqual(error_val, expected)


class Test_newton_refine(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(s, curve1, t, curve2):
        from bezier import _intersection_helpers

        return _intersection_helpers.newton_refine(s, curve1, t, curve2)

    def test_linear(self):
        import bezier

        nodes1 = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        nodes2 = np.array([
            [1.0, 0.0],
            [0.0, 3.0],
        ])
        curve1 = bezier.Curve(nodes1)
        curve2 = bezier.Curve(nodes2)

        known_s = 0.75
        known_t = 0.25
        self.assertEqual(curve1.evaluate(known_s),
                         curve2.evaluate(known_t))

        wrong_s = known_s - 0.125
        wrong_t = known_t + 0.125
        # NOTE: By construction, the Jacobian matrix will be
        #           [1, 1], [1, -3]
        #       which has determinant -4.0, hence there will
        #       be no round-off when solving.
        new_s, new_t = self._call_function_under_test(
            wrong_s, curve1, wrong_t, curve2)

        # Newton's method is exact on linear problems so will
        # always converge after one step.
        self.assertEqual(new_s, known_s)
        self.assertEqual(new_t, known_t)

    @staticmethod
    def _get_quadratics():
        import bezier

        nodes1 = np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ])
        nodes2 = np.array([
            [1.0, 0.75],
            [0.5, -0.25],
            [0.0, 0.75],
        ])
        curve1 = bezier.Curve(nodes1)
        curve2 = bezier.Curve(nodes2)
        return curve1, curve2

    def test_mixed_degree(self):
        import bezier

        curve1, _ = self._get_quadratics()
        nodes2 = np.array([
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        curve2 = bezier.Curve(nodes2)

        known_s = 0.5
        known_t = 0.5
        self.assertEqual(curve1.evaluate(known_s),
                         curve2.evaluate(known_t))

        wrong_s = 0.25
        wrong_t = 0.25
        # NOTE: By construction, the Jacobian matrix will be
        #           [1, 1], [1, -1]
        #       which has determinant -2.0, hence there will
        #       be no round-off when solving.
        new_s, new_t = self._call_function_under_test(
            wrong_s, curve1, wrong_t, curve2)

        self.assertEqual(new_s, 0.4375)
        self.assertEqual(new_t, 0.5625)

        # Make sure we have gotten closer to correct.
        self.assertLess(abs(known_s - new_s), abs(known_s - wrong_s))
        self.assertLess(abs(known_t - new_t), abs(known_t - wrong_t))

    def test_early_exit(self):
        curve1, curve2 = self._get_quadratics()

        known_s = 0.25
        known_t = 0.75

        self.assertEqual(curve1.evaluate(known_s),
                         curve2.evaluate(known_t))

        new_s, new_t = self._call_function_under_test(
            known_s, curve1, known_t, curve2)

        self.assertEqual(new_s, known_s)
        self.assertEqual(new_t, known_t)

    def test_quadratic(self):
        curve1, curve2 = self._get_quadratics()

        known_s = 0.25
        known_t = 0.75
        self.assertEqual(curve1.evaluate(known_s),
                         curve2.evaluate(known_t))

        wrong_s = known_s + 0.0625  # 1/16
        wrong_t = known_t + 0.0625  # 1/16
        # NOTE: By construction, the Jacobian matrix will be
        #           [1, 3/4], [1, -5/4]
        #       which has determinant -2.0, hence there will
        #       be no round-off when solving.
        new_s, new_t = self._call_function_under_test(
            wrong_s, curve1, wrong_t, curve2)

        self.assertEqual(new_s, 0.2421875)
        self.assertEqual(new_t, 0.7578125)

        # Make sure we have gotten closer to correct.
        self.assertLess(abs(known_s - new_s), abs(known_s - wrong_s))
        self.assertLess(abs(known_t - new_t), abs(known_t - wrong_t))

    def test_convergence(self):
        import six

        import bezier

        nodes1 = np.array([
            [0.0, 0.0],
            [0.25, 1.0],
            [0.5, -0.75],
            [0.75, 1.0],
            [1.0, 0.0],
        ])
        curve1 = bezier.Curve(nodes1)
        # Vertical line forces a unique solution.
        nodes2 = np.array([
            [0.5, 0.0],
            [0.5, 1.0],
        ])
        curve2 = bezier.Curve(nodes2)

        num_guess = 4
        parameters = np.zeros((num_guess, 2))
        # NOTE: This means our "first" guess is (s, t) = (0, 0).
        for guess in six.moves.xrange(1, num_guess):
            prev_s, prev_t = parameters[guess - 1, :]
            parameters[guess, :] = self._call_function_under_test(
                prev_s, curve1, prev_t, curve2)

        expected = np.array([
            [0.0, 0.0],
            [0.5, 2.0],
            [0.5, 0.21875],
            [0.5, 0.21875],
        ])
        self.assertEqual(parameters, expected)
        # Make sure that we've actually converged.
        exact_s, exact_t = parameters[-1, :]
        self.assertEqual(curve1.evaluate(exact_s),
                         curve2.evaluate(exact_t))


class Test_segment_intersection(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(start0, end0, start1, end1, **kwargs):
        from bezier import _intersection_helpers

        return _intersection_helpers.segment_intersection(
            start0, end0, start1, end1, **kwargs)

    def _helper(self, intersection, s_val, direction0,
                t_val, direction1, **kwargs):
        start0 = intersection + s_val * direction0
        end0 = intersection + (s_val - 1.0) * direction0
        start1 = intersection + t_val * direction1
        end1 = intersection + (t_val - 1.0) * direction1

        return self._call_function_under_test(
            start0, end0, start1, end1, **kwargs)

    def test_success(self):
        intersection = np.array([[1.0, 2.0]])
        s_val = 0.25
        t_val = 0.625
        direction0 = np.array([[3.0, 0.5]])
        direction1 = np.array([[-2.0, 1.0]])
        # D0 x D1 == 4.0, so there will be no round-off in answer.
        computed_s, computed_t, success = self._helper(
            intersection, s_val, direction0, t_val, direction1)

        self.assertEqual(computed_s, s_val)
        self.assertEqual(computed_t, t_val)
        self.assertTrue(success)

    def test_parallel(self):
        intersection = np.array([[0.0, 0.0]])
        s_val = 0.5
        t_val = 0.5
        direction0 = np.array([[0.0, 1.0]])
        direction1 = np.array([[0.0, 2.0]])
        computed_s, computed_t, success = self._helper(
            intersection, s_val,
            direction0, t_val, direction1)

        self.assertIsNone(computed_s)
        self.assertIsNone(computed_t)
        self.assertFalse(success)


class Test_parallel_different(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(start0, end0, start1, end1):
        from bezier import _intersection_helpers

        return _intersection_helpers.parallel_different(
            start0, end0, start1, end1)

    def test_same_line_no_overlap(self):
        start0 = np.array([[0.0, 0.0]])
        end0 = np.array([[3.0, 4.0]])
        start1 = np.array([[6.0, 8.0]])
        end1 = np.array([[9.0, 12.0]])
        self.assertTrue(
            self._call_function_under_test(start0, end0, start1, end1))

    def test_same_line_overlap_at_start(self):
        start0 = np.array([[6.0, -3.0]])
        end0 = np.array([[-7.0, 1.0]])
        start1 = np.array([[1.125, -1.5]])
        end1 = np.array([[-5.375, 0.5]])
        self.assertFalse(
            self._call_function_under_test(start0, end0, start1, end1))

    def test_same_line_overlap_at_end(self):
        start0 = np.array([[1.0, 2.0]])
        end0 = np.array([[3.0, 5.0]])
        start1 = np.array([[-0.5, -0.25]])
        end1 = np.array([[2.0, 3.5]])
        self.assertFalse(
            self._call_function_under_test(start0, end0, start1, end1))

    def test_same_line_contained(self):
        start0 = np.array([[-9.0, 0.0]])
        end0 = np.array([[4.0, 5.0]])
        start1 = np.array([[23.5, 12.5]])
        end1 = np.array([[-25.25, -6.25]])
        self.assertFalse(
            self._call_function_under_test(start0, end0, start1, end1))

    def test_different_line(self):
        start0 = np.array([[3.0, 2.0]])
        end0 = np.array([[3.0, 0.75]])
        start1 = np.array([[0.0, 0.0]])
        end1 = np.array([[0.0, 2.0]])
        self.assertTrue(
            self._call_function_under_test(start0, end0, start1, end1))


class Test_from_linearized(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(left, right, intersections):
        from bezier import _intersection_helpers

        return _intersection_helpers.from_linearized(
            left, right, intersections)

    def test_it(self):
        import bezier
        from bezier import _intersection_helpers

        nodes1 = np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 1.0],
        ])
        curve1 = bezier.Curve(nodes1)
        # NOTE: This curve isn't close to linear, but that's OK.
        lin1 = _intersection_helpers.Linearization(curve1)

        nodes2 = np.array([
            [0.0, 1.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ])
        curve2 = bezier.Curve(nodes2)
        # NOTE: This curve isn't close to linear, but that's OK.
        lin2 = _intersection_helpers.Linearization(curve2)

        intersections = []
        self.assertIsNone(
            self._call_function_under_test(lin1, lin2, intersections))
        self.assertEqual(len(intersections), 1)
        intersection = intersections[0]
        expected = curve1.evaluate(0.5)
        check_intersection(self, intersection, expected,
                           curve1, curve2, 0.5, 0.5)

    def test_no_intersection(self):
        import bezier
        from bezier import _intersection_helpers
        # The bounding boxes intersect but the lines do not.

        nodes1 = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve1 = bezier.Curve(nodes1)
        lin1 = _intersection_helpers.Linearization(curve1)

        nodes2 = np.array([
            [1.75, -0.75],
            [0.75, 0.25],
        ])
        curve2 = bezier.Curve(nodes2)
        lin2 = _intersection_helpers.Linearization(curve2)

        intersections = []
        self.assertIsNone(
            self._call_function_under_test(lin1, lin2, intersections))
        self.assertEqual(intersections, [])

    def test_parallel_intersection(self):
        import bezier
        from bezier import _intersection_helpers

        nodes1 = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve1 = bezier.Curve(nodes1)
        lin1 = _intersection_helpers.Linearization(curve1)

        nodes2 = np.array([
            [0.0, 1.0],
            [1.0, 2.0],
        ])
        curve2 = bezier.Curve(nodes2)
        lin2 = _intersection_helpers.Linearization(curve2)

        intersections = []
        self.assertIsNone(
            self._call_function_under_test(lin1, lin2, intersections))
        self.assertEqual(intersections, [])

    def test_same_line_intersection(self):
        import bezier
        from bezier import _intersection_helpers

        nodes1 = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve1 = bezier.Curve(nodes1)
        lin1 = _intersection_helpers.Linearization(curve1)

        nodes2 = np.array([
            [0.5, 0.5],
            [3.0, 3.0],
        ])
        curve2 = bezier.Curve(nodes2)
        lin2 = _intersection_helpers.Linearization(curve2)

        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(lin1, lin2, [])

    def test_parallel_non_degree_one(self):
        import bezier
        from bezier import _intersection_helpers

        nodes1 = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve1 = bezier.Curve(nodes1)
        lin1 = _intersection_helpers.Linearization(curve1)

        nodes2 = np.array([
            [2.0, 2.0],
            [2.5009765625, 2.5009765625],
            [3.0, 3.0],
        ])
        curve2 = bezier.Curve(nodes2)
        lin2 = _intersection_helpers.Linearization(curve2)

        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(lin1, lin2, [])


class Test__add_intersection(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection, intersections):
        from bezier import _intersection_helpers

        return _intersection_helpers._add_intersection(
            intersection, intersections)

    def test_new(self):
        from bezier import _intersection_helpers

        intersection1 = _intersection_helpers.Intersection(
            mock.sentinel.left, 0.5, mock.sentinel.right, 0.5)
        intersection2 = _intersection_helpers.Intersection(
            mock.sentinel.left, 0.75, mock.sentinel.right, 0.25)
        intersections = [intersection1]
        self.assertIsNone(
            self._call_function_under_test(intersection2, intersections))

        self.assertEqual(intersections, [intersection1, intersection2])

    def test_existing(self):
        from bezier import _intersection_helpers

        intersection1 = _intersection_helpers.Intersection(
            mock.sentinel.left, 0.0, mock.sentinel.right, 1.0)
        intersection2 = _intersection_helpers.Intersection(
            mock.sentinel.left, 0.0, mock.sentinel.right, 1.0)
        intersections = [intersection1]
        self.assertIsNone(
            self._call_function_under_test(intersection2, intersections))

        self.assertEqual(len(intersections), 1)
        self.assertIs(intersections[0], intersection1)
        self.assertIsNot(intersections[0], intersection2)


class Test__tangent_bbox_intersection(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(left, right, intersections):
        from bezier import _intersection_helpers

        return _intersection_helpers._tangent_bbox_intersection(
            left, right, intersections)

    def test_it(self):
        import bezier

        nodes1 = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
            [2.0, 0.0],
        ])
        curve1 = bezier.Curve(nodes1)
        nodes2 = np.array([
            [2.0, 0.0],
            [3.0, 2.0],
            [4.0, 0.0],
        ])
        curve2 = bezier.Curve(nodes2)

        intersections = []
        self.assertIsNone(
            self._call_function_under_test(curve1, curve2, intersections))
        self.assertEqual(len(intersections), 1)
        intersection = intersections[0]
        expected = curve1.evaluate(1.0)
        check_intersection(self, intersection, expected,
                           curve1, curve2, 1.0, 0.0)


class Test_intersect_one_round(utils.NumPyTestCase):

    # NOTE: NODES1 is a specialization of [0, 0], [1/2, 1], [1, 1]
    #       onto the interval [1/4, 1].
    NODES1 = np.array([
        [0.25, 0.4375],
        [0.625, 1.0],
        [1.0, 1.0],
    ])
    # NOTE: NODES2 is a specialization of [0, 1], [1/2, 1], [1, 0]
    #       onto the interval [0, 3/4].
    NODES2 = np.array([
        [0.0, 1.0],
        [0.375, 1.0],
        [0.75, 0.4375],
    ])
    LINE1 = np.array([
        [0.0, 0.0],
        [1.0, 1.0],
    ])
    LINE2 = np.array([
        [0.0, 1.0],
        [1.0, 0.0],
    ])

    @staticmethod
    def _call_function_under_test(candidates, intersections):
        from bezier import _intersection_helpers

        return _intersection_helpers.intersect_one_round(
            candidates, intersections)

    def test_simple(self):
        import bezier

        curve1 = bezier.Curve(self.NODES1)
        curve2 = bezier.Curve(self.NODES2)
        candidates = [(curve1, curve2)]
        accepted = self._call_function_under_test(
            candidates, [])

        self.assertEqual(accepted, [(curve1, curve2)])

    def test_left_linearized(self):
        import bezier
        from bezier import _intersection_helpers

        curve1 = bezier.Curve(self.LINE1)
        lin1 = _intersection_helpers.Linearization(curve1)
        curve2 = bezier.Curve(self.LINE2)

        intersections = []
        accepted = self._call_function_under_test(
            [(lin1, curve2)], intersections)

        self.assertEqual(intersections, [])
        self.assertEqual(accepted, [(lin1, curve2)])

    def test_right_linearized(self):
        import bezier
        from bezier import _intersection_helpers

        curve1 = bezier.Curve(self.LINE1)
        curve2 = bezier.Curve(self.LINE2)
        lin2 = _intersection_helpers.Linearization(curve2)

        intersections = []
        accepted = self._call_function_under_test(
            [(curve1, lin2)], intersections)

        self.assertEqual(intersections, [])
        self.assertEqual(accepted, [(curve1, lin2)])

    def test_both_linearized(self):
        import bezier
        from bezier import _intersection_helpers

        curve1 = bezier.Curve(self.LINE1)
        lin1 = _intersection_helpers.Linearization(curve1)
        curve2 = bezier.Curve(self.LINE2)
        lin2 = _intersection_helpers.Linearization(curve2)

        intersections = []
        accepted = self._call_function_under_test(
            [(lin1, lin2)], intersections)
        self.assertEqual(accepted, [])
        self.assertEqual(len(intersections), 1)
        intersection = intersections[0]
        expected = np.array([0.5, 0.5])
        check_intersection(self, intersection, expected,
                           curve1, curve2, 0.5, 0.5)


class Test__next_candidates(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(left, right):
        from bezier import _intersection_helpers

        return _intersection_helpers._next_candidates(left, right)

    def test_it(self):
        import types
        import bezier
        from bezier import _intersection_helpers

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve = bezier.Curve(nodes)
        lin = _intersection_helpers.Linearization(curve)

        result = self._call_function_under_test(lin, lin)

        self.assertIsInstance(result, types.GeneratorType)
        pairs = list(result)
        self.assertEqual(pairs, [(lin, lin)])


class Test_all_intersections(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(candidates):
        from bezier import _intersection_helpers

        return _intersection_helpers.all_intersections(candidates)

    def test_failure(self):
        patch = mock.patch(
            'bezier._intersection_helpers._MAX_INTERSECT_SUBDIVISIONS',
            new=-1)
        with patch:
            with self.assertRaises(ValueError):
                self._call_function_under_test([])

    def test_no_intersections(self):
        intersections = self._call_function_under_test([])
        self.assertEqual(intersections, [])

    def test_tangent(self):
        import itertools
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ])
        curve = bezier.Curve(nodes)

        # Start with two pieces.
        pieces = itertools.chain(*[curve.subdivide()])
        # Double the number of pieces 3 times.
        for unused_exponent in (2, 3, 4):
            pieces = itertools.chain(*[
                piece.subdivide() for piece in pieces])

        # Match each piece with itself.
        candidates = [(piece, piece) for piece in pieces]
        # But we need **more** than 16.
        candidates.append((curve, curve))

        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(candidates)

    def test_success(self):
        import bezier

        # NOTE: ``nodes1`` is a specialization of [0, 0], [1/2, 1], [1, 1]
        #       onto the interval [1/4, 1] and ``nodes`` is a specialization
        #       of [0, 1], [1/2, 1], [1, 0] onto the interval [0, 3/4].
        #       We expect them to intersect at s = 1/3, t = 2/3, which is
        #       the point [1/2, 3/4].
        nodes1 = np.array([
            [0.25, 0.4375],
            [0.625, 1.0],
            [1.0, 1.0],
        ])
        curve1 = bezier.Curve(nodes1)

        nodes2 = np.array([
            [0.0, 1.0],
            [0.375, 1.0],
            [0.75, 0.4375],
        ])
        curve2 = bezier.Curve(nodes2)

        candidates = [(curve1, curve2)]
        intersections = self._call_function_under_test(candidates)

        self.assertEqual(len(intersections), 1)
        intersection = intersections[0]
        expected = np.array([0.5, 0.75])

        s_val = 1.0 / 3.0
        # Due to round-off, the answer is wrong by a tiny wiggle.
        s_val += np.spacing(s_val)  # pylint: disable=no-member
        t_val = 2.0 / 3.0
        check_intersection(self, intersection, expected,
                           curve1, curve2, s_val, t_val)


class TestLinearization(utils.NumPyTestCase):

    NODES = np.array([
        [0.0, 0.0],
        [1.0, 1.0],
        [5.0, 6.0],
    ])

    @staticmethod
    def _get_target_class():
        from bezier import _intersection_helpers

        return _intersection_helpers.Linearization

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        linearization = self._make_one(mock.sentinel.curve)
        self.assertIs(linearization._curve, mock.sentinel.curve)
        self.assertIsNone(linearization._error)

    def test_constructor_with_error(self):
        error = 0.125
        linearization = self._make_one(mock.sentinel.curve, error=error)
        self.assertIs(linearization._curve, mock.sentinel.curve)
        self.assertEqual(linearization._error, error)

    def test_subdivide(self):
        linearization = self._make_one(None)
        self.assertEqual(linearization.subdivide(), (linearization,))

    def test_curve_property(self):
        linearization = self._make_one(mock.sentinel.curve)
        self.assertIs(linearization.curve, mock.sentinel.curve)

    def test_error_property(self):
        error = 0.0625
        linearization = self._make_one(None, error=error)
        self.assertEqual(linearization.error, error)

    def test_error_property_on_the_fly(self):
        linearization = self._make_one(mock.sentinel.curve)

        error = 0.09375
        patch = mock.patch(
            'bezier._intersection_helpers.linearization_error',
            return_value=error)
        with patch as mocked:
            self.assertEqual(linearization.error, error)
            mocked.assert_called_once_with(mock.sentinel.curve)

    def test_start_node_property(self):
        import bezier

        curve = bezier.Curve(self.NODES, _copy=False)
        linearization = self._make_one(curve)
        expected = self.NODES[[0], :]
        self.assertEqual(linearization.start_node, expected)

    def test_end_node_property(self):
        import bezier

        curve = bezier.Curve(self.NODES, _copy=False)
        linearization = self._make_one(curve)
        expected = self.NODES[[2], :]
        self.assertEqual(linearization.end_node, expected)

    def test_from_shape_factory_not_close_enough(self):
        import bezier

        curve = bezier.Curve(self.NODES, _copy=False)
        klass = self._get_target_class()
        new_shape = klass.from_shape(curve)
        self.assertIs(new_shape, curve)

    def test_from_shape_factory_close_enough(self):
        import bezier

        scale_factor = 2.0**(-27)
        nodes = self.NODES * scale_factor
        curve = bezier.Curve(nodes, _copy=False)
        klass = self._get_target_class()
        new_shape = klass.from_shape(curve)

        self.assertIsInstance(new_shape, klass)
        self.assertIs(new_shape._curve, curve)
        # NODES has constant second derivative equal to 2 * [3.0, 4.0].
        expected_error = 0.125 * 2 * 1 * 5.0 * scale_factor
        self.assertEqual(new_shape.error, expected_error)

    def test_from_shape_factory_no_error(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve = bezier.Curve(nodes, _copy=False)
        klass = self._get_target_class()
        new_shape = klass.from_shape(curve)
        self.assertIsInstance(new_shape, klass)
        self.assertIs(new_shape._curve, curve)
        # ``nodes`` is linear, so error is 0.0.
        self.assertEqual(new_shape.error, 0.0)

    def test_from_shape_factory_already_linearized(self):
        error = 0.078125
        linearization = self._make_one(None, error=error)

        klass = self._get_target_class()
        new_shape = klass.from_shape(linearization)
        self.assertIs(new_shape, linearization)
        self.assertEqual(new_shape.error, error)


class TestIntersection(unittest.TestCase):

    @staticmethod
    def _get_target_class():
        from bezier import _intersection_helpers

        return _intersection_helpers.Intersection

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def _constructor_helper(self, **kwargs):
        left = mock.sentinel.left
        s_val = 0.25
        right = mock.sentinel.right
        t_val = 0.75

        intersection = self._make_one(
            left, s_val, right, t_val, **kwargs)

        self.assertIs(intersection._left, left)
        self.assertEqual(intersection._s_val, s_val)
        self.assertIs(intersection._right, right)
        self.assertEqual(intersection._t_val, t_val)
        self.assertIsNone(intersection._interior_curve)
        return intersection

    def test_constructor(self):
        intersection = self._constructor_helper()
        self.assertIsNone(intersection._point)

    def test_constructor_with_point(self):
        intersection = self._constructor_helper(point=mock.sentinel.point)
        self.assertIs(intersection._point, mock.sentinel.point)

    def test_s_property(self):
        s = 0.375
        intersection = self._make_one(None, s, None, None)
        self.assertEqual(intersection.s, s)

    def test_t_property(self):
        t = 0.9990234375
        intersection = self._make_one(None, None, None, t)
        self.assertEqual(intersection.t, t)

    def test_left_property(self):
        intersection = self._make_one(
            mock.sentinel.left, None, None, None)
        self.assertIs(intersection.left, mock.sentinel.left)

    def test_right_property(self):
        intersection = self._make_one(
            None, None, mock.sentinel.right, None)
        self.assertIs(intersection.right, mock.sentinel.right)

    def test_point_property(self):
        s_val = 1.0
        t_val = 0.0
        intersection = self._make_one(
            mock.sentinel.left, s_val, mock.sentinel.right, t_val)

        patch = mock.patch(
            'bezier._intersection_helpers._check_close',
            return_value=mock.sentinel.point)
        with patch as mocked:
            self.assertIsNone(intersection._point)
            self.assertIs(intersection.point, mock.sentinel.point)
            self.assertIs(intersection._point, mock.sentinel.point)

            self.assertEqual(mocked.call_count, 1)
            # Make sure the cached value is used on future access.
            self.assertIs(intersection.point, mock.sentinel.point)
            self.assertEqual(mocked.call_count, 1)

    def test_interior_curve_get_unset(self):
        intersection = self._make_one(None, None, None, None)
        self.assertIsNone(intersection._interior_curve)
        with self.assertRaises(AttributeError):
            getattr(intersection, 'interior_curve')

    def test_interior_curve_get_when_set(self):
        intersection = self._make_one(None, None, None, None)
        intersection._interior_curve = mock.sentinel.interior
        self.assertIs(intersection.interior_curve, mock.sentinel.interior)

    def test_interior_curve_set_when_set(self):
        intersection = self._make_one(None, None, None, None)
        intersection.interior_curve = -1
        with self.assertRaises(AttributeError):
            setattr(intersection, 'interior_curve', 0)

    def test_interior_curve_set_valid(self):
        intersection = self._make_one(None, None, None, None)
        intersection.interior_curve = 0
        self.assertEqual(intersection.interior_curve, 0)
