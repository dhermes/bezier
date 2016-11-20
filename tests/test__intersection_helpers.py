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


class Test__vector_close(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(vec1, vec2):
        from bezier import _intersection_helpers

        return _intersection_helpers._vector_close(vec1, vec2)

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

    def test_near_zero(self):
        vec1 = np.array([0.0, 0.0])
        vec2 = np.array([3.0, 4.0]) / 2.0**45
        self.assertTrue(self._call_function_under_test(vec1, vec2))

    def test_near_zero_fail(self):
        vec1 = np.array([1.0, 0.0]) / 2.0**20
        vec2 = np.array([0.0, 0.0])
        self.assertFalse(self._call_function_under_test(vec1, vec2))


class Test__check_close(unittest.TestCase):

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
        self.assertTrue(np.all(result == expected))

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
        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(0.0, 0.5)

    def test_near_endpoint(self):
        with self.assertRaises(NotImplementedError):
            self._call_function_under_test(0.75, 1.0 + 2.0**(-20))

    def test_s_outside(self):
        with self.assertRaises(ValueError):
            self._call_function_under_test(-0.25, 0.5)

    def test_t_outside(self):
        with self.assertRaises(ValueError):
            self._call_function_under_test(0.25, 1.5)

    def test_valid(self):
        # We really just want to make sure it doesn't raise.
        self.assertIsNone(self._call_function_under_test(0.25, 0.5))


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
        nodes = self.UNIT_SQUARE + np.array([[0.5, 0.5]])
        self.assertTrue(self._call_function_under_test(
            self.UNIT_SQUARE, nodes))

    def test_far_apart(self):
        nodes = self.UNIT_SQUARE + np.array([[100.0, 100.0]])
        self.assertFalse(self._call_function_under_test(
            self.UNIT_SQUARE, nodes))

    def test_tangent(self):
        nodes = self.UNIT_SQUARE + np.array([[1.0, 0.0]])
        self.assertFalse(self._call_function_under_test(
            self.UNIT_SQUARE, nodes))

    def test_almost_tangent(self):
        x_val = 1.0 + np.spacing(1.0)  # pylint: disable=no-member
        nodes = self.UNIT_SQUARE + np.array([[x_val, 0.0]])
        self.assertFalse(self._call_function_under_test(
            self.UNIT_SQUARE, nodes))


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


class Test__evaluate_hodograph(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, s):
        from bezier import _intersection_helpers

        return _intersection_helpers._evaluate_hodograph(nodes, degree, s)

    def test_line(self):
        degree = 1
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])

        first_deriv1 = self._call_function_under_test(nodes, degree, 0.25)
        self.assertEqual(first_deriv1.ndim, 1)
        self.assertTrue(np.all(first_deriv1 == nodes[1, :] - nodes[0, :]))
        # Make sure it is the same elsewhere since
        # the derivative curve is degree 0.
        first_deriv2 = self._call_function_under_test(nodes, degree, 0.75)
        self.assertTrue(np.all(first_deriv1 == first_deriv2))

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
            self.assertEqual(first_deriv.shape, (2,))
            self.assertEqual(first_deriv[0], (2.0 + s_val) / 2.0)
            self.assertEqual(first_deriv[1], (4.0 - 7.0 * s_val) / 2.0)

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
            self.assertEqual(first_deriv.shape, (2,))
            x_prime = 3.0 * (1.0 + 2.0 * s_val - s_val * s_val) / 4.0
            self.assertEqual(first_deriv[0], x_prime)
            y_prime = 3.0 * (5.0 * s_val * s_val - 6.0 * s_val + 2.0) / 2.0
            self.assertEqual(first_deriv[1], y_prime)


class Test_newton_refine(unittest.TestCase):

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
        self.assertTrue(np.all(
            curve1.evaluate(known_s) == curve2.evaluate(known_t)))

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
        self.assertTrue(np.all(
            curve1.evaluate(known_s) == curve2.evaluate(known_t)))

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

        self.assertTrue(np.all(
            curve1.evaluate(known_s) == curve2.evaluate(known_t)))
        new_s, new_t = self._call_function_under_test(
            known_s, curve1, known_t, curve2)

        self.assertEqual(new_s, known_s)
        self.assertEqual(new_t, known_t)

    def test_quadratic(self):
        curve1, curve2 = self._get_quadratics()

        known_s = 0.25
        known_t = 0.75
        self.assertTrue(np.all(
            curve1.evaluate(known_s) == curve2.evaluate(known_t)))

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
        self.assertTrue(np.all(parameters == expected))
        # Make sure that we've actually converged.
        exact_s, exact_t = parameters[-1, :]
        self.assertTrue(np.all(
            curve1.evaluate(exact_s) == curve2.evaluate(exact_t)))


class Test__cross_product(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(vec0, vec1):
        from bezier import _intersection_helpers

        return _intersection_helpers._cross_product(vec0, vec1)

    def test_it(self):
        self.assertTrue(False)


class Test_segment_intersection(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(start0, end0, start1, end1):
        from bezier import _intersection_helpers

        return _intersection_helpers.segment_intersection(
            start0, end0, start1, end1)

    def test_it(self):
        self.assertTrue(False)


class Test_from_linearized(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(linearized_pairs):
        from bezier import _intersection_helpers

        return _intersection_helpers.from_linearized(linearized_pairs)

    def test_it(self):
        self.assertTrue(False)


class Test_intersect_one_round(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(candidates):
        from bezier import _intersection_helpers

        return _intersection_helpers.intersect_one_round(candidates)

    def test_it(self):
        self.assertTrue(False)


class Test_all_intersections(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(candidates):
        from bezier import _intersection_helpers

        return _intersection_helpers.all_intersections(candidates)

    def test_it(self):
        self.assertTrue(False)


class TestLinearization(unittest.TestCase):

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

    def test__nodes_property(self):
        import bezier

        curve = bezier.Curve(self.NODES, _copy=False)
        linearization = self._make_one(curve)
        self.assertIs(linearization._nodes, self.NODES)

    def test_start_property(self):
        import bezier

        curve = bezier.Curve(self.NODES, _copy=False)
        linearization = self._make_one(curve)
        self.assertTrue(np.all(linearization.start == self.NODES[[0], :]))

    def test_end_property(self):
        import bezier

        curve = bezier.Curve(self.NODES, _copy=False)
        linearization = self._make_one(curve)
        self.assertTrue(np.all(linearization.end == self.NODES[[2], :]))

    def test_from_shape_factory_not_close_enough(self):
        import bezier

        curve = bezier.Curve(self.NODES, _copy=False)
        klass = self._get_target_class()
        new_shape, new_error = klass.from_shape(curve)
        self.assertIs(new_shape, curve)
        # NODES has constant second derivative equal to 2 * [3.0, 4.0].
        self.assertEqual(new_error, 0.125 * 2 * 1 * 5.0)

    def test_from_shape_factory_close_enough(self):
        import bezier

        scale_factor = 2.0**(-27)
        nodes = self.NODES * scale_factor
        curve = bezier.Curve(nodes, _copy=False)
        klass = self._get_target_class()
        new_shape, new_error = klass.from_shape(curve)

        self.assertIsInstance(new_shape, klass)
        self.assertIs(new_shape._curve, curve)
        # NODES has constant second derivative equal to 2 * [3.0, 4.0].
        expected_error = 0.125 * 2 * 1 * 5.0 * scale_factor
        self.assertEqual(new_error, expected_error)

    def test_from_shape_factory_no_error(self):
        import bezier

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve = bezier.Curve(nodes, _copy=False)
        klass = self._get_target_class()
        new_shape, new_error = klass.from_shape(curve)
        self.assertIsInstance(new_shape, klass)
        self.assertIs(new_shape._curve, curve)
        # ``nodes`` is linear, so error is 0.0.
        self.assertEqual(new_error, 0.0)

    def test_from_shape_factory_already_linearized(self):
        error = 0.078125
        linearization = self._make_one(None, error=error)

        klass = self._get_target_class()
        new_shape, new_error = klass.from_shape(linearization)
        self.assertIs(new_shape, linearization)
        self.assertEqual(new_error, error)
