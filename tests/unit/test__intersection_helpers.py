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
import six

from tests.unit import utils


class Test__newton_refine(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(s, nodes1, t, nodes2):
        from bezier import _intersection_helpers

        return _intersection_helpers._newton_refine(s, nodes1, t, nodes2)

    def test_linear(self):
        import bezier

        nodes1 = np.asfortranarray([
            [0.0, 1.0],
            [0.0, 1.0],
        ])
        nodes2 = np.asfortranarray([
            [1.0, 0.0],
            [0.0, 3.0],
        ])
        curve1 = bezier.Curve(nodes1, degree=1)
        curve2 = bezier.Curve(nodes2, degree=1)

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
            wrong_s, nodes1, wrong_t, nodes2)

        # Newton's method is exact on linear problems so will
        # always converge after one step.
        self.assertEqual(new_s, known_s)
        self.assertEqual(new_t, known_t)

    @staticmethod
    def _get_quadratics():
        import bezier

        nodes1 = np.asfortranarray([
            [0.0, 0.5, 1.0],
            [0.0, 1.0, 0.0],
        ])
        nodes2 = np.asfortranarray([
            [1.0, 0.5, 0.0],
            [0.75, -0.25, 0.75],
        ])
        curve1 = bezier.Curve(nodes1, degree=2)
        curve2 = bezier.Curve(nodes2, degree=2)
        return curve1, curve2

    def test_mixed_degree(self):
        import bezier

        curve1, _ = self._get_quadratics()
        nodes2 = np.asfortranarray([
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        curve2 = bezier.Curve(nodes2, degree=1)

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
            wrong_s, curve1._nodes, wrong_t, nodes2)

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
            known_s, curve1._nodes, known_t, curve2._nodes)

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
            wrong_s, curve1._nodes, wrong_t, curve2._nodes)

        self.assertEqual(new_s, 0.2421875)
        self.assertEqual(new_t, 0.7578125)

        # Make sure we have gotten closer to correct.
        self.assertLess(abs(known_s - new_s), abs(known_s - wrong_s))
        self.assertLess(abs(known_t - new_t), abs(known_t - wrong_t))

    def test_convergence(self):
        import bezier

        nodes1 = np.asfortranarray([
            [0.0, 0.25, 0.5, 0.75, 1.0],
            [0.0, 1.0, -0.75, 1.0, 0.0],
        ])
        curve1 = bezier.Curve(nodes1, degree=4)
        # Vertical line forces a unique solution.
        nodes2 = np.asfortranarray([
            [0.5, 0.5],
            [0.0, 1.0],
        ])
        curve2 = bezier.Curve(nodes2, degree=1)

        num_guess = 4
        parameters = np.zeros((2, num_guess), order='F')
        # NOTE: This means our "first" guess is (s, t) = (0, 0).
        for guess in six.moves.xrange(1, num_guess):
            prev_s, prev_t = parameters[:, guess - 1]
            parameters[:, guess] = self._call_function_under_test(
                prev_s, nodes1, prev_t, nodes2)

        expected = np.asfortranarray([
            [0.0, 0.5, 0.5, 0.5],
            [0.0, 2.0, 0.21875, 0.21875],
        ])
        self.assertEqual(parameters, expected)
        # Make sure that we've actually converged.
        exact_s, exact_t = parameters[:, -1]
        self.assertEqual(curve1.evaluate(exact_s),
                         curve2.evaluate(exact_t))

    def test_singular_jacobian(self):
        nodes1 = np.asfortranarray([
            [0.5, 1.0, 1.5],
            [0.0, 1.0, 0.0],
        ])
        nodes2 = np.asfortranarray([
            [0.0, 1.0],
            [0.5, 0.5],
        ])

        with self.assertRaises(ValueError) as exc_info:
            self._call_function_under_test(
                0.5, nodes1, 0.5, nodes2)

        exc_args = exc_info.exception.args
        self.assertEqual(exc_args, ('Jacobian is singular.',))


class TestNewtonSimpleRoot(utils.NumPyTestCase):

    @staticmethod
    def _get_target_class():
        from bezier import _intersection_helpers

        return _intersection_helpers.NewtonSimpleRoot

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        nodes1 = unittest.mock.sentinel.nodes1
        first_deriv1 = unittest.mock.sentinel.first_deriv1
        nodes2 = unittest.mock.sentinel.nodes2
        first_deriv2 = unittest.mock.sentinel.first_deriv2
        evaluate_fn = self._make_one(
            nodes1, first_deriv1, nodes2, first_deriv2)
        self.assertIs(evaluate_fn.nodes1, nodes1)
        self.assertIs(evaluate_fn.first_deriv1, first_deriv1)
        self.assertIs(evaluate_fn.nodes2, nodes2)
        self.assertIs(evaluate_fn.first_deriv2, first_deriv2)

    def test___call__(self):
        # B1(s) = [s(s + 2) ]
        #         [4s(1 - s)]
        # B2(t) = [1 + 3t]
        #         [4 - 4t]
        # DF    = [2 + 2s, -3]
        #         [4 - 8s,  4]
        nodes1 = np.asfortranarray([
            [0.0, 1.0, 3.0],
            [0.0, 2.0, 0.0],
        ])
        first_deriv1 = np.asfortranarray([
            [2.0, 4.0],
            [4.0, -4.0],
        ])
        nodes2 = np.asfortranarray([
            [1.0, 4.0],
            [4.0, 0.0],
        ])
        first_deriv2 = np.asfortranarray([
            [3.0],
            [-4.0],
        ])

        evaluate_fn = self._make_one(
            nodes1, first_deriv1, nodes2, first_deriv2)
        jacobian, func_val = evaluate_fn(0.5, 0.25)

        expected_jacobian = np.asfortranarray([
            [3.0, -3.0],
            [0.0, 4.0],
        ])
        self.assertEqual(jacobian, expected_jacobian)
        expected_func_val = np.asfortranarray([
            [-0.5],
            [-2.0],
        ])
        self.assertEqual(func_val, expected_func_val)

    def test___call__exact_zero(self):
        # B1(s) = [2s(1 + s)]
        #         [6s(1 - s)]
        # B2(t) = [21t]
        #         [ 9t]
        # DF    = [2 +  4s, -21]
        #         [6 - 12s,  -9]
        nodes1 = np.asfortranarray([
            [0.0, 1.0, 4.0],
            [0.0, 3.0, 0.0],
        ])
        first_deriv1 = np.asfortranarray([
            [2.0, 6.0],
            [6.0, -6.0],
        ])
        nodes2 = np.asfortranarray([
            [0.0, 21.0],
            [0.0, 9.0],
        ])
        first_deriv2 = np.asfortranarray([
            [21.0],
            [9.0],
        ])

        evaluate_fn = self._make_one(
            nodes1, first_deriv1, nodes2, first_deriv2)
        jacobian, func_val = evaluate_fn(0.75, 0.125)

        self.assertIsNone(jacobian)
        expected_func_val = np.asfortranarray([
            [0.0],
            [0.0],
        ])
        self.assertEqual(func_val, expected_func_val)


class TestNewtonDoubleRoot(utils.NumPyTestCase):

    @staticmethod
    def _get_target_class():
        from bezier import _intersection_helpers

        return _intersection_helpers.NewtonDoubleRoot

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        nodes1 = unittest.mock.sentinel.nodes1
        first_deriv1 = unittest.mock.sentinel.first_deriv1
        second_deriv1 = unittest.mock.sentinel.second_deriv1
        nodes2 = unittest.mock.sentinel.nodes2
        first_deriv2 = unittest.mock.sentinel.first_deriv2
        second_deriv2 = unittest.mock.sentinel.second_deriv2

        evaluate_fn = self._make_one(
            nodes1, first_deriv1, second_deriv1,
            nodes2, first_deriv2, second_deriv2)
        self.assertIs(evaluate_fn.nodes1, nodes1)
        self.assertIs(evaluate_fn.first_deriv1, first_deriv1)
        self.assertIs(evaluate_fn.second_deriv1, second_deriv1)
        self.assertIs(evaluate_fn.nodes2, nodes2)
        self.assertIs(evaluate_fn.first_deriv2, first_deriv2)
        self.assertIs(evaluate_fn.second_deriv2, second_deriv2)

    def _make_default(self):
        # B1(s) = [4s(1 - s)]
        #         [       2s]
        # B2(t) = [2(2t^2 - 2t + 1)]
        #         [              2t]
        # B1'(s) x B2'(s) = -16(s + t - 1)
        # DG    = [4 - 8s, 4 - 8t]
        #         [     2,     -2]
        #         [   -16,    -16]
        # DG^T DG = [   4(16s^2 - 16s + 69), 4(16st - 8s - 8t + 67)]
        #           [4(16st - 8s - 8t + 67),    4(16t^2 - 16t + 69)]
        # DG^T G  = [4(8s^3 - 12s^2 + 8st^2 - 8st + 73s - 4t^2 + 67t - 66)]
        #           [4(8s^2t - 4s^2 - 8st + 67s + 8t^3 - 12t^2 + 73t - 66)]
        nodes1 = np.asfortranarray([
            [0.0, 2.0, 0.0],
            [0.0, 1.0, 2.0],
        ])
        first_deriv1 = np.asfortranarray([
            [4.0, -4.0],
            [2.0, 2.0],
        ])
        second_deriv1 = np.asfortranarray([
            [-8.0],
            [0.0],
        ])
        nodes2 = np.asfortranarray([
            [2.0, 0.0, 2.0],
            [0.0, 1.0, 2.0],
        ])
        first_deriv2 = np.asfortranarray([
            [-4.0, 4.0],
            [2.0, 2.0],
        ])
        second_deriv2 = np.asfortranarray([
            [8.0],
            [0.0],
        ])
        return self._make_one(
            nodes1, first_deriv1, second_deriv1,
            nodes2, first_deriv2, second_deriv2)

    def test___call__(self):
        evaluate_fn = self._make_default()
        jacobian, func_val = evaluate_fn(0.75, 0.25)

        expected_jacobian = np.asfortranarray([
            [264.0, 248.0],
            [248.0, 264.0],
        ])
        self.assertEqual(jacobian, expected_jacobian)
        expected_func_val = np.asfortranarray([
            [3.0],
            [-3.0],
        ])
        self.assertEqual(func_val, expected_func_val)

    def test___call__exact_zero(self):
        evaluate_fn = self._make_default()
        jacobian, func_val = evaluate_fn(0.5, 0.5)

        self.assertEqual(jacobian, None)
        expected_func_val = np.asfortranarray([
            [0.0],
            [0.0],
        ])
        self.assertEqual(func_val, expected_func_val)


class Test_newton_iterate(unittest.TestCase):

    HALF_EPS = 0.5**26

    @staticmethod
    def _call_function_under_test(evaluate_fn, s, t):
        from bezier import _intersection_helpers

        return _intersection_helpers.newton_iterate(evaluate_fn, s, t)

    @staticmethod
    def _simple_evaluate(quadratic1, quadratic2):
        from bezier import _intersection_helpers

        first_deriv1 = 2.0 * (quadratic1[:, 1:] - quadratic1[:, :-1])
        first_deriv2 = 2.0 * (quadratic2[:, 1:] - quadratic2[:, :-1])
        return _intersection_helpers.NewtonSimpleRoot(
            quadratic1, first_deriv1, quadratic2, first_deriv2)

    @staticmethod
    def _double_evaluate(quadratic1, quadratic2):
        from bezier import _intersection_helpers

        first_deriv1 = 2.0 * (quadratic1[:, 1:] - quadratic1[:, :-1])
        second_deriv1 = first_deriv1[:, 1:] - first_deriv1[:, :-1]
        first_deriv2 = 2.0 * (quadratic2[:, 1:] - quadratic2[:, :-1])
        second_deriv2 = first_deriv2[:, 1:] - first_deriv2[:, :-1]
        return _intersection_helpers.NewtonDoubleRoot(
            quadratic1, first_deriv1, second_deriv1,
            quadratic2, first_deriv2, second_deriv2)

    def test_rhs_exactly_zero(self):
        # B1([10922/32768, 10923/32768]) and B2([16383/16384, 1]) are
        # linearized and when the segments intersect they produce
        # t = 109217/109216 > 1.
        nodes1 = np.asfortranarray([
            [0.0, 4.5, 9.0],
            [0.0, 9.0, 0.0],
        ])
        s = 671023103.0 / 2013069312.0
        nodes2 = np.asfortranarray([
            [11.0, 7.0, 3.0],
            [8.0, 10.0, 4.0],
        ])
        t = 1789394945.0 / 1789394944.0

        evaluate_fn = self._simple_evaluate(nodes1, nodes2)
        converged, current_s, current_t = self._call_function_under_test(
            evaluate_fn, s, t)
        self.assertTrue(converged)
        self.assertEqual(3.0 * current_s, 1.0)
        self.assertEqual(current_t, 1.0)

    def test_singular_jacobian(self):
        # B1([5461/8192, 5462/8192]) and B2([2730/8192, 2731/8192]) are
        # linearized and the segments are parallel. The curves intersect
        # at the point B1(2/3) = [1/2, 1/2] = B2(1/3) and they have parallel
        # tangent vectors B1'(2/3) = [3/4, 0] = B2'(1/3).

        nodes1 = np.asfortranarray([
            [0.0, 0.375, 0.75],
            [0.0, 0.75, 0.375],
        ])
        s = 10923.0 / 16384.0
        nodes2 = np.asfortranarray([
            [0.25, 0.625, 1.0],
            [0.625, 0.25, 1.0],
        ])
        t = 5461.0 / 16384.0

        evaluate_fn = self._simple_evaluate(nodes1, nodes2)
        converged, current_s, current_t = self._call_function_under_test(
            evaluate_fn, s, t)
        self.assertFalse(converged)
        self.assertEqual(3.0 * current_s, 2.0 + 0.5**14)
        self.assertEqual(3.0 * current_t, 1 - 0.5**14)

    def _check_closer(
            self, s, current_s, expected_s, t, current_t, expected_t):
        # Make sure we are closer ...
        err_s = abs(expected_s - current_s)
        err_t = abs(expected_t - current_t)
        self.assertLess(err_s, abs(expected_s - s))
        self.assertLess(err_t, abs(expected_t - t))
        # ... but still not very close.
        self.assertGreater(err_s, self.HALF_EPS * expected_s)
        self.assertGreater(err_t, self.HALF_EPS * expected_t)

    def test_singular_jacobian_dg(self):
        # The curves are tangent and have the same curvature (i.e.
        # triple root).
        nodes1 = np.asfortranarray([
            [12.0, -4.0, -4.0],
            [4.0, -4.0, 4.0],
        ])
        s = float.fromhex('0x1.fffff4dad8308p-2')
        nodes2 = np.asfortranarray([
            [6.0, -2.0, -2.0],
            [1.0, -1.0, 1.0],
        ])
        t = float.fromhex('0x1.ffffe9b5a0f3ep-2')

        evaluate_fn = self._double_evaluate(nodes1, nodes2)
        converged, current_s, current_t = self._call_function_under_test(
            evaluate_fn, s, t)
        self.assertFalse(converged)
        self._check_closer(s, current_s, 0.5, t, current_t, 0.5)

    def test_convergence_linear(self):
        # B1([2730/8192, 2731/8192]) and B2([1/2, 2049/4096]) are
        # linearized and when the segments intersect they produce
        # t = -1/6 < 0.
        nodes1 = np.asfortranarray([
            [0.5, 1.25, 2.0],
            [0.125, -0.25, 0.5],
        ])
        s = 12287.0 / 36864.0
        nodes2 = np.asfortranarray([
            [0.5, 1.0, 1.5],
            [-0.125, 0.125, -0.125]
        ])
        t = 12287.0 / 24576.0

        # Due to linear convergence, this "bails out" after 5 iterations.
        evaluate_fn = self._simple_evaluate(nodes1, nodes2)
        converged, current_s, current_t = self._call_function_under_test(
            evaluate_fn, s, t)
        self.assertFalse(converged)
        self._check_closer(s, current_s, 1.0 / 3.0, t, current_t, 0.5)

    def test_convergence_linear_dg(self):
        # B1([16387/32768, 16388/32768]) and B2([8195/16384, 8196/16384]) are
        # linearized and when the segments intersect they produce
        # s = t = -9/7 < 0.
        nodes1 = np.asfortranarray([
            [12.0, -4.0, -4.0],
            [4.0, -4.0, 4.0],
        ])
        nodes2 = np.asfortranarray([
            [6.0, -2.0, -2.0],
            [1.0, -1.0, 1.0],
        ])
        # NOTE: These ``s-t`` values come after the simple root case exits
        #       due to linear convergence, having started from
        #       s = 28675 / 57344 and t = 14339 / 28672.
        s = float.fromhex('0x1.00006f0b2bb91p-1')
        t = float.fromhex('0x1.0000de165968ap-1')

        evaluate_fn = self._double_evaluate(nodes1, nodes2)
        converged, current_s, current_t = self._call_function_under_test(
            evaluate_fn, s, t)
        self.assertFalse(converged)
        self._check_closer(s, current_s, 0.5, t, current_t, 0.5)

    def test_below_error_ratio(self):
        # B1([12287/16384, 3/4]) and B2([2457/8192, 2458/8192]) are linearized
        # and when the segments intersect they produce
        # s = 33555797/33551701 > 1.
        nodes1 = np.asfortranarray([
            [1.0, -1.0, 1.0],
            [0.0, 0.25, 0.5],
        ])
        s = 25163776.0 / 33551701.0
        nodes2 = np.asfortranarray([
            [-0.125, 0.5, 1.125],
            [-0.28125, 1.28125, -0.28125],
        ])
        t = 41228331827.0 / 137427767296.0

        evaluate_fn = self._simple_evaluate(nodes1, nodes2)
        converged, current_s, current_t = self._call_function_under_test(
            evaluate_fn, s, t)

        self.assertTrue(converged)
        self.assertEqual(0.75, current_s)
        utils.almost(self, 3.0 / 10.0, current_t, 1)

    def test_below_error_ratio_dg(self):
        # B1([2730/8192, 2731/8192]) and B2([2047/4096, 1/2]) are
        # linearized and when the segments intersect they produce
        # t = 11/10 > 1.
        nodes1 = np.asfortranarray([
            [0.5, 1.25, 2.0],
            [0.125, -0.25, 0.5],
        ])
        nodes2 = np.asfortranarray([
            [0.5, 1.0, 1.5],
            [-0.125, 0.125, -0.125]
        ])
        # NOTE: These ``s-t`` values come after the simple root case exits
        #       due to linear convergence, having started from
        #       s = 6827 / 20480 and t = 20481 / 40960 and updating 4 times.
        s = 109227.0 / 327680.0
        t = 327681.0 / 655360.0

        evaluate_fn = self._double_evaluate(nodes1, nodes2)
        converged, current_s, current_t = self._call_function_under_test(
            evaluate_fn, s, t)

        self.assertTrue(converged)
        utils.almost(self, 1.0 / 3.0, current_s, 1)
        self.assertEqual(0.5, current_t)

    def test_all_iterations(self):
        # We start very far away from the root s = t = 0.5 by using
        # s, t >> 1. (There is also a root at s = t = 0.0.)
        nodes1 = np.asfortranarray([
            [0.0, 1.0, 2.0],
            [0.0, 2.0, 0.0],
        ])
        s = 64.0
        nodes2 = np.asfortranarray([
            [0.0, 2.0, 0.0],
            [0.0, 1.0, 2.0],
        ])
        t = 64.0

        evaluate_fn = self._simple_evaluate(nodes1, nodes2)
        # We "fake" MAX_NEWTON_ITERATIONS=3 because when we are sufficiently
        # far from a root, convergence **appears** linear. This is because
        # ``pn`` is moving significantly, so the change in ``||pn||`` tracks
        # the change in ``||p{n+1} - pn||``.
        patch = unittest.mock.patch(
            'bezier._intersection_helpers.MAX_NEWTON_ITERATIONS', new=3)
        with patch:
            converged, current_s, current_t = self._call_function_under_test(
                evaluate_fn, s, t)

        self.assertFalse(converged)
        self.assertEqual(8.221323371115176, current_s)
        self.assertEqual(8.221323371115176, current_t)


class Test_full_newton_nonzero(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(s, nodes1, t, nodes2):
        from bezier import _intersection_helpers

        return _intersection_helpers.full_newton_nonzero(s, nodes1, t, nodes2)

    def test_simple_root(self):
        # B1([4095/8192, 1/2]) and B2([1365/8192, 1366/8192]) are linearized
        # and when the segments intersect they produce s = 24580/24579 > 1.
        nodes1 = np.asfortranarray([
            [0.0, 0.375, 0.75],
            [0.0, 0.75, 0.375],
        ])
        s = 100675585.0 / 201351168.0
        nodes2 = np.asfortranarray([
            [0.25, 0.625, 1.0],
            [0.5625, 0.1875, 0.9375],
        ])
        t = 33558529.0 / 201351168.0

        computed_s, computed_t = self._call_function_under_test(
            s, nodes1, t, nodes2)
        utils.almost(self, 0.5, computed_s, 1)
        utils.almost(self, 1.0 / 6.0, computed_t, 3)

    def test_double_root(self):
        # B1([5461/8192, 5462/8192]) and B2([2730/8192, 2731/8192]) are
        # linearized and the segments are parallel. The curves intersect
        # at the point B1(2/3) = [1/2, 1/2] = B2(1/3) and they have parallel
        # tangent vectors B1'(2/3) = [3/4, 0] = B2'(1/3).

        nodes1 = np.asfortranarray([
            [0.0, 0.375, 0.75],
            [0.0, 0.75, 0.375],
        ])
        s = 10923.0 / 16384.0
        nodes2 = np.asfortranarray([
            [0.25, 0.625, 1.0],
            [0.625, 0.25, 1.0],
        ])
        t = 5461.0 / 16384.0

        computed_s, computed_t = self._call_function_under_test(
            s, nodes1, t, nodes2)
        utils.almost(self, 2.0 / 3.0, computed_s, 1)
        self.assertEqual(1.0 / 3.0, computed_t)

    def test_triple_root(self):
        from bezier import _intersection_helpers

        # B1([16382/32768, 16383/32768]) and B2([8190/16384, 8191/16384]) are
        # linearized and when the segments intersect they produce
        # s = t = 4/3 > 1.
        nodes1 = np.asfortranarray([
            [12.0, -4.0, -4.0],
            [4.0, -4.0, 4.0],
        ])
        s = 24575.0 / 49152.0
        nodes2 = np.asfortranarray([
            [6.0, -2.0, -2.0],
            [1.0, -1.0, 1.0],
        ])
        t = 12287.0 / 24576.0

        with self.assertRaises(NotImplementedError) as exc_info:
            self._call_function_under_test(s, nodes1, t, nodes2)

        expected = (_intersection_helpers.NEWTON_NO_CONVERGE,)
        self.assertEqual(exc_info.exception.args, expected)

    def test_line_and_curve(self):
        # B1([5461/16384, 5462/16384]) and B2([0, 1]) are linearized
        # and when the segments intersect they produce s = -1/3 < 0.
        nodes1 = np.asfortranarray([
            [0.0, 1.5, 3.0],
            [2.25, -2.25, 2.25],
        ])
        s = 8191.0 / 24576.0
        nodes2 = np.asfortranarray([
            [-0.5, 4.0],
            [1.75, -2.75],
        ])
        t = 12287.0 / 36864.0

        computed_s, computed_t = self._call_function_under_test(
            s, nodes1, t, nodes2)
        utils.almost(self, 1.0 / 3.0, computed_s, 1)
        self.assertEqual(1.0 / 3.0, computed_t)


class Test_full_newton(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(s, nodes1, t, nodes2):
        from bezier import _intersection_helpers

        return _intersection_helpers.full_newton(s, nodes1, t, nodes2)

    def test_both_near_zero(self):
        # B1([0, 1/8192]) and B2([0, 1/8192]) are linearized and the
        # segments are parallel, and the root is a double root.
        nodes1 = np.asfortranarray([
            [1.0, 1.0, 0.0],
            [0.0, 1.0, 1.0],
        ])
        s = 1.0 / 16384.0
        nodes2 = np.asfortranarray([
            [1.0, 1.0, -0.5],
            [0.0, 1.5, 1.5],
        ])
        t = 1.0 / 16384.0

        computed_s, computed_t = self._call_function_under_test(
            s, nodes1, t, nodes2)
        self.assertEqual(0.0, computed_s)
        self.assertEqual(0.0, computed_t)

    def _one_parameter_near_zero(self, swap=False):
        # B1([0, 1/8192]) and B2([1/4, 2049/8192]) are linearized and the
        # segments are parallel, and the root is a double root.
        nodes1 = np.asfortranarray([
            [1.0, 1.0, -0.5],
            [0.0, 1.5, 1.5],
        ])
        s = 1.0 / 16384.0
        nodes2 = np.asfortranarray([
            [0.9375, 1.1875, 0.4375],
            [-0.5625, 0.6875, 0.9375],
        ])
        t = 4097.0 / 16384.0

        if swap:
            nodes1, nodes2 = nodes2, nodes1
            s, t = t, s

        computed_s, computed_t = self._call_function_under_test(
            s, nodes1, t, nodes2)

        if swap:
            self.assertEqual(0.25, computed_s)
            self.assertEqual(0.0, computed_t)
        else:
            self.assertEqual(0.0, computed_s)
            self.assertEqual(0.25, computed_t)

    def test_s_near_zero(self):
        self._one_parameter_near_zero()

    def test_t_near_zero(self):
        self._one_parameter_near_zero(swap=True)

    def test_both_nonzero(self):
        # B1([6826/8192, 6827/8192]) and B2([1/2, 4097/8192]) are linearized
        # and when the segments intersect they produce t = -1/24579 < 0.
        # The root is a simple root.
        nodes1 = np.asfortranarray([
            [0.0, 0.375, 0.75],
            [0.0, 0.75, 0.375],
        ])
        s = 167792639.0 / 201351168.0
        nodes2 = np.asfortranarray([
            [0.25, 0.625, 1.0],
            [0.5625, 0.1875, 0.9375],
        ])
        t = 100675583.0 / 201351168.0

        computed_s, computed_t = self._call_function_under_test(
            s, nodes1, t, nodes2)
        utils.almost(self, 5.0 / 6.0, computed_s, 2)
        utils.almost(self, 0.5, computed_t, 1)


@utils.needs_speedup
class Test_speedup_newton_refine(Test__newton_refine):

    @staticmethod
    def _call_function_under_test(s, nodes1, t, nodes2):
        from bezier import _speedup

        return _speedup.newton_refine_curve_intersect(s, nodes1, t, nodes2)


class TestIntersection(unittest.TestCase):

    @staticmethod
    def _get_target_class():
        from bezier import _intersection_helpers

        return _intersection_helpers.Intersection

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def _constructor_helper(self, **kwargs):
        index_first = 2
        s_val = 0.25
        index_second = 1
        t_val = 0.75
        intersection = self._make_one(
            index_first, s_val, index_second, t_val, **kwargs)

        self.assertEqual(intersection.index_first, index_first)
        self.assertEqual(intersection.s, s_val)
        self.assertEqual(intersection.index_second, index_second)
        self.assertEqual(intersection.t, t_val)
        return intersection

    def test_constructor(self):
        intersection = self._constructor_helper()
        self.assertIsNone(intersection.interior_curve)

    def test_constructor_with_interior_curve(self):
        intersection = self._constructor_helper(
            interior_curve=unittest.mock.sentinel.interior_curve)
        self.assertIs(
            intersection.interior_curve,
            unittest.mock.sentinel.interior_curve)

    def test___dict___property(self):
        intersection = self._constructor_helper(
            interior_curve=unittest.mock.sentinel.interior_curve,
        )
        props_dict = intersection.__dict__
        expected = {
            'index_first': 2,
            's': 0.25,
            'index_second': 1,
            't': 0.75,
            'interior_curve': unittest.mock.sentinel.interior_curve,
        }
        self.assertEqual(props_dict, expected)
        # Check that modifying ``props_dict`` won't modify ``curve``.
        props_dict['s'] = 0.5
        self.assertNotEqual(intersection.s, props_dict['s'])
