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
            [0.0, 0.0],
            [1.0, 1.0],
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
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ])
        nodes2 = np.asfortranarray([
            [1.0, 0.75],
            [0.5, -0.25],
            [0.0, 0.75],
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
            [0.0, 0.0],
            [0.25, 1.0],
            [0.5, -0.75],
            [0.75, 1.0],
            [1.0, 0.0],
        ])
        curve1 = bezier.Curve(nodes1, degree=4)
        # Vertical line forces a unique solution.
        nodes2 = np.asfortranarray([
            [0.5, 0.0],
            [0.5, 1.0],
        ])
        curve2 = bezier.Curve(nodes2, degree=1)

        num_guess = 4
        parameters = np.zeros((num_guess, 2), order='F')
        # NOTE: This means our "first" guess is (s, t) = (0, 0).
        for guess in six.moves.xrange(1, num_guess):
            prev_s, prev_t = parameters[guess - 1, :]
            parameters[guess, :] = self._call_function_under_test(
                prev_s, nodes1, prev_t, nodes2)

        expected = np.asfortranarray([
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


@utils.needs_curve_intersection_speedup
class Test_speedup_newton_refine(Test__newton_refine):

    @staticmethod
    def _call_function_under_test(s, nodes1, t, nodes2):
        from bezier import _curve_intersection_speedup

        return _curve_intersection_speedup.newton_refine(
            s, nodes1, t, nodes2)


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
        expected['s'] = 0.5
        self.assertNotEqual(intersection.s, expected['s'])
