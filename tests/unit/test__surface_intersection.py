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


UNIT_TRIANGLE = np.asfortranarray([
    [0.0, 0.0],
    [1.0, 0.0],
    [0.0, 1.0],
])
SPACING = np.spacing  # pylint: disable=no-member


class Test_newton_refine_solve(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(jac_both, x_val, surf_x, y_val, surf_y):
        from bezier import _surface_intersection

        return _surface_intersection.newton_refine_solve(
            jac_both, x_val, surf_x, y_val, surf_y)

    def test_it(self):
        jac_both = np.asfortranarray([[1.0, 1.0, -2.0, 2.0]])
        delta_s, delta_t = self._call_function_under_test(
            jac_both, 0.5, 0.25, 0.75, 1.25)
        self.assertEqual(delta_s, -0.125)
        self.assertEqual(delta_t, -0.1875)


class Test__newton_refine(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, x_val, y_val, s, t):
        from bezier import _surface_intersection

        return _surface_intersection._newton_refine(
            nodes, degree, x_val, y_val, s, t)

    def test_improvement(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, -0.25],
            [1.0, 0.0],
            [0.0, 0.5],
            [0.5, 0.5],
            [-0.25, 0.875],
        ])
        # This surface is given by
        #     [(4 s - t^2) / 4, (4 s^2 + 4 s t - t^2 - 4 s + 8 t) / 8]
        s = 0.25
        t = 0.5
        # At our points, the Jacobian is
        #     [1, -1/4]
        #     [0,  1  ]
        # hence there will be no round-off when applying the inverse.
        # (x_val, y_val), = surface.evaluate_cartesian(0.5, 0.25)
        x_val = 0.484375
        y_val = 0.1796875
        new_s, new_t = self._call_function_under_test(
            nodes, 2, x_val, y_val, s, t)
        self.assertEqual(new_s, 247.0 / 512.0)
        self.assertEqual(new_t, 31.0 / 128.0)

    def test_at_solution(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.0],
            [0.0, 0.5],
            [0.5, 0.5],
            [0.0, 1.0],
        ])
        # This surface is given by [s, t].
        s = 0.375
        t = 0.75
        # Since x(s) = s and y(t) = t, we simply use the same x/y and s/t.
        x_val = s
        y_val = t
        new_s, new_t = self._call_function_under_test(
            nodes, 2, x_val, y_val, s, t)
        self.assertEqual(new_s, s)
        self.assertEqual(new_t, t)


@utils.needs_surface_intersection_speedup
class Test_speedup_newton_refine(Test__newton_refine):

    @staticmethod
    def _call_function_under_test(nodes, degree, x_val, y_val, s, t):
        from bezier import _surface_intersection_speedup

        return _surface_intersection_speedup.newton_refine(
            nodes, degree, x_val, y_val, s, t)


class Test_update_locate_candidates(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(
            candidate, next_candidates, x_val, y_val, degree):
        from bezier import _surface_intersection

        return _surface_intersection.update_locate_candidates(
            candidate, next_candidates, x_val, y_val, degree)

    @unittest.mock.patch(
        'bezier._surface_helpers.subdivide_nodes',
        return_value=(
            unittest.mock.sentinel.nodes_a,
            unittest.mock.sentinel.nodes_b,
            unittest.mock.sentinel.nodes_c,
            unittest.mock.sentinel.nodes_d,
        ))
    def test_contained(self, subdivide_nodes):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.25],
            [1.0, 0.0],
            [0.0, 0.5],
            [0.75, 0.75],
            [0.0, 1.0],
        ])
        candidate = (
            1.25,
            1.25,
            -0.25,
            nodes,
        )
        next_candidates = []

        ret_val = self._call_function_under_test(
            candidate, next_candidates, 0.5625, 0.375, 2)
        self.assertIsNone(ret_val)

        expected = [
            (1.375, 1.375, -0.125, unittest.mock.sentinel.nodes_a),
            (1.25, 1.25, 0.125, unittest.mock.sentinel.nodes_b),
            (1.0, 1.375, -0.125, unittest.mock.sentinel.nodes_c),
            (1.375, 1.0, -0.125, unittest.mock.sentinel.nodes_d),
        ]
        self.assertEqual(next_candidates, expected)
        subdivide_nodes.assert_called_once_with(nodes, 2)

    def test_not_contained(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [2.0, 3.0],
            [-1.0, 2.0],
        ])
        candidate = (
            2.0,
            0.5,
            0.5,
            nodes,
        )
        next_candidates = []

        ret_val = self._call_function_under_test(
            candidate, next_candidates, 9.0, 9.0, 1)
        self.assertIsNone(ret_val)

        self.assertEqual(next_candidates, [])


class Test_mean_centroid(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(candidates):
        from bezier import _surface_intersection

        return _surface_intersection.mean_centroid(candidates)

    def test_it(self):
        candidates = (
            (1.0, 1.0, None, None),
            (2.25, 1.5, None, None),
            (1.25, 4.25, None, None),
        )

        centroid_x, centroid_y = self._call_function_under_test(candidates)
        self.assertEqual(centroid_x, 0.5)
        self.assertEqual(centroid_y, 0.75)


class Test__locate_point(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes, degree, x_val, y_val):
        from bezier import _surface_intersection

        return _surface_intersection._locate_point(nodes, degree, x_val, y_val)

    def test_it(self):
        nodes = UNIT_TRIANGLE.copy(order='F')
        degree = 1
        x_val = 0.25
        y_val = 0.625
        s, t = self._call_function_under_test(nodes, degree, x_val, y_val)
        self.assertEqual(s, x_val)
        self.assertEqual(t, y_val)

    def test_extra_newton_step(self):
        # x(s, t) = -2 (s + 2 t) (t - 1)
        # y(s, t) =  2 (s + 1) t
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [2.0, 2.0],
            [0.0, 2.0],
        ])
        degree = 2
        x_val = 0.59375
        y_val = 0.25
        s, t = self._call_function_under_test(nodes, degree, x_val, y_val)
        # NOTE: We can use the resultant to find that the **true** answers
        #       are roots of the following polynomials.
        #           64 s^3 + 101 s^2 + 34 s - 5 = 0
        #           128 t^3 - 192 t^2 + 91 t - 8 = 0
        #       Using extended precision, we can find these values to more
        #       digits than what is supported by IEEE-754.
        expected_s = 0.109190958136897160638
        self.assertAlmostEqual(s, expected_s, delta=6 * SPACING(expected_s))
        expected_t = 0.11269475204698919699
        self.assertAlmostEqual(t, expected_t, delta=SPACING(expected_t))

    def test_no_match(self):
        nodes = UNIT_TRIANGLE.copy(order='F')
        degree = 1
        x_val = -0.125
        y_val = 0.25
        self.assertIsNone(
            self._call_function_under_test(nodes, degree, x_val, y_val))


@utils.needs_surface_intersection_speedup
class Test_speedup_locate_point(Test__locate_point):

    @staticmethod
    def _call_function_under_test(nodes, degree, x_val, y_val):
        from bezier import _surface_intersection_speedup

        return _surface_intersection_speedup.locate_point(
            nodes, degree, x_val, y_val)
