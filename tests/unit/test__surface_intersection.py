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

import threading
import unittest
import unittest.mock

import numpy as np

from tests import utils as base_utils
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


@utils.needs_speedup
class Test_speedup_newton_refine(Test__newton_refine):

    @staticmethod
    def _call_function_under_test(nodes, degree, x_val, y_val, s, t):
        from bezier import _speedup

        return _speedup.newton_refine_surface(
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


@utils.needs_speedup
class Test_speedup_locate_point(Test__locate_point):

    @staticmethod
    def _call_function_under_test(nodes, degree, x_val, y_val):
        from bezier import _speedup

        return _speedup.locate_point_surface(
            nodes, degree, x_val, y_val)


class Test_same_intersection(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(intersection1, intersection2, **kwargs):
        from bezier import _surface_intersection

        return _surface_intersection.same_intersection(
            intersection1, intersection2, **kwargs)

    def test_same(self):
        intersection = make_intersect(10, 0.5, 99, 0.75)
        result = self._call_function_under_test(intersection, intersection)
        self.assertTrue(result)

    def test_almost_same(self):
        intersection1 = make_intersect(10, 0.5, 99, 0.75)
        intersection2 = make_intersect(10, 0.5, 99, 0.875)
        result = self._call_function_under_test(intersection1, intersection2)
        self.assertFalse(result)
        result = self._call_function_under_test(
            intersection1, intersection2, wiggle=0.5)
        self.assertTrue(result)

    def test_different_edge(self):
        intersection1 = make_intersect(10, 0.5, 99, 0.5)
        intersection2 = make_intersect(10, 0.5, 98, 0.5)
        intersection3 = make_intersect(11, 0.5, 99, 0.5)
        self.assertFalse(
            self._call_function_under_test(intersection1, intersection2))
        self.assertFalse(
            self._call_function_under_test(intersection1, intersection3))

    def test_different_param(self):
        intersection1 = make_intersect(1, 0.5, 9, 0.5)
        intersection2 = make_intersect(1, 0.75, 9, 0.5)
        intersection3 = make_intersect(1, 0.5, 9, 0.75)
        self.assertFalse(
            self._call_function_under_test(intersection1, intersection2))
        self.assertFalse(
            self._call_function_under_test(intersection1, intersection3))


class Test_verify_duplicates(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(duplicates, uniques):
        from bezier import _surface_intersection

        return _surface_intersection.verify_duplicates(duplicates, uniques)

    def test_empty(self):
        self.assertIsNone(self._call_function_under_test([], []))

    def test_success(self):
        uniq = make_intersect(1, 0.0, 2, 0.25)
        self.assertIsNone(
            self._call_function_under_test([uniq], [uniq]))

    def test_success_triple(self):
        uniq = make_intersect(1, 0.0, 2, 0.0)
        self.assertIsNone(
            self._call_function_under_test([uniq, uniq, uniq], [uniq]))

    def test_failed_uniqueness(self):
        uniq = make_intersect(1, 0.375, 2, 0.75)
        with self.assertRaises(ValueError):
            self._call_function_under_test([], [uniq, uniq])

    def test_bad_duplicate(self):
        dupe = make_intersect(1, 0.75, 2, 0.25)
        uniq = make_intersect(1, 0.25, 2, 0.75)
        with self.assertRaises(ValueError):
            self._call_function_under_test([dupe], [uniq])

    def test_bad_single_corner(self):
        uniq = make_intersect(1, 0.125, 2, 0.125)
        with self.assertRaises(ValueError):
            self._call_function_under_test([uniq], [uniq])

    def test_bad_double_corner(self):
        uniq = make_intersect(1, 0.0, 2, 1.0)
        with self.assertRaises(ValueError):
            self._call_function_under_test([uniq, uniq, uniq], [uniq])

    def test_bad_count(self):
        uniq = make_intersect(1, 0.375, 2, 0.75)
        with self.assertRaises(ValueError):
            self._call_function_under_test([uniq, uniq], [uniq])


class Test_surface_intersections(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(edge_nodes1, edge_nodes2, all_intersections):
        from bezier import _surface_intersection

        return _surface_intersection.surface_intersections(
            edge_nodes1, edge_nodes2, all_intersections)

    def _check_intersection(
            self, intersection, index_first, s_val, index_second,
            t_val, interior_curve):
        from bezier import _intersection_helpers

        self.assertIsInstance(
            intersection, _intersection_helpers.Intersection)
        self.assertEqual(intersection.index_first, index_first)
        self.assertEqual(intersection.s, s_val)
        self.assertEqual(intersection.index_second, index_second)
        self.assertEqual(intersection.t, t_val)
        self.assertEqual(intersection.interior_curve, interior_curve)

    def test_with_unused(self):
        import bezier
        from bezier import _geometric_intersection

        surface1 = bezier.Surface.from_nodes(np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ]))
        surface2 = bezier.Surface.from_nodes(np.asfortranarray([
            [0.0, 0.0],
            [-2.0, 3.0],
            [1.0, -3.0],
        ]))

        edge_nodes1 = tuple(edge._nodes for edge in surface1.edges)
        edge_nodes2 = tuple(edge._nodes for edge in surface2.edges)
        all_intersections = _geometric_intersection.all_intersections
        result = self._call_function_under_test(
            edge_nodes1, edge_nodes2, all_intersections)
        intersections, duplicates, unused, all_types = result

        self.assertEqual(intersections, [])

        self.assertEqual(len(duplicates), 3)
        self._check_intersection(
            duplicates[0], 0, 0.0, 0, 0.0, None)
        self._check_intersection(
            duplicates[1], 0, 0.0, 0, 0.0, None)
        self._check_intersection(
            duplicates[2], 0, 0.0, 0, 0.0, None)

        self.assertEqual(len(unused), 1)
        enum_val = get_enum('IGNORED_CORNER')
        self._check_intersection(
            unused[0], 0, 0.0, 0, 0.0, enum_val)
        self.assertEqual(all_types, set([enum_val]))


class Test_generic_intersect(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(
            nodes1, degree1, nodes2, degree2, verify, all_intersections):
        from bezier import _surface_intersection

        return _surface_intersection.generic_intersect(
            nodes1, degree1, nodes2, degree2, verify, all_intersections)

    def test_disjoint_bbox(self):
        from bezier import _geometric_intersection

        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        nodes2 = np.asfortranarray([
            [10.0, 10.0],
            [11.0, 10.0],
            [10.0, 11.0],
        ])
        all_intersections = _geometric_intersection.all_intersections

        result = self._call_function_under_test(
            nodes1, 1, nodes2, 1, True, all_intersections)
        edge_infos, contained, all_edge_nodes = result

        self.assertEqual(edge_infos, [])
        self.assertIsNone(contained)
        self.assertEqual(all_edge_nodes, ())


class Test__geometric_intersect(utils.NumPyTestCase):

    NODES1 = np.asfortranarray([
        [-8.0, 0.0],
        [8.0, 0.0],
        [0.0, 8.0],
    ])
    NODES2 = np.asfortranarray([
        [4.0, 3.0],
        [0.0, -5.0],
        [-4.0, 3.0],
        [2.0, -3.0],
        [-2.0, -3.0],
        [0.0, -9.0],
    ])
    BAD_BOUNDARY_ARGS = ('Non-unique intersection',)
    BAD_BOUNDARY_TYPE = ValueError
    BAD_BOUNDARY_INCREASE_ULPS = True

    @staticmethod
    def _call_function_under_test(nodes1, degree1, nodes2, degree2, **kwargs):
        from bezier import _surface_intersection

        return _surface_intersection._geometric_intersect(
            nodes1, degree1, nodes2, degree2, **kwargs)

    @staticmethod
    def parallel_err():
        from bezier import _geometric_intersection

        return _geometric_intersection._SEGMENTS_PARALLEL

    def test_disjoint_bbox(self):
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        nodes2 = np.asfortranarray([
            [10.0, 10.0],
            [11.0, 10.0],
            [10.0, 11.0],
        ])

        result = self._call_function_under_test(
            nodes1, 1, nodes2, 1, verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self.assertEqual(curved_polygons, [])
        self.assertIsNone(contained)
        self.assertEqual(all_edge_nodes, ())

    def test_parallel(self):
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [5.0, 0.0],
            [0.0, 5.0],
        ])
        nodes2 = np.asfortranarray([
            [2.0, 0.0],
            [3.0, 0.0],
            [2.0, 1.0],
        ])

        with self.assertRaises(NotImplementedError) as exc_info:
            self._call_function_under_test(nodes1, 1, nodes2, 1, verify=True)

        exc_args = exc_info.exception.args
        self.assertEqual(exc_args, (self.parallel_err(),))

    def test_contained(self):
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [5.0, 0.0],
            [0.0, 5.0],
        ])
        nodes2 = np.asfortranarray([
            [2.0, 1.0],
            [3.0, 1.0],
            [2.0, 2.0],
        ])

        result = self._call_function_under_test(
            nodes1, 1, nodes2, 1, verify=False)
        curved_polygons, contained, all_edge_nodes = result
        self.assertIsNone(curved_polygons)
        self.assertIsNotNone(contained)
        self.assertFalse(contained)
        self.assertEqual(all_edge_nodes, ())  # NOT EMPTY GEOM.

        result = self._call_function_under_test(
            nodes2, 1, nodes1, 1, verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self.assertIsNone(curved_polygons)
        self.assertTrue(contained)
        self.assertEqual(all_edge_nodes, ())

    def test_disjoint_with_intersecting_bbox(self):
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [5.0, 0.0],
            [0.0, 5.0],
        ])
        nodes2 = np.asfortranarray([
            [4.0, 2.0],
            [4.0, 3.0],
            [3.0, 3.0],
        ])

        result = self._call_function_under_test(
            nodes1, 1, nodes2, 1, verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self.assertEqual(curved_polygons, [])
        self.assertIsNone(contained)
        self.assertEqual(all_edge_nodes, ())

    def _check_linear_intersection1(self, value):
        self.assertEqual(value, 0.125)

    def _check_linear_intersection(self, curved_polygons):
        self.assertEqual(len(curved_polygons), 1)
        edge_info = curved_polygons[0]
        self.assertIsInstance(edge_info, tuple)
        self.assertEqual(len(edge_info), 6)
        self.assertEqual(edge_info[0], (4, 0.5, 0.625))
        self.assertEqual(edge_info[2], (5, 0.375, 0.875))
        self.assertEqual(edge_info[3], (1, 0.5, 0.625))
        self.assertEqual(edge_info[4], (3, 0.125, 0.5))
        self.assertEqual(edge_info[5], (2, 0.375, 0.875))

        # Add special handling for edge 1 (algebraic messes this up).
        self.assertIsInstance(edge_info[1], tuple)
        self.assertEqual(len(edge_info[1]), 3)
        self.assertEqual(edge_info[1][0], 0)
        self.assertEqual(edge_info[1][2], 0.5)
        self._check_linear_intersection1(edge_info[1][1])

    def test_linear_intersection(self):
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [8.0, 0.0],
            [0.0, 8.0],
        ])
        nodes2 = np.asfortranarray([
            [4.0, 5.0],
            [-4.0, 5.0],
            [4.0, -3.0],
        ])

        result = self._call_function_under_test(
            nodes1, 1, nodes2, 1, verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self._check_linear_intersection(curved_polygons)
        self.assertIsNone(contained)
        check_edges(self, nodes1, 1, nodes2, 1, all_edge_nodes)

    def test_opposed_tangencies(self):
        nodes1 = np.asfortranarray([
            [4.0, 0.0],
            [2.0, 4.0],
            [0.0, 0.0],
            [3.0, -2.0],
            [1.0, -2.0],
            [2.0, -4.0],
        ])
        nodes2 = np.asfortranarray([
            [0.0, 4.0],
            [2.0, 0.0],
            [4.0, 4.0],
            [1.0, 6.0],
            [3.0, 6.0],
            [2.0, 8.0],
        ])

        result = self._call_function_under_test(
            nodes1, 2, nodes2, 2, verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self.assertEqual(curved_polygons, [])
        self.assertIsNone(contained)
        self.assertEqual(all_edge_nodes, ())

    def test_ignored_corner(self):
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [-2.0, 3.0],
            [1.0, -3.0],
        ])

        result = self._call_function_under_test(
            nodes1, 1, nodes2, 1, verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self.assertEqual(curved_polygons, [])
        self.assertIsNone(contained)
        self.assertEqual(all_edge_nodes, ())

    def test_tangent_contained(self):
        nodes1 = np.asfortranarray([
            [2.0, 1.25],
            [4.0, 0.75],
            [6.0, 1.25],
            [3.0, 3.0],
            [5.0, 3.0],
            [4.0, 5.0],
        ])
        nodes2 = np.asfortranarray([
            [0.0, 0.0],
            [4.0, 2.0],
            [8.0, 0.0],
            [1.5, 7.5],
            [11.0, 6.0],
            [8.0, 8.0],
        ])

        result = self._call_function_under_test(
            nodes1, 2, nodes2, 2, verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self.assertIsNone(curved_polygons)
        self.assertTrue(contained)
        self.assertEqual(all_edge_nodes, ())

        result = self._call_function_under_test(
            nodes2, 2, nodes1, 2, verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self.assertIsNone(curved_polygons)
        self.assertIsNotNone(contained)
        self.assertFalse(contained)
        self.assertEqual(all_edge_nodes, ())

    def _two_curved_polygons(self, **kwargs):
        expected = [
            (
                (5, 0.75, 1.0),
                (3, 0.0, 0.25),
                (0, 0.625, 0.6875),
            ), (
                (0, 0.3125, 0.375),
                (3, 0.75, 1.0),
                (4, 0.0, 0.25),
            ),
        ]
        return (
            self._call_function_under_test(
                self.NODES1, 1, self.NODES2, 2, **kwargs),
            expected,
        )

    def test_two_curved_polygons(self):
        result, expected = self._two_curved_polygons(verify=True)
        curved_polygons, contained, all_edge_nodes = result
        self.assertEqual(curved_polygons, expected)
        self.assertIsNone(contained)
        check_edges(self, self.NODES1, 1, self.NODES2, 2, all_edge_nodes)

    def _almost(self, actual, expected, num_ulps):
        delta = num_ulps * np.spacing(expected)
        self.assertAlmostEqual(actual, expected, delta=delta)

    def test_bad_boundary(self):
        from bezier import _geometric_intersection

        nodes1 = np.asfortranarray([
            [-0.519247936646441, 0.008196262233806585],
            [-0.5206153113582636, 0.01587839530234961],
            [-0.5220361984817866, 0.023564120023660037],
            [-0.5234919947680754, 0.03126920999203689],
            [-0.5306013153328831, 0.010994185722894958],
            [-0.5320185860555878, 0.018720567953255777],
            [-0.533471521483344, 0.026443532602826243],
            [-0.5418564772120571, 0.013842126833196683],
            [-0.5433055706504355, 0.02160135650021924],
            [-0.5530139922771271, 0.016726767154940626],
        ])
        nodes2 = np.asfortranarray([
            [-0.5492475273303934, 0.004806678750684627],
            [-0.5507539103531026, 0.011215423321262663],
            [-0.552283211157318, 0.017577411042302475],
            [-0.553868947686436, 0.023906392050982415],
            [-0.5543287770635862, 0.0032189095236835417],
            [-0.5558905319225977, 0.00960140358887212],
            [-0.5574932750785362, 0.015938552164569394],
            [-0.5596130517429451, 0.0014977481523963628],
            [-0.5612419767391262, 0.007849282849502192],
            [-0.5650788564504334, -0.0003406439314109452],
        ])
        with self.assertRaises(self.BAD_BOUNDARY_TYPE) as exc_info:
            self._call_function_under_test(
                nodes1, 3, nodes2, 3, verify=True)

        self.assertEqual(exc_info.exception.args, self.BAD_BOUNDARY_ARGS)

        if not self.BAD_BOUNDARY_INCREASE_ULPS:
            return

        # Increase ``SIMILAR_ULPS`` so that the intersection succeeds.
        similar_ulps = _geometric_intersection.get_similar_ulps()
        if (base_utils.IS_MAC_OS_X and
                not base_utils.IS_64_BIT):  # pragma: NO COVER
            _geometric_intersection.set_similar_ulps(601)
        else:
            _geometric_intersection.set_similar_ulps(418)

        try:
            result = self._call_function_under_test(
                nodes1, 3, nodes2, 3, verify=True)
        finally:
            # Restore the original value.
            _geometric_intersection.set_similar_ulps(similar_ulps)

        curved_polygons, contained, all_edge_nodes = result
        self.assertIsNone(contained)
        check_edges(self, nodes1, 3, nodes2, 3, all_edge_nodes)
        self.assertEqual(len(curved_polygons), 1)
        self.assertEqual(len(curved_polygons[0]), 3)
        # First triplet ("edge info")
        index, s_val, t_val = curved_polygons[0][0]
        self.assertEqual(index, 3)
        self._almost(s_val, 0.60937510406326101, 16)
        self._almost(t_val, 0.64417397312262581, 32)
        # Second triplet ("edge info")
        index, s_val, t_val = curved_polygons[0][1]
        self.assertEqual(index, 1)
        self._almost(s_val, 0.97192953004411253, 32)
        self.assertEqual(t_val, 1.0)
        # Third and final triplet ("edge info")
        index, s_val, t_val = curved_polygons[0][2]
        self.assertEqual(index, 2)
        self.assertEqual(s_val, 0.0)
        self._almost(t_val, 0.029255079571203865, 1024)


class Test_algebraic_intersect(Test__geometric_intersect):

    BAD_BOUNDARY_ARGS = ('Coincident curves not currently supported',)
    BAD_BOUNDARY_TYPE = RuntimeError
    BAD_BOUNDARY_INCREASE_ULPS = False

    @staticmethod
    def _call_function_under_test(nodes1, degree1, nodes2, degree2, **kwargs):
        from bezier import _surface_intersection

        return _surface_intersection.algebraic_intersect(
            nodes1, degree1, nodes2, degree2, **kwargs)

    @staticmethod
    def parallel_err():
        from bezier import _algebraic_intersection

        return _algebraic_intersection._COINCIDENT_ERR

    def _check_linear_intersection1(self, value):
        expected = 0.125
        delta = SPACING(expected)
        self.assertAlmostEqual(value, expected, delta=delta)

    def test_opposed_tangencies(self):
        from bezier import _algebraic_intersection

        with self.assertRaises(NotImplementedError) as exc_info:
            super(Test_algebraic_intersect, self).test_opposed_tangencies()

        exc_args = exc_info.exception.args
        self.assertEqual(len(exc_args), 2)
        self.assertEqual(exc_args[0], _algebraic_intersection._NON_SIMPLE_ERR)
        self.assertIsInstance(exc_args[1], np.ndarray)
        self.assertEqual(exc_args[1].shape, (3,))

    def test_tangent_contained(self):
        from bezier import _algebraic_intersection

        with self.assertRaises(NotImplementedError) as exc_info:
            super(Test_algebraic_intersect, self).test_tangent_contained()

        exc_args = exc_info.exception.args
        self.assertEqual(len(exc_args), 2)
        self.assertEqual(exc_args[0], _algebraic_intersection._NON_SIMPLE_ERR)
        self.assertIsInstance(exc_args[1], np.ndarray)
        self.assertEqual(exc_args[1].shape, (3,))


@utils.needs_speedup
class Test_speedup_geometric_intersect(Test__geometric_intersect):

    BAD_BOUNDARY_ARGS = ('Unknown error has occured.',)
    BAD_BOUNDARY_TYPE = RuntimeError
    BAD_BOUNDARY_INCREASE_ULPS = True

    @staticmethod
    def _call_function_under_test(nodes1, degree1, nodes2, degree2, **kwargs):
        from bezier import _surface_intersection

        return _surface_intersection.geometric_intersect(
            nodes1, degree1, nodes2, degree2, **kwargs)

    def test_two_curved_polygons(self):
        # Make sure there is enough space so that no resize is needed.
        sizes = surface_workspace_sizes()
        segment_ends_size, segments_size = sizes
        self.assertGreaterEqual(segment_ends_size, 2)
        self.assertGreaterEqual(segments_size, 6)

        super_ = super(Test_speedup_geometric_intersect, self)
        super_.test_two_curved_polygons()

        # Make sure the workspace was **not** resized.
        self.assertEqual(surface_workspace_sizes(), sizes)

    def test_resize_both(self):
        reset_surface_workspaces(segment_ends_size=1, segments_size=1)

        super_ = super(Test_speedup_geometric_intersect, self)
        super_.test_two_curved_polygons()

        # Make sure the sizes were resized from (1, 1).
        self.assertEqual(surface_workspace_sizes(), (2, 6))

    def test_insufficient_segment_ends(self):
        from bezier import _speedup

        reset_surface_workspaces(segment_ends_size=1)
        sizes = surface_workspace_sizes()

        with self.assertRaises(ValueError) as exc_info:
            self._two_curved_polygons(resizes_allowed=0)

        exc_args = exc_info.exception.args
        template = _speedup.SEGMENT_ENDS_TOO_SMALL
        self.assertEqual(exc_args, (template.format(2, 1),))
        # Make sure the workspace was **not** resized.
        self.assertEqual(surface_workspace_sizes(), sizes)

    def test_insufficient_segments(self):
        from bezier import _speedup

        reset_surface_workspaces(segment_ends_size=2, segments_size=2)
        sizes = surface_workspace_sizes()

        with self.assertRaises(ValueError) as exc_info:
            self._two_curved_polygons(resizes_allowed=0)

        exc_args = exc_info.exception.args
        template = _speedup.SEGMENTS_TOO_SMALL
        self.assertEqual(exc_args, (template.format(6, 2),))
        # Make sure the workspace was **not** resized.
        self.assertEqual(surface_workspace_sizes(), sizes)


@utils.needs_speedup
class Test_reset_surface_workspaces(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(**kwargs):
        return reset_surface_workspaces(**kwargs)

    def test_it(self):
        return_value = self._call_function_under_test(
            segment_ends_size=1, segments_size=2)
        self.assertIsNone(return_value)
        self.assertEqual(surface_workspace_sizes(), (1, 2))

    @unittest.expectedFailure
    def test_threadsafe(self):
        sizes_main = (4, 3)
        self._call_function_under_test(
            segment_ends_size=sizes_main[0], segments_size=sizes_main[1])

        worker = WorkspaceThreadedAccess()
        self.assertIsNone(worker.sizes1)
        self.assertIsNone(worker.sizes2)

        sizes1 = (1, 3)
        sizes2 = (2, 2)
        thread1 = threading.Thread(target=worker.task1, args=(sizes1,))
        thread2 = threading.Thread(target=worker.task2, args=(sizes2,))
        thread1.start()
        thread2.start()
        thread1.join()
        thread2.join()

        # This check demonstrates the **broken-ness** of the implementation.
        # The sizes for each thread should be the sizes actually **set** in
        # the given thread and the workspace in the main thread should be
        # unchanged (i.e. should have ``sizes_main``). What we'll actually
        # observe is ``(sizes2, sizes1, sizes2)``.
        expected = (sizes1, sizes2, sizes_main)
        actual = (
            worker.sizes1,
            worker.sizes2,
            surface_workspace_sizes(),
        )
        self.assertEqual(actual, expected)


@utils.needs_speedup
class Test_surface_workspace_sizes(unittest.TestCase):

    @staticmethod
    def _call_function_under_test():
        return surface_workspace_sizes()

    def test_it(self):
        reset_surface_workspaces(segment_ends_size=3, segments_size=5)
        self.assertEqual(self._call_function_under_test(), (3, 5))
        reset_surface_workspaces(segment_ends_size=1)
        self.assertEqual(self._call_function_under_test(), (1, 5))
        reset_surface_workspaces(segments_size=2)
        self.assertEqual(self._call_function_under_test(), (1, 2))


@utils.needs_speedup
class Test_speedup__type_info(unittest.TestCase):

    @staticmethod
    def _call_function_under_test():
        from bezier import _speedup

        return _speedup._type_info()

    def test_it(self):
        result = self._call_function_under_test()
        is_native, item_size, dtype_num, size_of_struct = result

        self.assertTrue(is_native)
        self.assertEqual(dtype_num, 20)
        if base_utils.IS_64_BIT or base_utils.IS_WINDOWS:
            self.assertEqual(item_size, 24)
            self.assertEqual(size_of_struct, 24)
        else:  # pragma: NO COVER
            self.assertEqual(item_size, 20)
            self.assertEqual(size_of_struct, 20)


def make_intersect(*args, **kwargs):
    from bezier import _intersection_helpers

    return _intersection_helpers.Intersection(*args, **kwargs)


def get_enum(str_val):
    from bezier import _surface_helpers

    return _surface_helpers.IntersectionClassification[str_val]


def check_edges(test_case, nodes1, degree1, nodes2, degree2, all_edge_nodes):
    from bezier import _surface_helpers

    test_case.assertIsInstance(all_edge_nodes, tuple)
    test_case.assertEqual(len(all_edge_nodes), 6)

    edge_nodes1 = _surface_helpers.compute_edge_nodes(nodes1, degree1)
    edge_nodes2 = _surface_helpers.compute_edge_nodes(nodes2, degree2)

    test_case.assertEqual(edge_nodes1[0], all_edge_nodes[0])
    test_case.assertEqual(edge_nodes1[1], all_edge_nodes[1])
    test_case.assertEqual(edge_nodes1[2], all_edge_nodes[2])
    test_case.assertEqual(edge_nodes2[0], all_edge_nodes[3])
    test_case.assertEqual(edge_nodes2[1], all_edge_nodes[4])
    test_case.assertEqual(edge_nodes2[2], all_edge_nodes[5])


def reset_surface_workspaces(**kwargs):
    from bezier import _speedup

    return _speedup.reset_surface_workspaces(**kwargs)


def surface_workspace_sizes():
    from bezier import _speedup

    return _speedup.surface_workspace_sizes()


class WorkspaceThreadedAccess(object):

    def __init__(self):
        self.barrier1 = threading.Event()
        self.barrier2 = threading.Event()
        self.barrier3 = threading.Event()
        self.sizes1 = None
        self.sizes2 = None

    def event1(self, sizes):
        # NOTE: There is no need to ``wait`` since this is the first event.
        reset_surface_workspaces(
            segment_ends_size=sizes[0], segments_size=sizes[1])
        self.barrier1.set()

    def event2(self):
        self.barrier1.wait()
        result = surface_workspace_sizes()
        self.barrier2.set()
        return result

    def event3(self, sizes):
        self.barrier2.wait()
        reset_surface_workspaces(
            segment_ends_size=sizes[0], segments_size=sizes[1])
        self.barrier3.set()

    def event4(self):
        self.barrier3.wait()
        # NOTE: There is no barrier to ``set`` since this is the last event.
        return surface_workspace_sizes()

    def task1(self, sizes):
        self.event1(sizes)
        self.sizes1 = self.event4()

    def task2(self, sizes):
        self.sizes2 = self.event2()
        self.event3(sizes)
