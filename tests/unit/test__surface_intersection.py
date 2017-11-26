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
            nodes1, 1, nodes2, 1, verify=True)
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


class Test_algebraic_intersect(Test__geometric_intersect):

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


@utils.needs_surface_intersection_speedup
class Test_speedup_geometric_intersect(Test__geometric_intersect):

    @staticmethod
    def _call_function_under_test(nodes1, degree1, nodes2, degree2, **kwargs):
        from bezier import _surface_intersection

        return _surface_intersection.geometric_intersect(
            nodes1, degree1, nodes2, degree2, **kwargs)

    @staticmethod
    def reset_workspaces(**kwargs):
        from bezier import _surface_intersection_speedup

        return _surface_intersection_speedup.reset_workspaces(**kwargs)

    @staticmethod
    def workspace_sizes():
        from bezier import _surface_intersection_speedup

        return _surface_intersection_speedup.workspace_sizes()

    def test_two_curved_polygons(self):
        # Make sure there is enough space so that no resize is needed.
        sizes = self.workspace_sizes()
        segment_ends_size, segments_size = sizes
        self.assertGreaterEqual(segment_ends_size, 2)
        self.assertGreaterEqual(segments_size, 6)

        super_ = super(Test_speedup_geometric_intersect, self)
        super_.test_two_curved_polygons()

        # Make sure the workspace was **not** resized.
        self.assertEqual(self.workspace_sizes(), sizes)

    def test_resize_both(self):
        self.reset_workspaces(segment_ends_size=1, segments_size=1)

        super_ = super(Test_speedup_geometric_intersect, self)
        super_.test_two_curved_polygons()

        # Make sure the sizes were resized from (1, 1).
        self.assertEqual(self.workspace_sizes(), (2, 6))

    def test_insufficient_segment_ends(self):
        from bezier import _surface_intersection_speedup

        self.reset_workspaces(segment_ends_size=1)
        sizes = self.workspace_sizes()

        with self.assertRaises(ValueError) as exc_info:
            self._two_curved_polygons(resizes_allowed=0)

        exc_args = exc_info.exception.args
        template = _surface_intersection_speedup.SEGMENT_ENDS_TOO_SMALL
        self.assertEqual(exc_args, (template.format(2, 1),))
        # Make sure the workspace was **not** resized.
        self.assertEqual(self.workspace_sizes(), sizes)

    def test_insufficient_segments(self):
        from bezier import _surface_intersection_speedup

        self.reset_workspaces(segment_ends_size=2, segments_size=2)
        sizes = self.workspace_sizes()

        with self.assertRaises(ValueError) as exc_info:
            self._two_curved_polygons(resizes_allowed=0)

        exc_args = exc_info.exception.args
        template = _surface_intersection_speedup.SEGMENTS_TOO_SMALL
        self.assertEqual(exc_args, (template.format(6, 2),))
        # Make sure the workspace was **not** resized.
        self.assertEqual(self.workspace_sizes(), sizes)


@utils.needs_surface_intersection_speedup
class Test_speedup__type_info(unittest.TestCase):

    @staticmethod
    def _call_function_under_test():
        from bezier import _surface_intersection_speedup

        return _surface_intersection_speedup._type_info()

    def test_it(self):
        expected = (True, 24, 20, 24)
        result = self._call_function_under_test()
        self.assertEqual(expected, result)


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
