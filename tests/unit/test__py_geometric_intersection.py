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

from tests import utils as base_utils
from tests.unit import utils


SPACING = np.spacing  # pylint: disable=no-member
UNIT_SQUARE = np.asfortranarray([[0.0, 1.0, 1.0, 0.0], [0.0, 0.0, 1.0, 1.0]])


class Test_bbox_intersect(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.bbox_intersect(nodes1, nodes2)

    def test_intersect(self):
        from bezier import _py_geometric_intersection

        nodes = UNIT_SQUARE + np.asfortranarray([[0.5], [0.5]])
        result = self._call_function_under_test(UNIT_SQUARE, nodes)
        expected = _py_geometric_intersection.BoxIntersectionType.INTERSECTION
        self.assertEqual(result, expected)

    def test_far_apart(self):
        from bezier import _py_geometric_intersection

        nodes = UNIT_SQUARE + np.asfortranarray([[100.0], [100.0]])
        result = self._call_function_under_test(UNIT_SQUARE, nodes)
        expected = _py_geometric_intersection.BoxIntersectionType.DISJOINT
        self.assertEqual(result, expected)

    def test_disjoint_but_aligned(self):
        from bezier import _py_geometric_intersection

        nodes = UNIT_SQUARE + np.asfortranarray([[1.0], [2.0]])
        result = self._call_function_under_test(UNIT_SQUARE, nodes)
        expected = _py_geometric_intersection.BoxIntersectionType.DISJOINT
        self.assertEqual(result, expected)

    def test_tangent(self):
        from bezier import _py_geometric_intersection

        nodes = UNIT_SQUARE + np.asfortranarray([[1.0], [0.0]])
        result = self._call_function_under_test(UNIT_SQUARE, nodes)
        expected = _py_geometric_intersection.BoxIntersectionType.TANGENT
        self.assertEqual(result, expected)

    def test_almost_tangent(self):
        from bezier import _py_geometric_intersection

        x_val = 1.0 + SPACING(1.0)
        nodes = UNIT_SQUARE + np.asfortranarray([[x_val], [0.0]])
        result = self._call_function_under_test(UNIT_SQUARE, nodes)
        expected = _py_geometric_intersection.BoxIntersectionType.DISJOINT
        self.assertEqual(result, expected)


class Test_linearization_error(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.linearization_error(nodes)

    def test_linear(self):
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 2.0]])
        error_val = self._call_function_under_test(nodes)
        self.assertEqual(error_val, 0.0)

    def test_degree_elevated_linear(self):
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 1.0, 2.0]])
        error_val = self._call_function_under_test(nodes)
        self.assertEqual(error_val, 0.0)
        nodes = np.asfortranarray(
            [[0.0, 0.25, 0.5, 0.75, 1.0], [0.0, 0.5, 1.0, 1.5, 2.0]]
        )
        error_val = self._call_function_under_test(nodes)
        self.assertEqual(error_val, 0.0)

    def test_hidden_linear(self):
        # NOTE: This is the line 3 y = 4 x, but with the parameterization
        #       x(s) = 3 s (4 - 3 s).
        nodes = np.asfortranarray([[0.0, 6.0, 3.0], [0.0, 8.0, 4.0]])
        error_val = self._call_function_under_test(nodes)
        # D^2 v = [-9, -12]
        expected = 0.125 * 2 * 1 * 15.0
        self.assertEqual(error_val, expected)

    def test_quadratic(self):
        from bezier.hazmat import curve_helpers

        nodes = np.asfortranarray([[0.0, 1.0, 5.0], [0.0, 1.0, 6.0]])
        # NOTE: This is hand picked so that
        #             d Nodes = [1, 1], [4, 5]
        #           d^2 Nodes = [3, 4]
        #       so that sqrt(3^2 + 4^2) = 5.0
        error_val = self._call_function_under_test(nodes)
        expected = 0.125 * 2 * 1 * 5.0
        self.assertEqual(error_val, expected)
        # For a degree two curve, the 2nd derivative is constant
        # so by subdividing, our error should drop by a factor
        # of (1/2)^2 = 4.
        left_nodes, right_nodes = curve_helpers.subdivide_nodes(nodes)
        error_left = self._call_function_under_test(left_nodes)
        error_right = self._call_function_under_test(right_nodes)
        self.assertEqual(error_left, 0.25 * expected)
        self.assertEqual(error_right, 0.25 * expected)

    def test_higher_dimension(self):
        nodes = np.asfortranarray(
            [[1.5, 3.5, 8.5], [0.0, -5.0, 2.0], [6.25, 10.25, 10.25]]
        )
        # NOTE: This is hand picked so that
        #             d Nodes = [2, -5, 4], [5, 7, 0]
        #           d^2 Nodes = [3, 12, -4]
        #       so that sqrt(3^2 + 12^2 + 4^2) = 13.0
        error_val = self._call_function_under_test(nodes)
        expected = 0.125 * 2 * 1 * 13.0
        self.assertEqual(error_val, expected)

    def test_hidden_quadratic(self):
        # NOTE: This is the quadratic y = 1 + x^2 / 4, but with the
        #       parameterization x(s) = (3 s - 1)^2.
        nodes = np.asfortranarray(
            [[1.0, -0.5, -0.5, 1.0, 4.0], [1.25, 0.5, 2.0, -1.0, 5.0]]
        )
        error_val = self._call_function_under_test(nodes)
        # D^2 v = [1.5, 2.25], [1.5, -4.5], [1.5, 9]
        expected = 0.125 * 4 * 3 * np.sqrt(1.5 ** 2 + 9.0 ** 2)
        self.assertEqual(expected, error_val)

    def test_cubic(self):
        nodes = np.asfortranarray([[0.0, 1.0, 5.0, 6.0], [0.0, 1.0, 6.0, 7.0]])
        # NOTE: This is hand picked so that
        #             d Nodes = [1, 1], [4, 5], [1, 1]
        #           d^2 Nodes = [3, 4], [-3, -4]
        #       so that sqrt(3^2 + 4^2) = 5.0
        error_val = self._call_function_under_test(nodes)
        expected = 0.125 * 3 * 2 * 5.0
        self.assertEqual(error_val, expected)

    def test_quartic(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0, 5.0, 6.0, 4.0], [0.0, 1.0, 6.0, 7.0, 7.0]]
        )
        # NOTE: This is hand picked so that
        #             d Nodes = [1, 1], [4, 5], [1, 1], [-2, 0]
        #           d^2 Nodes = [3, 4], [-3, -4], [-3, -1]
        #       so that sqrt(3^2 + 4^2) = 5.0
        error_val = self._call_function_under_test(nodes)
        expected = 0.125 * 4 * 3 * 5.0
        self.assertEqual(error_val, expected)

    def test_degree_weights_on_the_fly(self):
        nodes = np.asfortranarray(
            [
                [0.0, 1.0, 7.0, 11.0, 15.0, 16.0],
                [0.0, 1.0, 3.0, 8.0, 1.0, -3.0],
            ]
        )
        # NOTE: This is hand picked so that
        #             d Nodes = [1, 1], [6, 2], [4, 5], [4, -7], [1, -4]
        #           d^2 Nodes = [5, 1], [-2, 3], [0, -12], [-3, 3]
        #       so that sqrt(5^2 + 12^2) = 13.0
        error_val = self._call_function_under_test(nodes)
        expected = 0.125 * 5 * 4 * 13.0
        self.assertEqual(error_val, expected)


class Test_segment_intersection(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(start0, end0, start1, end1):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.segment_intersection(
            start0, end0, start1, end1
        )

    def _helper(
        self, intersection, s_val, direction0, t_val, direction1, **kwargs
    ):
        start0 = intersection + s_val * direction0
        end0 = intersection + (s_val - 1.0) * direction0
        start1 = intersection + t_val * direction1
        end1 = intersection + (t_val - 1.0) * direction1
        return self._call_function_under_test(
            start0, end0, start1, end1, **kwargs
        )

    def test_success(self):
        intersection = np.asfortranarray([1.0, 2.0])
        s_val = 0.25
        t_val = 0.625
        direction0 = np.asfortranarray([3.0, 0.5])
        direction1 = np.asfortranarray([-2.0, 1.0])
        # D0 x D1 == 4.0, so there will be no round-off in answer.
        computed_s, computed_t, success = self._helper(
            intersection, s_val, direction0, t_val, direction1
        )
        self.assertEqual(computed_s, s_val)
        self.assertEqual(computed_t, t_val)
        self.assertTrue(success)

    def test_parallel(self):
        intersection = np.asfortranarray([0.0, 0.0])
        s_val = 0.5
        t_val = 0.5
        direction0 = np.asfortranarray([0.0, 1.0])
        direction1 = np.asfortranarray([0.0, 2.0])
        computed_s, computed_t, success = self._helper(
            intersection, s_val, direction0, t_val, direction1
        )
        self.assertIsNone(computed_s)
        self.assertIsNone(computed_t)
        self.assertFalse(success)


class Test_parallel_lines_parameters(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(start0, end0, start1, end1):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.parallel_lines_parameters(
            start0, end0, start1, end1
        )

    def test_same_line_no_overlap(self):
        start0 = np.asfortranarray([0.0, 0.0])
        end0 = np.asfortranarray([3.0, 4.0])
        start1 = np.asfortranarray([6.0, 8.0])
        end1 = np.asfortranarray([9.0, 12.0])
        disjoint, params = self._call_function_under_test(
            start0, end0, start1, end1
        )
        self.assertTrue(disjoint)
        self.assertIsNone(params)
        # Do the same, but reverse the direction of line 1.
        disjoint, params = self._call_function_under_test(
            start0, end0, end1, start1
        )
        self.assertTrue(disjoint)
        self.assertIsNone(params)
        # Do the same, but swap the lines.
        disjoint, params = self._call_function_under_test(
            start1, end1, start0, end0
        )
        self.assertTrue(disjoint)
        self.assertIsNone(params)
        # Do the same, but reverse the direction of line 0.
        disjoint, params = self._call_function_under_test(
            end0, start0, start1, end1
        )
        self.assertTrue(disjoint)
        self.assertIsNone(params)

    def test_same_line_overlap_at_start(self):
        start0 = np.asfortranarray([6.0, -3.0])
        end0 = np.asfortranarray([-7.0, 1.0])
        start1 = np.asfortranarray([1.125, -1.5])
        end1 = np.asfortranarray([-5.375, 0.5])
        disjoint, params = self._call_function_under_test(
            start0, end0, start1, end1
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[0.375, 0.875], [0.0, 1.0]])
        self.assertEqual(params, expected)
        # Do the same, but reverse the direction of line 1.
        disjoint, params = self._call_function_under_test(
            start0, end0, end1, start1
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[0.875, 0.375], [0.0, 1.0]])
        self.assertEqual(params, expected)
        # Do the same, but swap the lines.
        disjoint, params = self._call_function_under_test(
            start1, end1, start0, end0
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[0.0, 1.0], [0.375, 0.875]])
        self.assertEqual(params, expected)

    def test_same_line_overlap_at_end(self):
        start0 = np.asfortranarray([1.0, 2.0])
        end0 = np.asfortranarray([3.0, 5.0])
        start1 = np.asfortranarray([-2.0, -2.5])
        end1 = np.asfortranarray([2.0, 3.5])
        disjoint, params = self._call_function_under_test(
            start0, end0, start1, end1
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[0.0, 0.5], [0.75, 1.0]])
        self.assertEqual(params, expected)
        # Do the same, but reverse the direction of line 1.
        disjoint, params = self._call_function_under_test(
            start0, end0, end1, start1
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[0.5, 0.0], [0.0, 0.25]])
        self.assertEqual(params, expected)
        # Do the same, but swap the lines.
        disjoint, params = self._call_function_under_test(
            start1, end1, start0, end0
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[0.75, 1.0], [0.0, 0.5]])
        self.assertEqual(params, expected)

    def test_same_line_contained(self):
        start0 = np.asfortranarray([-9.0, 0.0])
        end0 = np.asfortranarray([4.0, 5.0])
        start1 = np.asfortranarray([23.5, 12.5])
        end1 = np.asfortranarray([-28.5, -7.5])
        disjoint, params = self._call_function_under_test(
            start0, end0, start1, end1
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[1.0, 0.0], [0.375, 0.625]])
        self.assertEqual(params, expected)
        # Do the same, but reverse the direction of line 1.
        disjoint, params = self._call_function_under_test(
            start0, end0, end1, start1
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[0.0, 1.0], [0.375, 0.625]])
        self.assertEqual(params, expected)
        # Do the same, but swap the lines.
        disjoint, params = self._call_function_under_test(
            start1, end1, start0, end0
        )
        self.assertFalse(disjoint)
        expected = np.asfortranarray([[0.625, 0.375], [0.0, 1.0]])
        self.assertEqual(params, expected)

    def test_different_line(self):
        start0 = np.asfortranarray([3.0, 2.0])
        end0 = np.asfortranarray([3.0, 0.75])
        start1 = np.asfortranarray([0.0, 0.0])
        end1 = np.asfortranarray([0.0, 2.0])
        disjoint, params = self._call_function_under_test(
            start0, end0, start1, end1
        )
        self.assertTrue(disjoint)
        self.assertIsNone(params)


class Test_line_line_collide(unittest.TestCase):
    FIRST = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])

    @staticmethod
    def _call_function_under_test(line1, line2):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.line_line_collide(line1, line2)

    def test_general_position_collide(self):
        # I.e. not parallel.
        line2 = np.asfortranarray([[0.0, 1.0], [1.0, 0.0]])
        do_collide = self._call_function_under_test(self.FIRST, line2)
        self.assertTrue(do_collide)

    def test_general_position_miss(self):
        # I.e. not parallel.
        line2 = np.asfortranarray([[3.0, 4.0], [3.0, 5.0]])
        do_collide = self._call_function_under_test(self.FIRST, line2)
        self.assertFalse(do_collide)

    def test_parallel(self):
        line2 = np.asfortranarray([[0.0, 1.0], [1.0, 2.0]])
        do_collide = self._call_function_under_test(self.FIRST, line2)
        self.assertFalse(do_collide)

    def test_parallel_disjoint_same_line(self):
        line2 = np.asfortranarray([[2.0, 3.0], [2.0, 3.0]])
        do_collide = self._call_function_under_test(self.FIRST, line2)
        self.assertFalse(do_collide)

    def test_parallel_overlap(self):
        line2 = np.asfortranarray([[0.5, 1.5], [0.5, 1.5]])
        do_collide = self._call_function_under_test(self.FIRST, line2)
        self.assertTrue(do_collide)


class Test_convex_hull_collide(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.convex_hull_collide(nodes1, nodes2)

    def test_line_line(self):
        nodes1 = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 2.0]])
        nodes2 = np.asfortranarray([[0.0, 0.0], [1.0, 2.0]])
        do_collide = self._call_function_under_test(nodes1, nodes2)
        self.assertFalse(do_collide)

    def test_hull_hull(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        nodes2 = np.asfortranarray([[0.0, 1.0, 2.0], [0.5, 1.0, 0.5]])
        do_collide = self._call_function_under_test(nodes1, nodes2)
        self.assertTrue(do_collide)


class Test_from_linearized(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(first, second, intersections):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.from_linearized(
            first, second, intersections
        )

    def test_does_intersect(self):
        nodes1 = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 1.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        # NOTE: This curve isn't close to linear, but that's OK.
        lin1 = make_linearization(curve1, error=np.nan)
        nodes2 = np.asfortranarray([[0.0, 0.5, 1.0], [1.0, 1.0, 0.0]])
        curve2 = subdivided_curve(nodes2)
        # NOTE: This curve isn't close to linear, but that's OK.
        lin2 = make_linearization(curve2, error=np.nan)
        intersections = []
        self.assertIsNone(
            self._call_function_under_test(lin1, lin2, intersections)
        )
        self.assertEqual(intersections, [(0.5, 0.5)])

    def test_lines_outside_unit_interval_disjoint_hulls(self):
        # The bounding boxes intersect but the lines do not.
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.0)
        nodes2 = np.asfortranarray([[1.75, 0.75], [-0.75, 0.25]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        intersections = []
        self.assertIsNone(
            self._call_function_under_test(lin1, lin2, intersections)
        )
        self.assertEqual(intersections, [])

    def test_curved_outside_unit_interval_disjoint_hulls(self):
        # The bounding boxes intersect but the lines do not.
        nodes1 = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.25)
        nodes2 = np.asfortranarray([[1.75, 1.25, 0.75], [-0.75, -0.75, 0.25]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.25)
        intersections = []
        return_value = self._call_function_under_test(
            lin1, lin2, intersections
        )
        self.assertIsNone(return_value)
        self.assertEqual(intersections, [])

    def test_unhandled_parallel_lines(self):
        from bezier import _py_geometric_intersection

        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.0)
        nodes2 = np.asfortranarray([[0.0, 1.0], [1.0, 2.0]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        intersections = []
        with self.assertRaises(ValueError) as exc_info:
            self._call_function_under_test(lin1, lin2, intersections)
        expected_args = (_py_geometric_intersection._UNHANDLED_LINES,)
        self.assertEqual(exc_info.exception.args, expected_args)
        self.assertEqual(intersections, [])

    def test_curved_parallel_segments_disjoint_hulls(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.0)
        # This is a "bad" parameterization of ``y == x``.
        nodes2 = np.asfortranarray(
            [[2.0, 2.5009765625, 3.0], [2.0, 2.5009765625, 3.0]]
        )
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=np.nan)
        intersections = []
        return_value = self._call_function_under_test(
            lin1, lin2, intersections
        )
        self.assertIsNone(return_value)
        self.assertEqual(intersections, [])

    def test_curved_parallel_segments_tangent(self):
        # NOTE: There is no corresponding "enable", but the disable only
        #       applies in this lexical scope.
        # pylint: disable=too-many-locals
        from bezier.hazmat import curve_helpers

        start1 = 5461.0 / 8192.0
        end1 = 5462.0 / 8192.0
        original_nodes1 = np.asfortranarray(
            [[0.0, 0.375, 0.75], [0.0, 0.75, 0.375]]
        )
        nodes1 = curve_helpers.specialize_curve(original_nodes1, start1, end1)
        curve1 = subdivided_curve(
            nodes1, original_nodes=original_nodes1, start=start1, end=end1
        )
        lin1 = make_linearization(curve1)
        start2 = 2730.0 / 8192.0
        end2 = 2731.0 / 8192.0
        original_nodes2 = np.asfortranarray(
            [[0.25, 0.625, 1.0], [0.625, 0.25, 1.0]]
        )
        nodes2 = curve_helpers.specialize_curve(original_nodes2, start2, end2)
        curve2 = subdivided_curve(
            nodes2, original_nodes=original_nodes2, start=start2, end=end2
        )
        lin2 = make_linearization(curve2)
        intersections = []
        return_value = self._call_function_under_test(
            lin1, lin2, intersections
        )
        self.assertIsNone(return_value)
        self.assertEqual(len(intersections), 1)
        s, t = intersections[0]
        utils.almost(self, 2.0 / 3.0, s, 1)
        utils.almost(self, 1.0 / 3.0, t, 1)

    def test_line_and_curve_tangent(self):
        # NOTE: There is no corresponding "enable", but the disable only
        #       applies in this lexical scope.
        # pylint: disable=too-many-locals
        from bezier.hazmat import curve_helpers

        # B1([5461/16384, 5462/16384]) and B2([0, 1]) are linearized
        # and when the segments intersect they produce s = -1/3 < 0.
        start1 = 5461.0 / 16384.0
        end1 = 5462.0 / 16384.0
        original_nodes1 = np.asfortranarray(
            [[0.0, 1.5, 3.0], [2.25, -2.25, 2.25]]
        )
        nodes1 = curve_helpers.specialize_curve(original_nodes1, start1, end1)
        curve1 = subdivided_curve(
            nodes1, original_nodes=original_nodes1, start=start1, end=end1
        )
        lin1 = make_linearization(curve1)
        nodes2 = np.asfortranarray([[-0.5, 4.0], [1.75, -2.75]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        intersections = []
        return_value = self._call_function_under_test(
            lin1, lin2, intersections
        )
        self.assertIsNone(return_value)
        self.assertEqual(len(intersections), 1)
        s, t = intersections[0]
        utils.almost(self, 1.0 / 3.0, s, 1)
        utils.almost(self, 1.0 / 3.0, t, 1)

    def _wiggle_outside_helper(self, swap=False):
        from bezier.hazmat import curve_helpers

        # B1([127/128, 1]) and B2([0, 1]) are linearized and when the segments
        # intersect they produce s = 0.9999999999940287 and
        # t = 0.35294117647054574. After performing Newton's method the final
        # parameter values are s = 1.000000000000334, t = 0.3529411764708877,
        # the first of which is outside of ``[-2^{-44}, 1 + 2^{-44}]``.
        start1 = 127.0 / 128.0
        end1 = 1.0
        original_nodes1 = np.asfortranarray(
            [
                [
                    -0.8558466524001767,
                    -0.8384755863215392,
                    -0.8214285714285716,
                ],
                [-0.9648818559703933, -0.9825154456957814, -1.0],
            ]
        )
        nodes1 = curve_helpers.specialize_curve(original_nodes1, start1, end1)
        curve1 = subdivided_curve(
            nodes1, original_nodes=original_nodes1, start=start1, end=end1
        )
        lin1 = make_linearization(curve1)
        nodes2 = np.asfortranarray(
            [
                [
                    -0.8348214285714286,
                    -0.8158482142857144,
                    -0.7968750000000001,
                ],
                [-0.986607142857143, -1.0055803571428572, -1.0245535714285714],
            ]
        )
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2)
        intersections = []
        if swap:
            lin1, lin2 = lin2, lin1
        return_value = self._call_function_under_test(
            lin1, lin2, intersections
        )
        self.assertIsNone(return_value)
        self.assertEqual(intersections, [])

    @unittest.skipIf(not base_utils.IS_64_BIT, "32-bit is skipped")
    def test_s_wiggle_outside(self):
        self._wiggle_outside_helper()

    @unittest.skipIf(not base_utils.IS_64_BIT, "32-bit is skipped")
    def test_t_wiggle_outside(self):
        self._wiggle_outside_helper(swap=True)


class Test_add_intersection(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(s, t, intersections):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.add_intersection(s, t, intersections)

    def test_new(self):
        intersections = [(0.5, 0.5)]
        self.assertIsNone(
            self._call_function_under_test(0.75, 0.25, intersections)
        )
        expected = [(0.5, 0.5), (0.75, 0.25)]
        self.assertEqual(intersections, expected)

    def test_existing(self):
        intersections = [(0.25, 1.0)]
        self.assertIsNone(
            self._call_function_under_test(0.25, 1.0, intersections)
        )
        self.assertEqual(intersections, [(0.25, 1.0)])

    def test_s_under_zero(self):
        intersections = [(0.0, 0.75)]
        candidate_s = 0.5 ** 43
        self.assertIsNone(
            self._call_function_under_test(candidate_s, 0.75, intersections)
        )
        self.assertEqual(intersections, [(0.0, 0.75)])

    def test_t_under_zero(self):
        intersections = [(0.125, 0.0)]
        candidate_t = 0.5 ** 42
        self.assertIsNone(
            self._call_function_under_test(0.125, candidate_t, intersections)
        )
        self.assertEqual(intersections, [(0.125, 0.0)])


class Test_endpoint_check(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(
        first, node_first, s, second, node_second, t, intersections
    ):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.endpoint_check(
            first, node_first, s, second, node_second, t, intersections
        )

    def test_not_close(self):
        node_first = np.asfortranarray([0.0, 0.0])
        node_second = np.asfortranarray([1.0, 1.0])
        intersections = []
        self._call_function_under_test(
            None, node_first, None, None, node_second, None, intersections
        )
        self.assertEqual(intersections, [])

    def test_same(self):
        nodes_first = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        first = subdivided_curve(nodes_first)
        nodes_second = np.asfortranarray([[1.0, 2.0], [1.0, 1.0]])
        second = subdivided_curve(nodes_second)
        s_val = 1.0
        node_first = first.nodes[:, 1]
        t_val = 0.0
        node_second = second.nodes[:, 0]
        intersections = []
        self._call_function_under_test(
            first, node_first, s_val, second, node_second, t_val, intersections
        )
        self.assertEqual(intersections, [(s_val, t_val)])

    def test_subcurves_middle(self):
        nodes1 = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 1.0, 0.0]])
        root1 = subdivided_curve(nodes1)
        first, _ = root1.subdivide()
        nodes2 = np.asfortranarray([[1.0, 0.0, 1.0], [1.5, 0.5, -0.5]])
        root2 = subdivided_curve(nodes2)
        _, second = root2.subdivide()
        s_val = 1.0
        node_first = first.nodes[:, 2]
        t_val = 0.0
        node_second = second.nodes[:, 0]
        intersections = []
        self._call_function_under_test(
            first, node_first, s_val, second, node_second, t_val, intersections
        )
        self.assertEqual(intersections, [(0.5, 0.5)])


class Test_tangent_bbox_intersection(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(first, second, intersections):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.tangent_bbox_intersection(
            first, second, intersections
        )

    def test_one_endpoint(self):
        nodes1 = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 2.0, 0.0]])
        curve1 = subdivided_curve(nodes1)
        nodes2 = np.asfortranarray([[2.0, 3.0, 4.0], [0.0, 2.0, 0.0]])
        curve2 = subdivided_curve(nodes2)
        intersections = []
        self.assertIsNone(
            self._call_function_under_test(curve1, curve2, intersections)
        )
        self.assertEqual(intersections, [(1.0, 0.0)])

    def test_two_endpoints(self):
        nodes1 = np.asfortranarray([[0.0, -1.0, 0.0], [0.0, 0.5, 1.0]])
        curve1 = subdivided_curve(nodes1)
        nodes2 = np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.5, 1.0]])
        curve2 = subdivided_curve(nodes2)
        intersections = []
        self.assertIsNone(
            self._call_function_under_test(curve1, curve2, intersections)
        )
        expected = [(0.0, 0.0), (1.0, 1.0)]
        self.assertEqual(intersections, expected)

    def test_no_endpoints(self):
        # Lines have tangent bounding boxes but don't intersect.
        nodes1 = np.asfortranarray([[0.0, 2.0], [0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        nodes2 = np.asfortranarray([[0.5, 2.5], [1.0, 2.0]])
        curve2 = subdivided_curve(nodes2)
        intersections = []
        self.assertIsNone(
            self._call_function_under_test(curve1, curve2, intersections)
        )
        self.assertEqual(intersections, [])


class Test_bbox_line_intersect(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes, line_start, line_end):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.bbox_line_intersect(
            nodes, line_start, line_end
        )

    def test_start_in_bbox(self):
        from bezier import _py_geometric_intersection

        line_start = np.asfortranarray([0.5, 0.5])
        line_end = np.asfortranarray([0.5, 1.5])
        result = self._call_function_under_test(
            UNIT_SQUARE, line_start, line_end
        )
        expected = _py_geometric_intersection.BoxIntersectionType.INTERSECTION
        self.assertEqual(result, expected)

    def test_end_in_bbox(self):
        from bezier import _py_geometric_intersection

        line_start = np.asfortranarray([-1.0, 0.5])
        line_end = np.asfortranarray([0.5, 0.5])
        result = self._call_function_under_test(
            UNIT_SQUARE, line_start, line_end
        )
        expected = _py_geometric_intersection.BoxIntersectionType.INTERSECTION
        self.assertEqual(result, expected)

    def test_segment_intersect_bbox_bottom(self):
        from bezier import _py_geometric_intersection

        line_start = np.asfortranarray([0.5, -0.5])
        line_end = np.asfortranarray([0.5, 1.5])
        result = self._call_function_under_test(
            UNIT_SQUARE, line_start, line_end
        )
        expected = _py_geometric_intersection.BoxIntersectionType.INTERSECTION
        self.assertEqual(result, expected)

    def test_segment_intersect_bbox_right(self):
        from bezier import _py_geometric_intersection

        line_start = np.asfortranarray([-0.5, 0.5])
        line_end = np.asfortranarray([1.5, 0.5])
        result = self._call_function_under_test(
            UNIT_SQUARE, line_start, line_end
        )
        expected = _py_geometric_intersection.BoxIntersectionType.INTERSECTION
        self.assertEqual(result, expected)

    def test_segment_intersect_bbox_top(self):
        from bezier import _py_geometric_intersection

        line_start = np.asfortranarray([-0.25, 0.5])
        line_end = np.asfortranarray([0.5, 1.25])
        result = self._call_function_under_test(
            UNIT_SQUARE, line_start, line_end
        )
        expected = _py_geometric_intersection.BoxIntersectionType.INTERSECTION
        self.assertEqual(result, expected)

    def test_disjoint(self):
        from bezier import _py_geometric_intersection

        line_start = np.asfortranarray([2.0, 2.0])
        line_end = np.asfortranarray([2.0, 5.0])
        result = self._call_function_under_test(
            UNIT_SQUARE, line_start, line_end
        )
        expected = _py_geometric_intersection.BoxIntersectionType.DISJOINT
        self.assertEqual(result, expected)


class Test_intersect_one_round(utils.NumPyTestCase):
    # NOTE: QUADRATIC1 is a specialization of [0, 0], [1/2, 1], [1, 1]
    #       onto the interval [1/4, 1].
    QUADRATIC1 = np.asfortranarray([[0.25, 0.625, 1.0], [0.4375, 1.0, 1.0]])
    # NOTE: QUADRATIC2 is a specialization of [0, 1], [1/2, 1], [1, 0]
    #       onto the interval [0, 3/4].
    QUADRATIC2 = np.asfortranarray([[0.0, 0.375, 0.75], [1.0, 1.0, 0.4375]])
    LINE1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
    LINE2 = np.asfortranarray([[0.0, 1.0], [1.0, 0.0]])

    @staticmethod
    def _call_function_under_test(candidates, intersections):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.intersect_one_round(
            candidates, intersections
        )

    def _curves_compare(self, curve1, curve2):
        from bezier import _py_geometric_intersection

        if isinstance(curve1, _py_geometric_intersection.Linearization):
            self.assertIsInstance(
                curve2, _py_geometric_intersection.Linearization
            )
            # We just check identity, since we assume a ``Linearization``
            # can't be subdivided.
            self.assertIs(curve1, curve2)
        else:
            self.assertIsInstance(
                curve1, _py_geometric_intersection.SubdividedCurve
            )
            self.assertIsInstance(
                curve2, _py_geometric_intersection.SubdividedCurve
            )
            self.assertIs(curve1.original_nodes, curve2.original_nodes)
            self.assertEqual(curve1.start, curve2.start)
            self.assertEqual(curve1.end, curve2.end)
            self.assertEqual(curve1.nodes, curve2.nodes)

    def _candidates_compare(self, actual, expected):
        self.assertEqual(len(actual), len(expected))
        for first, second in zip(actual, expected):
            self.assertEqual(len(first), 2)
            self.assertEqual(len(second), 2)
            self._curves_compare(first[0], second[0])
            self._curves_compare(first[1], second[1])

    def test_simple(self):
        curve1 = subdivided_curve(self.QUADRATIC1)
        curve2 = subdivided_curve(self.QUADRATIC2)
        candidates = [(curve1, curve2)]
        next_candidates = self._call_function_under_test(candidates, [])
        left1, right1 = curve1.subdivide()
        left2, right2 = curve2.subdivide()
        expected = [
            (left1, left2),
            (left1, right2),
            (right1, left2),
            (right1, right2),
        ]
        self._candidates_compare(next_candidates, expected)

    def test_first_linearized(self):
        curve1 = subdivided_curve(self.LINE1)
        lin1 = make_linearization(curve1, error=0.0)
        curve2 = subdivided_curve(self.QUADRATIC2)
        intersections = []
        next_candidates = self._call_function_under_test(
            [(lin1, curve2)], intersections
        )
        self.assertEqual(intersections, [])
        left2, right2 = curve2.subdivide()
        expected = [(lin1, left2), (lin1, right2)]
        self._candidates_compare(next_candidates, expected)

    def test_second_linearized(self):
        curve1 = subdivided_curve(self.QUADRATIC1)
        curve2 = subdivided_curve(self.LINE2)
        lin2 = make_linearization(curve2, error=0.0)
        intersections = []
        next_candidates = self._call_function_under_test(
            [(curve1, lin2)], intersections
        )
        self.assertEqual(intersections, [])
        left1, right1 = curve1.subdivide()
        expected = [(left1, lin2), (right1, lin2)]
        self._candidates_compare(next_candidates, expected)

    def test_both_linearized(self):
        curve1 = subdivided_curve(self.LINE1)
        lin1 = make_linearization(curve1, error=0.0)
        curve2 = subdivided_curve(self.LINE2)
        lin2 = make_linearization(curve2, error=0.0)
        intersections = []
        next_candidates = self._call_function_under_test(
            [(lin1, lin2)], intersections
        )
        self.assertEqual(next_candidates, [])
        self.assertEqual(intersections, [(0.5, 0.5)])

    def test_failure_due_to_unhandled_lines(self):
        from bezier import _py_geometric_intersection

        curve1 = subdivided_curve(self.LINE1)
        lin1 = make_linearization(curve1, error=0.0)
        nodes2 = np.asfortranarray([[0.5, 4.5], [0.5, 4.5]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        intersections = []
        with self.assertRaises(ValueError) as exc_info:
            self._call_function_under_test([(lin1, lin2)], intersections)
        expected_args = (_py_geometric_intersection._UNHANDLED_LINES,)
        self.assertEqual(exc_info.exception.args, expected_args)
        self.assertEqual(intersections, [])

    def test_disjoint_bboxes(self):
        curve1 = subdivided_curve(self.QUADRATIC1)
        nodes2 = np.asfortranarray([[1.0, 0.0], [1.25, 2.0]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        intersections = []
        next_candidates = self._call_function_under_test(
            [(curve1, lin2)], intersections
        )
        self.assertEqual(next_candidates, [])
        self.assertEqual(intersections, [])

    def test_tangent_bboxes(self):
        nodes1 = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 1.0, 0.0]])
        curve1 = subdivided_curve(nodes1)
        nodes2 = np.asfortranarray([[1.0, 1.5, 2.0], [0.0, 0.5, -0.25]])
        curve2 = subdivided_curve(nodes2)
        intersections = []
        next_candidates = self._call_function_under_test(
            [(curve1, curve2)], intersections
        )
        self.assertEqual(next_candidates, [])
        self.assertEqual(intersections, [(1.0, 0.0)])


class Test_prune_candidates(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(candidates):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.prune_candidates(candidates)

    def test_curved(self):
        nodes1 = np.asfortranarray([[0.0, 2.0, 2.0], [0.0, 0.0, 2.0]])
        curve1 = subdivided_curve(nodes1)
        nodes2 = np.asfortranarray([[0.0, 0.0, 2.0], [1.0, 3.0, 3.0]])
        curve2 = subdivided_curve(nodes2)
        candidates = [(curve1, curve2)]
        pruned = self._call_function_under_test(candidates)
        self.assertEqual(pruned, [])

    def test_linear(self):
        nodes1 = np.asfortranarray([[0.0, 4.0], [0.0, 4.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.0)
        nodes2 = np.asfortranarray([[2.0, 6.0], [1.0, 1.0]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        candidates = [(lin1, lin2)]
        pruned = self._call_function_under_test(candidates)
        self.assertEqual(pruned, [])


class Test_make_same_degree(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.make_same_degree(nodes1, nodes2)

    def test_same_degree(self):
        nodes1 = np.asfortranarray([[0.0, 1.0, 2.0], [1.0, 2.0, 1.0]])
        nodes2 = np.asfortranarray([[1.0, 2.0, 4.0], [2.0, 1.0, 2.0]])
        elevated1, elevated2 = self._call_function_under_test(nodes1, nodes2)
        self.assertIs(elevated1, nodes1)
        self.assertIs(elevated2, nodes2)

    def test_elevate_once(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        nodes2 = np.asfortranarray([[1.0, 2.0, 0.0], [2.0, 2.0, 0.0]])
        # Elevate as the first argument.
        elevated1, elevated2 = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 0.5, 1.0]])
        self.assertEqual(elevated1, expected)
        self.assertIs(elevated2, nodes2)
        # Elevate as the second argument.
        elevated1, elevated2 = self._call_function_under_test(nodes2, nodes1)
        self.assertIs(elevated1, nodes2)
        self.assertEqual(elevated2, expected)

    def test_elevate_twice(self):
        nodes1 = np.asfortranarray([[0.0, 3.0], [0.0, 3.0]])
        nodes2 = np.asfortranarray(
            [[0.0, 1.0, 3.0, 4.0], [1.0, 2.0, 2.0, 2.0]]
        )
        # Elevate as the first argument.
        elevated1, elevated2 = self._call_function_under_test(nodes1, nodes2)
        expected = np.asfortranarray(
            [[0.0, 1.0, 2.0, 3.0], [0.0, 1.0, 2.0, 3.0]]
        )
        self.assertEqual(elevated1, expected)
        self.assertIs(elevated2, nodes2)
        # Elevate as the second argument.
        elevated1, elevated2 = self._call_function_under_test(nodes2, nodes1)
        self.assertIs(elevated1, nodes2)
        self.assertEqual(elevated2, expected)


class Test_coincident_parameters(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.coincident_parameters(nodes1, nodes2)

    def test_different_degree(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        nodes2 = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 0.0]])
        result = self._call_function_under_test(nodes1, nodes2)
        self.assertIsNone(result)

    def test_elevated_degree(self):
        from bezier.hazmat import curve_helpers

        nodes1 = np.asfortranarray([[0.0, 3.0, 6.0], [0.0, 3.0, 0.0]])
        nodes2 = curve_helpers.elevate_nodes(nodes1)
        result = self._call_function_under_test(nodes1, nodes2)
        expected = ((0.0, 0.0), (1.0, 1.0))
        self.assertEqual(result, expected)

    def test_invalid_point(self):
        nodes1 = np.asfortranarray(
            [[0.0, -8.0, 8.0, -6.0], [16.0, 0.0, 8.0, 13.0]]
        )
        nodes2 = np.asfortranarray([[-2.0, -2.0], [11.0, 16.0]])
        with self.assertRaises(ValueError) as exc_info1:
            self._call_function_under_test(nodes1, nodes2)
        with self.assertRaises(ValueError) as exc_info2:
            self._call_function_under_test(nodes2, nodes1)
        self.assertEqual(
            exc_info1.exception.args[0],
            "Parameters not close enough to one another",
        )
        self.assertEqual(
            exc_info2.exception.args[0],
            "Parameters not close enough to one another",
        )

    def test_touch_no_intersect_same_curve(self):
        from bezier.hazmat import curve_helpers

        nodes1 = np.asfortranarray([[0.0, 1.0, 3.0], [0.0, 2.0, 2.0]])
        new_params = ((1.0, 2.0), (2.0, 1.0), (-1.0, 0.0), (0.0, -1.0))
        for start, end in new_params:
            nodes2 = curve_helpers.specialize_curve(nodes1, start, end)
            result = self._call_function_under_test(nodes1, nodes2)
            self.assertIsNone(result)

    def test_disjoint_segments_same_curve(self):
        from bezier.hazmat import curve_helpers

        nodes1 = np.asfortranarray(
            [[1.0, 1.0, 3.0, 4.0], [0.0, 2.0, 2.0, 0.0]]
        )
        new_params = ((1.5, 2.0), (2.0, 1.5), (-1.0, -0.5), (-0.5, -1.0))
        for start, end in new_params:
            nodes2 = curve_helpers.specialize_curve(nodes1, start, end)
            result = self._call_function_under_test(nodes1, nodes2)
            self.assertIsNone(result)

    def test_ends_touch_different_curve(self):
        nodes1 = np.asfortranarray([[1.0, 3.0, 3.0], [8.0, 6.0, 5.0]])
        nodes2 = np.asfortranarray([[3.0, 3.0, 4.0], [5.0, -1.0, -4.0]])
        result = self._call_function_under_test(nodes1, nodes2)
        self.assertIsNone(result)

    def test_identical_curves(self):
        nodes = np.asfortranarray([[0.0, 2.0, 5.0], [0.0, 3.0, 5.0]])
        result = self._call_function_under_test(nodes, nodes)
        expected = ((0.0, 0.0), (1.0, 1.0))
        self.assertEqual(result, expected)

    def test_contained_and_touching(self):
        from bezier.hazmat import curve_helpers

        nodes1 = np.asfortranarray([[4.0, 6.0, 2.0], [1.0, 3.0, 1.0]])
        new_params = (
            (0.0, 0.5),
            (0.5, 0.0),
            (0.5, 1.0),
            (1.0, 0.5),
            (0.0, 2.0),
            (2.0, 0.0),
            (-1.0, 1.0),
            (1.0, -1.0),
        )
        for start, end in new_params:
            nodes2 = curve_helpers.specialize_curve(nodes1, start, end)
            result = self._call_function_under_test(nodes1, nodes2)
            if start == 2.0 or end == 2.0:
                expected = (
                    (0.0, -start / (end - start)),
                    (1.0, (end - 1.0) / (end - start)),
                )
            elif start == -1.0 or end == -1.0:
                expected = (
                    (0.0, -start / (end - start)),
                    (1.0, (1.0 - start) / (end - start)),
                )
            else:
                expected = ((start, 0.0), (end, 1.0))
            self.assertEqual(result, expected)

    def test_fully_contained(self):
        from bezier.hazmat import curve_helpers

        nodes1 = np.asfortranarray([[-1.0, 0.0, 2.0], [-1.0, 2.0, 0.0]])
        new_params = ((0.25, 0.75), (0.75, 0.25), (-0.5, 1.5), (1.5, -0.5))
        for start, end in new_params:
            nodes2 = curve_helpers.specialize_curve(nodes1, start, end)
            result = self._call_function_under_test(nodes1, nodes2)
            if start == -0.5 or end == -0.5:
                expected = (
                    (0.0, (end - 1.0) / (end - start)),
                    (1.0, end / (end - start)),
                )
            else:
                expected = ((start, 0.0), (end, 1.0))
            self.assertEqual(result, expected)

    def test_staggered_overlap(self):
        from bezier.hazmat import curve_helpers

        nodes1 = np.asfortranarray(
            [[0.0, 1.0, 1.0, 3.0], [-1.0, 2.0, 0.0, 2.0]]
        )
        new_params = ((0.5, 1.5), (1.5, 0.5), (-0.5, 0.5), (0.5, -0.5))
        for start, end in new_params:
            nodes2 = curve_helpers.specialize_curve(nodes1, start, end)
            result = self._call_function_under_test(nodes1, nodes2)
            if start == 1.5 or end == 1.5:
                expected = ((0.5, start - 0.5), (1.0, 0.5))
            else:
                expected = ((0.0, 0.5), (0.5, end + 0.5))
            self.assertEqual(result, expected)

    def test_one_endpoint_from_each_not_coincident(self):
        nodes1 = np.asfortranarray([[0.0, 16.0, 16.0], [0.0, 0.0, 16.0]])
        nodes2 = np.asfortranarray([[7.0, 7.0, 23.0], [1.0, 17.0, 17.0]])
        result = self._call_function_under_test(nodes1, nodes2)
        self.assertIsNone(result)

    def test_both_endpoints_on_other_not_coincident(self):
        nodes1 = np.asfortranarray(
            [[0.0, 32.0, 96.0, 128.0], [0.0, 32.0, 32.0, 0.0]]
        )
        nodes2 = np.asfortranarray(
            [[29.0, 49.0, 79.0, 99.0], [18.0, 38.0, 38.0, 18.0]]
        )
        result = self._call_function_under_test(nodes1, nodes2)
        self.assertIsNone(result)
        result = self._call_function_under_test(nodes2, nodes1)
        self.assertIsNone(result)

    def test_endpoint_in_middle_not_coincident(self):
        nodes1 = np.asfortranarray([[0.0, 2.0, 4.0], [0.0, 2.0, 0.0]])
        nodes2 = np.asfortranarray([[2.0, 4.0, 6.0], [1.0, 3.0, 2.0]])
        result = self._call_function_under_test(nodes1, nodes2)
        self.assertIsNone(result)
        result = self._call_function_under_test(nodes2, nodes1)
        self.assertIsNone(result)


class Test_check_lines(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(first, second):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.check_lines(first, second)

    def test_not_linearized(self):
        nodes1 = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 0.0]])
        curve1 = subdivided_curve(nodes1)
        nodes2 = np.asfortranarray([[0.0, 1.0, 2.0], [1.0, 0.0, 1.0]])
        curve2 = subdivided_curve(nodes2)
        both_linear, result = self._call_function_under_test(curve1, curve2)
        self.assertFalse(both_linear)
        self.assertIsNone(result)

    def test_nonzero_linearization_error(self):
        nodes1 = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 0.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1)
        self.assertGreater(lin1.error, 0.0)
        nodes2 = np.asfortranarray([[0.0, 1.0, 2.0], [1.0, 0.0, 1.0]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2)
        self.assertGreater(lin2.error, 0.0)
        both_linear, result = self._call_function_under_test(lin1, lin2)
        self.assertFalse(both_linear)
        self.assertIsNone(result)

    def test_do_intersect(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.0)
        nodes2 = np.asfortranarray([[0.0, 1.0], [1.0, 0.0]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        both_linear, result = self._call_function_under_test(lin1, lin2)
        self.assertTrue(both_linear)
        parameters, coincident = result
        expected = np.asfortranarray([[0.5], [0.5]])
        self.assertEqual(parameters, expected)
        self.assertFalse(coincident)

    def test_do_not_intersect(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.0)
        nodes2 = np.asfortranarray([[2.0, 3.0], [3.0, 2.0]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        both_linear, result = self._call_function_under_test(lin1, lin2)
        self.assertTrue(both_linear)
        parameters, coincident = result
        self.assertEqual(parameters.shape, (2, 0))
        self.assertFalse(coincident)

    def test_parallel_disjoint(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.0)
        nodes2 = np.asfortranarray([[0.0, 1.0], [1.0, 2.0]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        both_linear, result = self._call_function_under_test(lin1, lin2)
        self.assertTrue(both_linear)
        parameters, coincident = result
        self.assertEqual(parameters.shape, (2, 0))
        self.assertFalse(coincident)

    def test_parallel_coincident(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = subdivided_curve(nodes1)
        lin1 = make_linearization(curve1, error=0.0)
        nodes2 = np.asfortranarray([[0.5, 4.5], [0.5, 4.5]])
        curve2 = subdivided_curve(nodes2)
        lin2 = make_linearization(curve2, error=0.0)
        both_linear, result = self._call_function_under_test(lin1, lin2)
        self.assertTrue(both_linear)
        parameters, coincident = result
        expected = np.asfortranarray([[0.5, 1.0], [0.0, 0.125]])
        self.assertEqual(parameters, expected)
        self.assertTrue(coincident)


class Test_all_intersections(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(nodes_first, nodes_second, **kwargs):
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.all_intersections(
            nodes_first, nodes_second, **kwargs
        )

    def test_no_intersections(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        nodes2 = np.asfortranarray([[3.0, 4.0], [3.0, 3.0]])
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        self.assertEqual(intersections.shape, (2, 0))
        self.assertFalse(coincident)

    def test_quadratics_intersect_once(self):
        # NOTE: ``nodes1`` is a specialization of [0, 0], [1/2, 1], [1, 1]
        #       onto the interval [1/4, 1] and ``nodes`` is a specialization
        #       of [0, 1], [1/2, 1], [1, 0] onto the interval [0, 3/4].
        #       We expect them to intersect at s = 1/3, t = 2/3, which is
        #       the point [1/2, 3/4].
        nodes1 = np.asfortranarray([[0.25, 0.625, 1.0], [0.4375, 1.0, 1.0]])
        nodes2 = np.asfortranarray([[0.0, 0.375, 0.75], [1.0, 1.0, 0.4375]])
        s_val = 1.0 / 3.0
        t_val = 2.0 / 3.0
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        # Due to round-off, the answer may be wrong by a tiny wiggle.
        self.assertEqual(intersections.shape, (2, 1))
        utils.almost(self, intersections[0, 0], s_val, 1)
        self.assertEqual(intersections[1, 0], t_val)
        self.assertFalse(coincident)

    def test_tangent_curves_parallel_when_linearized(self):
        nodes1 = np.asfortranarray([[0.0, 0.375, 0.75], [0.0, 0.75, 0.375]])
        nodes2 = np.asfortranarray([[0.25, 0.625, 1.0], [0.625, 0.25, 1.0]])
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        # Due to round-off, the answer may be wrong by a tiny wiggle.
        self.assertEqual(intersections.shape, (2, 1))
        utils.almost(self, 2.0 / 3.0, intersections[0, 0], 1)
        utils.almost(self, 1.0 / 3.0, intersections[1, 0], 1)
        self.assertFalse(coincident)

    def test_pruned_candidates(self):
        nodes1 = np.asfortranarray([[0.0, -0.5, 1.0], [0.0, 1.5, 1.0]])
        nodes2 = np.asfortranarray([[-1.0, 0.5, 0.0], [1.0, 0.5, 2.0]])
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        expected = np.asfortranarray([[0.5], [0.5]])
        self.assertEqual(intersections, expected)
        self.assertFalse(coincident)

    def test_pruned_candidates_odd_index(self):
        # Same as ``test_pruned_candidates`` but re-parameterized on [0, 2]
        # and [-1, 1], respectively.
        nodes1 = np.asfortranarray([[0.0, -1.0, 6.0], [0.0, 3.0, -2.0]])
        nodes2 = np.asfortranarray([[-6.0, 1.0, 0.0], [4.0, -1.0, 2.0]])
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        expected = np.asfortranarray([[0.25], [0.75]])
        self.assertEqual(intersections, expected)
        self.assertFalse(coincident)

    def test_non_convergence(self):
        from bezier import _py_geometric_intersection

        multiplier = 16384.0
        nodes1 = multiplier * np.asfortranarray(
            [[0.0, 4.5, 9.0], [0.0, 9.0, 0.0]]
        )
        nodes2 = multiplier * np.asfortranarray([[0.0, 6.0], [8.0, 0.0]])
        with self.assertRaises(ValueError) as exc_info:
            self._call_function_under_test(nodes1, nodes2)
        exc_args = exc_info.exception.args
        expected = _py_geometric_intersection._NO_CONVERGE_TEMPLATE.format(
            _py_geometric_intersection._MAX_INTERSECT_SUBDIVISIONS
        )
        self.assertEqual(exc_args, (expected,))

    def test_duplicates(self):
        # After three subdivisions, there are 8 pairs of curve segments
        # which have bounding boxes that touch at corners (these corners are
        # also intersections). This test makes sure the duplicates are
        # de-duplicated.
        nodes1 = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 1.0, 0.0]])
        nodes2 = np.asfortranarray([[0.0, 0.5, 1.0], [0.75, -0.25, 0.75]])
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        expected = np.asfortranarray([[0.25, 0.75], [0.25, 0.75]])
        self.assertEqual(intersections, expected)
        self.assertFalse(coincident)

    def test_hull_failure(self):
        nodes1 = np.asfortranarray(
            [
                [
                    -0.7838204403623438,
                    -0.7894577677825452,
                    -0.7946421067207265,
                    -0.799367666650849,
                ],
                [
                    -0.25519640597397464,
                    -0.24259531488131633,
                    -0.22976394420044136,
                    -0.21671303774854855,
                ],
            ]
        )
        nodes2 = np.asfortranarray(
            [
                [
                    -0.7993236103108717,
                    -0.8072986524226636,
                    -0.8152736945344552,
                    -0.8232487366462472,
                ],
                [
                    -0.21683567278362156,
                    -0.21898490744674426,
                    -0.2211341421098668,
                    -0.2232833767729893,
                ],
            ]
        )
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        self.assertEqual(intersections.shape, (2, 0))
        self.assertFalse(coincident)
        intersections, coincident = self._call_function_under_test(
            nodes2, nodes1
        )
        self.assertEqual(intersections.shape, (2, 0))
        self.assertFalse(coincident)

    def test_coincident(self):
        nodes1 = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 0.25, 0.0]])
        nodes2 = np.asfortranarray([[0.5, 0.75, 1.0], [0.125, 0.125, 0.0]])
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        expected = np.asfortranarray([[0.5, 1.0], [0.0, 1.0]])
        self.assertEqual(intersections, expected)
        self.assertTrue(coincident)

    def test_triple_root(self):
        from bezier import _py_intersection_helpers

        # Curves intersect and are tangent with the same curvature.
        nodes1 = np.asfortranarray([[12.0, -4.0, -4.0], [4.0, -4.0, 4.0]])
        nodes2 = np.asfortranarray([[6.0, -2.0, -2.0], [1.0, -1.0, 1.0]])
        with self.assertRaises(NotImplementedError) as exc_info:
            self._call_function_under_test(nodes1, nodes2)
        expected = (_py_intersection_helpers.NEWTON_NO_CONVERGE,)
        self.assertEqual(exc_info.exception.args, expected)

    def test_too_many_candidates(self):
        from bezier import _py_geometric_intersection

        nodes1 = np.asfortranarray([[12.0, -12.0, 0.0], [4.0, -8.0, 16.0]])
        nodes2 = np.asfortranarray([[6.0, -6.0, 0.0], [1.0, -2.0, 4.0]])
        with self.assertRaises(NotImplementedError) as exc_info:
            self._call_function_under_test(nodes1, nodes2)
        expected = (_py_geometric_intersection._TOO_MANY_TEMPLATE.format(74),)
        self.assertEqual(exc_info.exception.args, expected)


class TestSubdividedCurve(utils.NumPyTestCase):
    @staticmethod
    def _get_target_class():
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.SubdividedCurve

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor_defaults(self):
        curve = self._make_one(
            unittest.mock.sentinel.nodes, unittest.mock.sentinel.original_nodes
        )
        self.assertIs(curve.nodes, unittest.mock.sentinel.nodes)
        self.assertIs(
            curve.original_nodes, unittest.mock.sentinel.original_nodes
        )
        self.assertEqual(curve.start, 0.0)
        self.assertEqual(curve.end, 1.0)

    def test_constructor_explicit(self):
        curve = self._make_one(
            unittest.mock.sentinel.nodes,
            unittest.mock.sentinel.original_nodes,
            start=3.0,
            end=5.0,
        )
        self.assertIs(curve.nodes, unittest.mock.sentinel.nodes)
        self.assertIs(
            curve.original_nodes, unittest.mock.sentinel.original_nodes
        )
        self.assertEqual(curve.start, 3.0)
        self.assertEqual(curve.end, 5.0)

    def test___dict___property(self):
        curve = self._make_one(
            unittest.mock.sentinel.nodes,
            unittest.mock.sentinel.original_nodes,
            start=0.25,
            end=0.75,
        )
        props_dict = curve.__dict__
        expected = {
            "nodes": unittest.mock.sentinel.nodes,
            "original_nodes": unittest.mock.sentinel.original_nodes,
            "start": 0.25,
            "end": 0.75,
        }
        self.assertEqual(props_dict, expected)
        # Check that modifying ``props_dict`` won't modify ``curve``.
        props_dict["start"] = -1.0
        self.assertNotEqual(curve.start, props_dict["start"])

    def test_subdivide(self):
        klass = self._get_target_class()
        nodes = np.asfortranarray([[0.0, 2.0], [2.0, 0.0]])
        curve = self._make_one(nodes, nodes)
        left, right = curve.subdivide()
        self.assertIsInstance(left, klass)
        self.assertIs(left.original_nodes, nodes)
        self.assertEqual(left.start, 0.0)
        self.assertEqual(left.end, 0.5)
        expected = np.asfortranarray([[0.0, 1.0], [2.0, 1.0]])
        self.assertEqual(left.nodes, expected)
        self.assertIsInstance(right, klass)
        self.assertIs(right.original_nodes, nodes)
        self.assertEqual(right.start, 0.5)
        self.assertEqual(right.end, 1.0)
        expected = np.asfortranarray([[1.0, 2.0], [1.0, 0.0]])
        self.assertEqual(right.nodes, expected)


class TestLinearization(utils.NumPyTestCase):
    NODES = np.asfortranarray([[0.0, 1.0, 5.0], [0.0, 1.0, 6.0]])

    @staticmethod
    def _get_target_class():
        from bezier import _py_geometric_intersection

        return _py_geometric_intersection.Linearization

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def _simple_curve(self):
        return subdivided_curve(self.NODES)

    def test_constructor(self):
        nodes = np.asfortranarray([[4.0, 0.0], [-5.0, 7.0]])
        curve = subdivided_curve(nodes)
        error = 0.125
        linearization = self._make_one(curve, error)
        self.assertIs(linearization.curve, curve)
        self.assertEqual(linearization.error, error)
        self.assertEqual(linearization.start_node, nodes[:, 0])
        self.assertEqual(linearization.end_node, nodes[:, 1])

    def test___dict___property(self):
        nodes = np.asfortranarray([[0.0, 0.0], [1.0, 2.0]])
        curve = subdivided_curve(nodes)
        error = 0.0
        linearization = self._make_one(curve, error)
        props_dict = linearization.__dict__
        # NOTE: We cannot use dictionary equality check because of
        #       the comparison of NumPy arrays.
        self.assertEqual(len(props_dict), 4)
        self.assertIs(props_dict["curve"], curve)
        self.assertEqual(props_dict["error"], error)
        self.assertEqual(props_dict["start_node"], nodes[:, 0])
        self.assertEqual(props_dict["end_node"], nodes[:, 1])
        # Check that modifying ``props_dict`` won't modify ``linearization``.
        props_dict["error"] = 0.5
        self.assertNotEqual(linearization.error, props_dict["error"])

    def test_subdivide(self):
        linearization = self._make_one(self._simple_curve(), np.nan)
        self.assertEqual(linearization.subdivide(), (linearization,))

    def test_start_node_attr(self):
        curve = self._simple_curve()
        linearization = self._make_one(curve, np.nan)
        expected = self.NODES[:, 0]
        self.assertEqual(linearization.start_node, expected)
        # Make sure the data is not a copy (here is a view of the nodes).
        self.assertIs(linearization.start_node.base, curve.nodes)

    def test_end_node_attr(self):
        curve = self._simple_curve()
        linearization = self._make_one(curve, np.nan)
        expected = self.NODES[:, 2]
        self.assertEqual(linearization.end_node, expected)
        # Make sure the data is not a copy (here is a view of the nodes).
        self.assertIs(linearization.end_node.base, curve.nodes)

    def test_from_shape_factory_not_close_enough(self):
        curve = self._simple_curve()
        klass = self._get_target_class()
        new_shape = klass.from_shape(curve)
        self.assertIs(new_shape, curve)

    def test_from_shape_factory_close_enough(self):
        scale_factor = 0.5 ** 27
        nodes = self.NODES * scale_factor
        curve = subdivided_curve(nodes)
        klass = self._get_target_class()
        new_shape = klass.from_shape(curve)
        self.assertIsInstance(new_shape, klass)
        self.assertIs(new_shape.curve, curve)
        # NODES has constant second derivative equal to 2 * [3.0, 4.0].
        expected_error = 0.125 * 2 * 1 * 5.0 * scale_factor
        self.assertEqual(new_shape.error, expected_error)

    def test_from_shape_factory_no_error(self):
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve = subdivided_curve(nodes)
        klass = self._get_target_class()
        new_shape = klass.from_shape(curve)
        self.assertIsInstance(new_shape, klass)
        self.assertIs(new_shape.curve, curve)
        # ``nodes`` is linear, so error is 0.0.
        self.assertEqual(new_shape.error, 0.0)

    def test_from_shape_factory_already_linearized(self):
        error = 0.078125
        linearization = self._make_one(self._simple_curve(), error)
        klass = self._get_target_class()
        new_shape = klass.from_shape(linearization)
        self.assertIs(new_shape, linearization)
        self.assertEqual(new_shape.error, error)


def subdivided_curve(nodes, **kwargs):
    from bezier import _py_geometric_intersection

    if "original_nodes" not in kwargs:
        kwargs["original_nodes"] = nodes
    return _py_geometric_intersection.SubdividedCurve(nodes, **kwargs)


def make_linearization(curve, error=None):
    from bezier import _py_geometric_intersection

    if error is None:
        error = _py_geometric_intersection.linearization_error(curve.nodes)
    return _py_geometric_intersection.Linearization(curve, error)
