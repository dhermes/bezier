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

import numpy as np


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
