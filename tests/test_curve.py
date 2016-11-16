# See the License for the specific language governing permissions and
# limitations under the License.
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

from tests import utils


class Test__make_subdivision_matrix(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(degree):
        from bezier import curve

        return curve._make_subdivision_matrix(degree)

    def _helper(self, degree, expected):
        result = self._call_function_under_test(degree)
        self.assertTrue(np.all(result == expected))

    def test_linear(self):
        from bezier import curve

        self._helper(1, curve._LINEAR_SUBDIVIDE)

    def test_quadratic(self):
        from bezier import curve

        self._helper(2, curve._QUADRATIC_SUBDIVIDE)

    def test_cubic(self):
        from bezier import curve

        self._helper(3, curve._CUBIC_SUBDIVIDE)


class TestCurve(unittest.TestCase):

    @staticmethod
    def _get_target_class():
        from bezier import curve

        return curve.Curve

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        nodes = np.array([
            [0.0, 0.0],
            [0.625, 0.5],
            [1.0, 0.75],
        ])
        curve = self._make_one(nodes, _copy=False)
        self.assertEqual(curve._degree, 2)
        self.assertEqual(curve._dimension, 2)
        self.assertIs(curve._nodes, nodes)
        self.assertIsNone(curve._length)

    def test_constructor_wrong_dimension(self):
        nodes = np.array([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

        nodes = np.zeros((2, 2, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test_constructor_bad_degree(self):
        nodes = np.array([
            [1.0, 2.0],
        ])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

        nodes = np.zeros((0, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test___repr__(self):
        degree = 4
        dimension = 3
        nodes = np.zeros((degree + 1, dimension))
        curve = self._make_one(nodes)
        expected = '<Curve (degree={:d}, dimension={:d})>'.format(
            degree, dimension)
        self.assertEqual(repr(curve), expected)

    def test_length_property_not_cached(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        curve = self._make_one(nodes)
        self.assertIsNone(curve._length)
        with self.assertRaises(NotImplementedError):
            getattr(curve, 'length')

    def test_length_property(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        curve = self._make_one(nodes)
        length = 2.817281728
        curve._length = length
        self.assertEqual(curve.length, length)

    def test_evaluate(self):
        s = 0.25
        nodes = np.array([
            [0.0, 0.0],
            [0.5, 0.5],
            [1.0, 1.25],
        ])
        curve = self._make_one(nodes)
        expected = np.array([0.25, 0.265625])
        result = curve.evaluate(s)
        self.assertTrue(np.all(expected == result))

    def test_evaluate_multi(self):
        s_vals = np.array([0.0, 0.25, 0.5, 1.0, 1.25])
        nodes = np.array([
            [0.0, 0.0],
            [0.375, 0.375],
            [1.0, 1.0],
        ])
        curve = self._make_one(nodes)
        expected = np.array([
            [0.0, 0.0],
            [0.203125, 0.203125],
            [0.4375, 0.4375],
            [1.0, 1.0],
            [1.328125, 1.328125],
        ])
        result = curve.evaluate_multi(s_vals)
        self.assertTrue(np.all(expected == result))

    def test_evaluate_multi_calls_evaluate(self):
        import mock

        s_val1 = 3.14159
        s_val2 = 2.817281728
        s_vals = np.array([s_val1, s_val2])
        num_pts = len(s_vals)
        curve = self._make_one(np.zeros((2, 1)))
        ret_vals = [10.0, -1.0]
        curve.evaluate = mock.Mock(side_effect=ret_vals)

        result = curve.evaluate_multi(s_vals)
        self.assertEqual(result.shape, (num_pts, 1))
        # pylint: disable=no-member
        expected = np.array([ret_vals]).T
        # pylint: enable=no-member
        self.assertTrue(np.all(result == expected))

        curve.evaluate.assert_any_call(s_val1)
        curve.evaluate.assert_any_call(s_val2)
        self.assertEqual(curve.evaluate.call_count, 2)

    def _plot_helper(self, show=False):
        import mock

        nodes = np.array([
            [0.0, 1.0],
            [1.0, 3.0],
        ])
        curve = self._make_one(nodes)
        plt = mock.Mock()

        figure = mock.Mock()
        plt.figure.return_value = figure
        ax = mock.Mock()
        figure.gca.return_value = ax

        with mock.patch('bezier.curve.plt', new=plt):
            if show:
                result = curve.plot(2, show=True)
            else:
                result = curve.plot(2)

        self.assertIs(result, figure)

        # Check mocks.
        plt.figure.assert_called_once_with()
        figure.gca.assert_called_once_with()
        # Can't use nodes[:, col] since == breaks on array.
        self.assertEqual(ax.plot.call_count, 1)
        call = ax.plot.mock_calls[0]
        # Unpack the call as name, positional args, keyword args
        _, positional, keyword = call
        self.assertEqual(keyword, {})
        self.assertEqual(len(positional), 2)
        self.assertTrue(np.all(positional[0] == nodes[:, 0]))
        self.assertTrue(np.all(positional[1] == nodes[:, 1]))
        if show:
            plt.show.assert_called_once_with()
        else:
            plt.show.assert_not_called()

    def test_plot(self):
        self._plot_helper()

    def test_plot_show(self):
        self._plot_helper(show=True)

    def test_plot_wrong_dimension(self):
        nodes = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 3.0, 4.0],
        ])
        curve = self._make_one(nodes)
        with self.assertRaises(NotImplementedError):
            curve.plot(32)

    def _subdivide_helper(self, nodes, expected_l, expected_r):
        klass = self._get_target_class()

        curve = self._make_one(nodes)
        left, right = curve.subdivide()

        self.assertIsInstance(left, klass)
        self.assertTrue(np.all(left._nodes == expected_l))
        self.assertIsInstance(right, klass)
        self.assertTrue(np.all(right._nodes == expected_r))

    def _subdivide_points_check(self, curve, pts_exponent=5):
        # Using the exponent means that ds = 1/2**exp, which
        # can be computed without roundoff.
        num_pts = 2**pts_exponent + 1
        left, right = curve.subdivide()

        left_half = np.linspace(0.0, 0.5, num_pts)
        right_half = np.linspace(0.5, 1.0, num_pts)
        unit_interval = np.linspace(0.0, 1.0, num_pts)

        pairs = [
            (left, left_half),
            (right, right_half),
        ]
        for sub_curve, half in pairs:
            # Make sure sub_curve([0, 1]) == curve(half)
            main_vals = curve.evaluate_multi(half)
            sub_vals = sub_curve.evaluate_multi(unit_interval)
            self.assertTrue(np.all(main_vals == sub_vals))

    def test_subdivide_line(self):
        nodes = np.array([
            [0.0, 1.0],
            [4.0, 6.0],
        ])
        expected_l = np.array([
            [0.0, 1.0],
            [2.0, 3.5],
        ])
        expected_r = np.array([
            [2.0, 3.5],
            [4.0, 6.0],
        ])
        self._subdivide_helper(nodes, expected_l, expected_r)

    def test_subdivide_line_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(2, 2), seed=88991, num_bits=8)

        curve = self._make_one(nodes)
        self.assertEqual(curve.degree, 1)
        self._subdivide_points_check(curve)

    def test_subdivide_quadratic(self):
        nodes = np.array([
            [0.0, 1.0],
            [4.0, 6.0],
            [7.0, 3.0],
        ])
        expected_l = np.array([
            [0.0, 1.0],
            [2.0, 3.5],
            [3.75, 4.0],
        ])
        expected_r = np.array([
            [3.75, 4.0],
            [5.5, 4.5],
            [7.0, 3.0],
        ])
        self._subdivide_helper(nodes, expected_l, expected_r)

    def test_subdivide_quadratic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(3, 2), seed=10764, num_bits=8)

        curve = self._make_one(nodes)
        self.assertEqual(curve.degree, 2)
        self._subdivide_points_check(curve)

    def test_subdivide_cubic(self):
        nodes = np.array([
            [0.0, 1.0],
            [4.0, 6.0],
            [7.0, 3.0],
            [6.0, 5.0],
        ])
        expected_l = np.array([
            [0.0, 1.0],
            [2.0, 3.5],
            [3.75, 4.0],
            [4.875, 4.125],
        ])
        expected_r = np.array([
            [4.875, 4.125],
            [6.0, 4.25],
            [6.5, 4.0],
            [6.0, 5.0],
        ])
        self._subdivide_helper(nodes, expected_l, expected_r)

    def test_subdivide_cubic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(4, 2), seed=990077, num_bits=8)

        curve = self._make_one(nodes)
        self.assertEqual(curve.degree, 3)
        self._subdivide_points_check(curve)

    def test_subdivide_dynamic_subdivision_matrix(self):
        degree = 4
        shape = (degree + 1, 2)
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=shape, seed=103, num_bits=8)

        curve = self._make_one(nodes)
        self.assertEqual(curve.degree, degree)
        self._subdivide_points_check(curve)
