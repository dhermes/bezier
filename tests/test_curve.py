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

import mock
import numpy as np

from tests import utils


class TestCurve(utils.NumPyTestCase):

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
        self.assertIsNone(curve._edge_index)
        self.assertIsNone(curve._next_edge)
        self.assertIsNone(curve._previous_edge)
        self.assertEqual(curve._start, 0.0)
        self.assertEqual(curve._end, 1.0)
        self.assertIs(curve._root, curve)

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

    def test___repr__custom_endpoints(self):
        from bezier import curve as curve_mod

        degree = 4
        dimension = 3
        start = 0.25
        end = 0.75
        nodes = np.zeros((degree + 1, dimension))
        curve = self._make_one(nodes, start=start, end=end)
        expected = curve_mod._REPR_TEMPLATE.format(
            'Curve', degree, dimension, start, end)
        self.assertEqual(repr(curve), expected)

    def test__get_degree(self):
        klass = self._get_target_class()

        self.assertEqual(0, klass._get_degree(1))
        self.assertEqual(1, klass._get_degree(2))

    def test_length_property_not_cached(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        curve = self._make_one(nodes)
        self.assertIsNone(curve._length)

        patch = mock.patch(
            'bezier._curve_helpers.compute_length',
            return_value=mock.sentinel.length)
        with patch as mocked:
            self.assertEqual(curve.length, mock.sentinel.length)

            self.assertEqual(mocked.call_count, 1)
            call = mocked.mock_calls[0]
            _, positional, keyword = call
            self.assertEqual(keyword, {})
            self.assertEqual(len(positional), 2)
            self.assertEqual(positional[0], nodes)
            self.assertEqual(positional[1], 1)

    def test_length_property(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        curve = self._make_one(nodes)
        length = 2.817281728
        curve._length = length
        self.assertEqual(curve.length, length)

    def test_start_property(self):
        curve = self._make_one(np.zeros((2, 2)),
                               start=mock.sentinel.start)
        self.assertIs(curve.start, mock.sentinel.start)

    def test_end_property(self):
        curve = self._make_one(np.zeros((2, 2)),
                               end=mock.sentinel.end)
        self.assertIs(curve.end, mock.sentinel.end)

    def test_root_property(self):
        curve = self._make_one(np.zeros((2, 2)),
                               root=mock.sentinel.root)
        self.assertIs(curve.root, mock.sentinel.root)

    def test_edge_index_property(self):
        curve = self._make_one(np.zeros((2, 2)))
        curve._edge_index = 7
        self.assertEqual(curve.edge_index, 7)

    def test_next_edge_property(self):
        curve = self._make_one(np.zeros((2, 2)))
        curve._next_edge = mock.sentinel.next_edge
        self.assertIs(curve.next_edge, mock.sentinel.next_edge)

    def test_previous_edge_property(self):
        curve = self._make_one(np.zeros((2, 2)))
        curve._previous_edge = mock.sentinel.previous
        self.assertIs(curve.previous_edge, mock.sentinel.previous)

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
        self.assertEqual(expected, result)

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
        self.assertEqual(expected, result)

    def _check_plot_call(self, call, expected, **kwargs):
        # Unpack the call as name, positional args, keyword args
        _, positional, keyword = call
        self.assertEqual(keyword, kwargs)
        self.assertEqual(len(positional), 2)
        self.assertEqual(positional[0], expected[:, 0])
        self.assertEqual(positional[1], expected[:, 1])

    def _plot_helper(self, show=False):
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

        self.assertIs(result, ax)

        # Check mocks.
        plt.figure.assert_called_once_with()
        figure.gca.assert_called_once_with()
        # Can't use nodes[:, col] since == breaks on array.
        self.assertEqual(ax.plot.call_count, 1)
        call = ax.plot.mock_calls[0]
        self._check_plot_call(call, nodes)
        if show:
            plt.show.assert_called_once_with()
        else:
            plt.show.assert_not_called()

    def test_plot(self):
        self._plot_helper()

    def test_plot_show(self):
        self._plot_helper(show=True)

    def test_plot_existing_axis(self):
        import matplotlib.lines

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve = self._make_one(nodes)
        plt = mock.Mock()

        ax = mock.Mock()
        color = (0.25, 0.25, 0.125)
        line = matplotlib.lines.Line2D([], [], color=color)
        ax.plot.return_value = (line,)

        with mock.patch('bezier.surface.plt', new=plt):
            result = curve.plot(2, ax=ax)

        self.assertIs(result, ax)

        # Check mocks.
        plt.figure.assert_not_called()
        plt.show.assert_not_called()

        # Can't use nodes[:, col] since == breaks on array.
        self.assertEqual(ax.plot.call_count, 1)
        call = ax.plot.mock_calls[0]
        self._check_plot_call(call, nodes)

    def test_plot_wrong_dimension(self):
        nodes = np.array([
            [0.0, 0.0, 0.0],
            [1.0, 3.0, 4.0],
        ])
        curve = self._make_one(nodes)
        with self.assertRaises(NotImplementedError):
            curve.plot(32)

    def test_subdivide_multilevel_root(self):
        curve = self._make_one(np.zeros((2, 2)))
        left, right = curve.subdivide()
        self.assertIs(left.root, curve)
        self.assertIs(right.root, curve)

        one, two = left.subdivide()
        three, four = right.subdivide()
        self.assertIs(one.root, curve)
        self.assertIs(two.root, curve)
        self.assertIs(three.root, curve)
        self.assertIs(four.root, curve)

    def _subdivide_helper(self, nodes, expected_l, expected_r):
        klass = self._get_target_class()

        curve = self._make_one(nodes)
        left, right = curve.subdivide()
        self.assertIs(left.root, curve)
        self.assertIs(right.root, curve)

        self.assertIsInstance(left, klass)
        self.assertEqual(left._nodes, expected_l)
        self.assertIsInstance(right, klass)
        self.assertEqual(right._nodes, expected_r)

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
            self.assertEqual(main_vals, sub_vals)

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

    def test_intersect_empty(self):
        nodes1 = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve1 = self._make_one(nodes1)

        nodes2 = np.array([
            [3.0, 0.0],
            [2.0, 1.0],
        ])
        curve2 = self._make_one(nodes2)

        result = curve1.intersect(curve2)
        self.assertEqual(result.shape, (0, 2))

    def test_intersect_at_boundary(self):
        nodes = np.array([
            [0.0, 0.0],
            [0.5, -0.25],
            [1.0, 0.0],
        ])
        curve = self._make_one(nodes)
        left, right = curve.subdivide()

        result = left.intersect(right)
        expected = np.array([[0.5, -0.125]])
        self.assertEqual(result, expected)

    def test_intersect(self):
        import bezier

        # NOTE: ``nodes1`` is a specialization of [0, 0], [1/2, 1], [1, 1]
        #       onto the interval [1/4, 1] and ``nodes`` is a specialization
        #       of [0, 1], [1/2, 1], [1, 0] onto the interval [0, 3/4].
        #       We expect them to intersect at s = 1/3, t = 2/3, which is
        #       the point [1/2, 3/4].
        nodes_left = np.array([
            [0.25, 0.4375],
            [0.625, 1.0],
            [1.0, 1.0],
        ])
        left = bezier.Curve(nodes_left)

        nodes_right = np.array([
            [0.0, 1.0],
            [0.375, 1.0],
            [0.75, 0.4375],
        ])
        right = bezier.Curve(nodes_right)

        result = left.intersect(right)
        expected = np.array([[0.5, 0.75]])
        self.assertEqual(result, expected)

    def test_intersect_non_curve(self):
        nodes = np.array([
            [0.0, 0.0],
            [0.5, -0.25],
            [1.0, 0.0],
        ])
        curve = self._make_one(nodes)
        with self.assertRaises(TypeError):
            curve.intersect(object())

    def test_intersect_unsupported_dimension(self):
        nodes = np.array([
            [0.0, 0.0, 0.0],
            [0.5, -0.25, 0.75],
            [1.0, 0.0, 1.25],
        ])
        curve1 = self._make_one(nodes)
        curve2 = self._make_one(nodes[:, :2])

        with self.assertRaises(NotImplementedError):
            curve1.intersect(curve2)
        with self.assertRaises(NotImplementedError):
            curve2.intersect(curve1)

    def test_elevate(self):
        nodes = np.array([
            [0.0, 0.5],
            [1.0, 1.0],
            [3.0, 2.0],
            [3.5, 4.0],
        ])
        curve = self._make_one(nodes)
        self.assertEqual(curve.degree, 3)
        elevated = curve.elevate()
        self.assertEqual(elevated.degree, 4)

        s_vals = np.linspace(0.0, 1.0, 64 + 1)
        orig_vals = curve.evaluate_multi(s_vals)
        new_vals = elevated.evaluate_multi(s_vals)
        self.assertEqual(orig_vals, new_vals)

    def test_specialize(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 6.0],
            [5.0, 2.0],
        ])
        curve = self._make_one(nodes)
        start = 0.25
        end = 0.875
        new_curve = curve.specialize(start, end)

        self.assertEqual(new_curve.start, start)
        self.assertEqual(new_curve.end, end)
        self.assertIs(new_curve.root, curve)
        expected = np.array([
            [0.6875, 2.375],
            [1.78125, 4.5625],
            [4.046875, 2.84375],
        ])
        self.assertEqual(new_curve.nodes, expected)

    def test_locate_wrong_shape(self):
        curve = self._make_one(np.array([
            [0.0, 0.0],
            [1.0, 1.0],
        ]))
        point = np.array([[0.0, 1.0, 2.0]])
        with self.assertRaises(ValueError):
            curve.locate(point)

    def test_locate(self):
        curve = self._make_one(np.array([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, -1.0],
            [5.0, 1.0],
        ]))
        s_val = 0.75
        point = curve.evaluate_multi(np.array([s_val]))
        result = curve.locate(point)
        self.assertEqual(result, s_val)
