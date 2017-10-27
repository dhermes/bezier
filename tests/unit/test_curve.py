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

import unittest.mock

import numpy as np

from tests.unit import utils


class TestCurve(utils.NumPyTestCase):

    ZEROS = np.zeros((2, 2), order='F')

    @staticmethod
    def _get_target_class():
        from bezier import curve

        return curve.Curve

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.625, 0.5],
            [1.0, 0.75],
        ])
        curve = self._make_one(nodes, 2, _copy=False)
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
        nodes = np.asfortranarray([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes, None)

        nodes = np.zeros((2, 2, 2), order='F')
        with self.assertRaises(ValueError):
            self._make_one(nodes, None)

    def test_from_nodes_factory(self):
        nodes = np.asfortranarray([
            [1.0, 2.0, 0.0],
            [1.0, 3.0, 0.0],
        ])
        klass = self._get_target_class()

        curve = klass.from_nodes(
            nodes, start=0.25, end=0.75, root=unittest.mock.sentinel.root)
        self.assertIsInstance(curve, klass)
        self.assertEqual(curve._degree, 1)
        self.assertEqual(curve._dimension, 3)
        self.assertEqual(curve._nodes, nodes)
        self.assertIsNone(curve._length)
        self.assertIsNone(curve._edge_index)
        self.assertIsNone(curve._next_edge)
        self.assertIsNone(curve._previous_edge)
        self.assertEqual(curve._start, 0.25)
        self.assertEqual(curve._end, 0.75)
        self.assertIs(curve._root, unittest.mock.sentinel.root)

    def test___repr__(self):
        degree = 4
        dimension = 3
        nodes = np.zeros((degree + 1, dimension), order='F')
        curve = self._make_one(nodes, degree)
        expected = '<Curve (degree={:d}, dimension={:d})>'.format(
            degree, dimension)
        self.assertEqual(repr(curve), expected)

    def test___repr__custom_endpoints(self):
        from bezier import curve as curve_mod

        degree = 4
        dimension = 3
        start = 0.25
        end = 0.75
        nodes = np.zeros((degree + 1, dimension), order='F')
        curve = self._make_one(nodes, degree, start=start, end=end)
        expected = curve_mod._REPR_TEMPLATE.format(
            'Curve', degree, dimension, start, end)
        self.assertEqual(repr(curve), expected)

    def test__get_degree(self):
        klass = self._get_target_class()

        self.assertEqual(0, klass._get_degree(1))
        self.assertEqual(1, klass._get_degree(2))

    def test_length_property_not_cached(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        curve = self._make_one(nodes, 1)
        self.assertIsNone(curve._length)

        patch = unittest.mock.patch(
            'bezier._curve_helpers.compute_length',
            return_value=unittest.mock.sentinel.length)
        with patch as mocked:
            self.assertEqual(curve.length, unittest.mock.sentinel.length)

            self.assertEqual(mocked.call_count, 1)
            call = mocked.mock_calls[0]
            _, positional, keyword = call
            self.assertEqual(keyword, {})
            self.assertEqual(len(positional), 1)
            self.assertEqual(positional[0], nodes)

    def test_length_property(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        curve = self._make_one(nodes, 1)
        length = 2.817281728
        curve._length = length
        self.assertEqual(curve.length, length)

    def test_start_property(self):
        curve = self._make_one(self.ZEROS, 1,
                               start=unittest.mock.sentinel.start)
        self.assertIs(curve.start, unittest.mock.sentinel.start)

    def test_end_property(self):
        curve = self._make_one(self.ZEROS, 1,
                               end=unittest.mock.sentinel.end)
        self.assertIs(curve.end, unittest.mock.sentinel.end)

    def test_root_property(self):
        curve = self._make_one(self.ZEROS, 1,
                               root=unittest.mock.sentinel.root)
        self.assertIs(curve.root, unittest.mock.sentinel.root)

    def test_edge_index_property(self):
        curve = self._make_one(self.ZEROS, 1)
        curve._edge_index = 7
        self.assertEqual(curve.edge_index, 7)

    def test_next_edge_property(self):
        curve = self._make_one(self.ZEROS, 1)
        curve._next_edge = unittest.mock.sentinel.next_edge
        self.assertIs(curve.next_edge, unittest.mock.sentinel.next_edge)

    def test_previous_edge_property(self):
        curve = self._make_one(self.ZEROS, 1)
        curve._previous_edge = unittest.mock.sentinel.previous
        self.assertIs(curve.previous_edge, unittest.mock.sentinel.previous)

    def test___dict___property(self):
        curve = self._make_one(self.ZEROS, 1, _copy=False)
        props_dict = curve.__dict__
        expected = {
            '_nodes': self.ZEROS,
            '_dimension': 2,
            '_degree': 1,
            '_start': 0.0,
            '_end': 1.0,
            '_root': curve,
            '_length': None,
            '_edge_index': None,
            '_next_edge': None,
            '_previous_edge': None,
        }
        self.assertEqual(props_dict, expected)
        # Check that modifying ``props_dict`` won't modify ``curve``.
        expected['_end'] = 1.5
        self.assertNotEqual(curve._end, expected['_end'])

    def _copy_helper(self, **kwargs):
        np_shape = (2, 2)
        length = kwargs.pop('_length', None)
        curve = self._make_one(np.zeros(np_shape, order='F'), 1, **kwargs)
        curve._length = length
        fake_nodes = unittest.mock.Mock(
            ndim=2, shape=np_shape, spec=['ndim', 'shape', 'copy'])
        curve._nodes = fake_nodes

        copied_nodes = np.zeros(np_shape, order='F')
        fake_nodes.copy.return_value = copied_nodes

        new_curve = curve._copy()
        self.assertIsInstance(new_curve, self._get_target_class())
        self.assertEqual(new_curve._degree, curve._degree)
        self.assertEqual(new_curve._dimension, curve._dimension)
        self.assertIs(new_curve._nodes, copied_nodes)
        self.assertEqual(new_curve._start, curve._start)
        self.assertEqual(new_curve._end, curve._end)
        if 'root' in kwargs:
            self.assertIs(new_curve._root, kwargs['root'])
        else:
            self.assertIs(new_curve._root, new_curve)
        self.assertEqual(new_curve._length, curve._length)
        self.assertIsNone(new_curve._edge_index)
        self.assertIsNone(new_curve._next_edge)
        self.assertIsNone(new_curve._previous_edge)

        fake_nodes.copy.assert_called_once_with(order='F')

    def test__copy(self):
        self._copy_helper()

    def test__copy_with_root(self):
        self._copy_helper(
            start=0.25, end=0.5,
            root=unittest.mock.sentinel.root, _length=2.25)

    def test_evaluate(self):
        s = 0.25
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 0.5],
            [1.0, 1.25],
        ])
        curve = self._make_one(nodes, 2)
        expected = np.asfortranarray([[0.25, 0.265625]])
        result = curve.evaluate(s)
        self.assertEqual(expected, result)

    def test_evaluate_multi(self):
        s_vals = np.asfortranarray([0.0, 0.25, 0.5, 1.0, 1.25])
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.375, 0.375],
            [1.0, 1.0],
        ])
        curve = self._make_one(nodes, 2)
        expected = np.asfortranarray([
            [0.0, 0.0],
            [0.203125, 0.203125],
            [0.4375, 0.4375],
            [1.0, 1.0],
            [1.328125, 1.328125],
        ])
        result = curve.evaluate_multi(s_vals)
        self.assertEqual(expected, result)

    def test_plot_wrong_dimension(self):
        nodes = np.asfortranarray([
            [0.0, 0.0, 0.0],
            [1.0, 3.0, 4.0],
        ])
        curve = self._make_one(nodes, 1)
        with self.assertRaises(NotImplementedError):
            curve.plot(32)

    @unittest.mock.patch('bezier._plot_helpers.new_axis')
    def test_plot_defaults(self, new_axis_mock):
        ax = unittest.mock.Mock(spec=['plot'])
        new_axis_mock.return_value = ax

        nodes = np.asfortranarray([
            [0.0, 1.0],
            [1.0, 3.0],
        ])
        curve = self._make_one(nodes, 1, _copy=False)

        num_pts = 2  # This value is crucial for the plot call.
        result = curve.plot(num_pts)
        self.assertIs(result, ax)

        # Verify mocks.
        new_axis_mock.assert_called_once_with()
        # Check the call to ax.plot(). We can't assert_any_call()
        # since == breaks on NumPy arrays.
        self.assertEqual(ax.plot.call_count, 1)
        call = ax.plot.mock_calls[0]
        utils.check_plot_call(self, call, nodes, color=None, alpha=None)

    @unittest.mock.patch('bezier._plot_helpers.new_axis')
    def test_plot_explicit(self, new_axis_mock):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve = self._make_one(nodes, 1, _copy=False)

        num_pts = 2  # This value is crucial for the plot call.
        ax = unittest.mock.Mock(spec=['plot'])
        color = (0.75, 1.0, 1.0)
        alpha = 0.625
        result = curve.plot(num_pts, color=color, alpha=alpha, ax=ax)
        self.assertIs(result, ax)

        # Verify mocks.
        new_axis_mock.assert_not_called()
        # Check the call to ax.plot(). We can't assert_any_call()
        # since == breaks on NumPy arrays.
        self.assertEqual(ax.plot.call_count, 1)
        call = ax.plot.mock_calls[0]
        utils.check_plot_call(self, call, nodes, color=color, alpha=alpha)

    def test_subdivide(self):
        nodes = np.asfortranarray([
            [0.0, 1.0],
            [4.0, 6.0],
        ])
        klass = self._get_target_class()
        curve = klass.from_nodes(nodes, start=0.25, end=0.75)

        # Call ``subdivide()`` and then compare.
        left, right = curve.subdivide()

        # Check the "left" sub-curve.
        self.assertIs(left._root, curve)
        self.assertEqual(left._degree, 1)
        self.assertEqual(left._start, 0.25)
        self.assertEqual(left._end, 0.5)
        self.assertIsInstance(left, klass)
        expected_l = np.asfortranarray([
            [0.0, 1.0],
            [2.0, 3.5],
        ])
        self.assertEqual(left._nodes, expected_l)

        # Check the "right" sub-curve.
        self.assertIs(right._root, curve)
        self.assertEqual(right._start, 0.5)
        self.assertEqual(right._end, 0.75)
        self.assertIsInstance(right, klass)
        expected_r = np.asfortranarray([
            [2.0, 3.5],
            [4.0, 6.0],
        ])
        self.assertEqual(right._nodes, expected_r)

    def test_subdivide_multilevel_root(self):
        curve = self._make_one(self.ZEROS, 1)
        left, right = curve.subdivide()
        self.assertIs(left.root, curve)
        self.assertIs(right.root, curve)

        one, two = left.subdivide()
        three, four = right.subdivide()
        self.assertIs(one.root, curve)
        self.assertIs(two.root, curve)
        self.assertIs(three.root, curve)
        self.assertIs(four.root, curve)

    def test_intersect_bad_strategy(self):
        curve = self._make_one(self.ZEROS, 1)
        strategy = unittest.mock.sentinel.bad_strategy
        with self.assertRaises(ValueError) as exc_info:
            curve.intersect(curve, strategy=strategy)

        exc_args = exc_info.exception.args
        self.assertEqual(exc_args, ('Unexpected strategy.', strategy))

    def test_intersect_algebraic(self):
        from bezier import _intersection_helpers

        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve1 = self._make_one(nodes1, 1)
        nodes2 = np.asfortranarray([
            [0.0, 1.0],
            [1.0, 0.0],
        ])
        curve2 = self._make_one(nodes2, 1)
        strategy = _intersection_helpers.IntersectionStrategy.ALGEBRAIC

        intersections = curve1.intersect(curve2, strategy=strategy)
        expected = np.asfortranarray([
            [0.5, 0.5],
        ])
        self.assertEqual(intersections, expected)

    def test_intersect_empty(self):
        nodes1 = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve1 = self._make_one(nodes1, 1)

        nodes2 = np.asfortranarray([
            [3.0, 0.0],
            [2.0, 1.0],
        ])
        curve2 = self._make_one(nodes2, 1)

        result = curve1.intersect(curve2)
        self.assertEqual(result.shape, (0, 2))

    def test_intersect_at_boundary(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, -0.25],
            [1.0, 0.0],
        ])
        curve = self._make_one(nodes, 2)
        left, right = curve.subdivide()
        # Patch the ``root``-s because intersection assumes (correctly or
        # incorrectly) that ``root is self``.
        left._root = left
        right._root = right

        result = left.intersect(right)
        expected = np.asfortranarray([[0.5, -0.125]])
        self.assertEqual(result, expected)

    def _intersect_helper(self, **kwargs):
        # NOTE: ``nodes1`` is a specialization of [0, 0], [1/2, 1], [1, 1]
        #       onto the interval [1/4, 1] and ``nodes`` is a specialization
        #       of [0, 1], [1/2, 1], [1, 0] onto the interval [0, 3/4].
        #       We expect them to intersect at s = 1/3, t = 2/3, which is
        #       the point [1/2, 3/4].
        nodes_left = np.asfortranarray([
            [0.25, 0.4375],
            [0.625, 1.0],
            [1.0, 1.0],
        ])
        left = self._make_one(nodes_left, 2)

        nodes_right = np.asfortranarray([
            [0.0, 1.0],
            [0.375, 1.0],
            [0.75, 0.4375],
        ])
        right = self._make_one(nodes_right, 2)

        result = left.intersect(right, **kwargs)
        expected = np.asfortranarray([[0.5, 0.75]])
        self.assertEqual(result, expected)

    def test_intersect(self):
        self._intersect_helper()

    def test_intersect_no_verify(self):
        self._intersect_helper(_verify=False)

    def test_intersect_non_curve(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, -0.25],
            [1.0, 0.0],
        ])
        curve = self._make_one(nodes, 2)
        with self.assertRaises(TypeError):
            curve.intersect(object())

    def test_intersect_unsupported_dimension(self):
        nodes = np.asfortranarray([
            [0.0, 0.0, 0.0],
            [0.5, -0.25, 0.75],
            [1.0, 0.0, 1.25],
        ])
        curve1 = self._make_one(nodes, 2)
        curve2 = self._make_one(nodes[:, :2], 2)

        with self.assertRaises(NotImplementedError):
            curve1.intersect(curve2)
        with self.assertRaises(NotImplementedError):
            curve2.intersect(curve1)

    def test_elevate(self):
        nodes = np.asfortranarray([
            [0.0, 0.5],
            [1.0, 1.0],
            [3.0, 2.0],
            [3.5, 4.0],
        ])
        curve = self._make_one(nodes, 3)
        self.assertEqual(curve.degree, 3)
        elevated = curve.elevate()
        self.assertEqual(elevated.degree, 4)
        self.assertIs(elevated.root, elevated)
        self.assertEqual(elevated.start, curve.start)
        self.assertEqual(elevated.end, curve.end)

        s_vals = np.linspace(0.0, 1.0, 64 + 1)
        orig_vals = curve.evaluate_multi(s_vals)
        new_vals = elevated.evaluate_multi(s_vals)
        self.assertEqual(orig_vals, new_vals)

    def test_reduce_(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 3.0],
            [2.0, 3.0],
            [3.0, 0.0],
        ])
        curve = self._make_one(nodes, 3)
        self.assertEqual(curve.degree, 3)
        reduced = curve.reduce_()

        expected = np.asfortranarray([
            [0.0, 0.0],
            [1.5, 4.5],
            [3.0, 0.0],
        ])
        self.assertEqual(reduced.nodes, expected)
        self.assertIs(reduced.root, reduced)
        self.assertEqual(reduced.start, curve.start)
        self.assertEqual(reduced.end, curve.end)

        s_vals = np.linspace(0.0, 1.0, 64 + 1)
        orig_vals = curve.evaluate_multi(s_vals)
        new_vals = reduced.evaluate_multi(s_vals)
        self.assertEqual(orig_vals, new_vals)

    def test_specialize(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 6.0],
            [5.0, 2.0],
        ])
        curve = self._make_one(nodes, 2)
        start = 0.25
        end = 0.875
        new_curve = curve.specialize(start, end)

        self.assertEqual(new_curve.start, start)
        self.assertEqual(new_curve.end, end)
        self.assertIs(new_curve.root, curve)
        expected = np.asfortranarray([
            [0.6875, 2.375],
            [1.78125, 4.5625],
            [4.046875, 2.84375],
        ])
        self.assertEqual(new_curve.nodes, expected)

    def test_locate_wrong_shape(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
        ])
        curve = self._make_one(nodes, 1)
        point = np.asfortranarray([[0.0, 1.0, 2.0]])
        with self.assertRaises(ValueError):
            curve.locate(point)

    def test_locate(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, -1.0],
            [5.0, 1.0],
        ])
        curve = self._make_one(nodes, 3)
        s_val = 0.75
        point = curve.evaluate_multi(np.asfortranarray([s_val]))
        result = curve.locate(point)
        self.assertEqual(result, s_val)

    def test_locate_non_default_endpoints(self):
        nodes = np.asfortranarray([
            [0.0, 0.0],
            [0.5, 1.0],
            [1.0, 0.0],
        ])
        curve = self._make_one(
            nodes, degree=2, start=0.25, end=1.0, _copy=False)
        # C(1/2) = p
        point = np.asfortranarray([[0.5, 0.5]])
        result = curve.locate(point)
        self.assertEqual(result, 0.5)
