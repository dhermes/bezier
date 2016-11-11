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


class TestSurface(unittest.TestCase):

    @staticmethod
    def _get_target_class():
        from bezier import surface

        return surface.Surface

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        nodes = np.array([
            [0.0, 0.0],
            [0.625, 0.5],
            [1.0, 0.75],
        ])
        surface = self._make_one(nodes)
        self.assertEqual(surface._degree, 1)
        self.assertEqual(surface._dimension, 2)
        self.assertIs(surface._nodes, nodes)
        self.assertIsNone(surface._area)
        self.assertIsNone(surface._edges)

    def test_constructor_wrong_dimension(self):
        nodes = np.array([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

        nodes = np.zeros((2, 2, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test_constructor_bad_degree(self):
        nodes = np.array([
            [0.0, 0.0],
        ])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test__get_degree_valid(self):
        klass = self._get_target_class()

        self.assertEqual(0, klass._get_degree(1))
        self.assertEqual(1, klass._get_degree(3))
        self.assertEqual(2, klass._get_degree(6))
        self.assertEqual(3, klass._get_degree(10))
        self.assertEqual(11, klass._get_degree(78))

    def test__get_degree_invalid(self):
        klass = self._get_target_class()

        with self.assertRaises(ValueError):
            klass._get_degree(2)

        with self.assertRaises(ValueError):
            klass._get_degree(9)

    def test___repr__(self):
        nodes = np.zeros((15, 3))
        surface = self._make_one(nodes)
        expected = '<Surface (degree=4, dimension=3)>'
        self.assertEqual(repr(surface), expected)

    def test_area_property_not_cached(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
            [2.0, 3.0],
        ])
        surface = self._make_one(nodes)
        self.assertIsNone(surface._area)
        with self.assertRaises(NotImplementedError):
            getattr(surface, 'area')

    def test_area_property(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
            [2.0, 3.0],
        ])
        surface = self._make_one(nodes)
        area = 3.14159
        surface._area = area
        self.assertEqual(surface.area, area)

    def _edges_helper(self, edge1, edge2, edge3,
                      nodes1, nodes2, nodes3):
        import bezier

        self.assertIsInstance(edge1, bezier.Curve)
        self.assertTrue(np.all(edge1.nodes == nodes1))

        self.assertIsInstance(edge2, bezier.Curve)
        self.assertTrue(np.all(edge2.nodes == nodes2))

        self.assertIsInstance(edge3, bezier.Curve)
        self.assertTrue(np.all(edge3.nodes == nodes3))

    def test__compute_edges_linear(self):
        nodes = np.array([
            [0.0, 0.0],
            [2.0, 1.0],
            [-3.0, 3.0],
        ])
        p100, p010, p001 = nodes
        surface = self._make_one(nodes)

        edge1, edge2, edge3 = surface._compute_edges()
        self._edges_helper(
            edge1, edge2, edge3,
            np.vstack([p100, p010]),
            np.vstack([p010, p001]),
            np.vstack([p001, p100]))

    def test__compute_edges_quadratic(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.25, 0.5],
            [2.0, 1.0],
            [-1.5, 0.75],
            [0.0, 2.0],
            [-3.0, 3.0],
        ])
        p200, p110, p020, p101, p011, p002 = nodes
        surface = self._make_one(nodes)

        edge1, edge2, edge3 = surface._compute_edges()
        self._edges_helper(
            edge1, edge2, edge3,
            np.vstack([p200, p110, p020]),
            np.vstack([p020, p011, p002]),
            np.vstack([p002, p101, p200]))

    def test__compute_edges_cubic(self):
        nodes = np.array([
            [0.0, 0.0],
            [0.328125, 0.1484375],
            [0.65625, 0.1484375],
            [1.0, 0.0],
            [0.1484375, 0.328125],
            [0.5, 0.5],
            [1.0, 0.53125],
            [0.1484375, 0.65625],
            [0.53125, 1.0],
            [0.0, 1.0],
        ])
        (p300, p210, p120, p030, p201,
         unused_p111, p021, p102, p012, p003) = nodes
        surface = self._make_one(nodes)

        edges = surface._compute_edges()
        self._edges_helper(
            edges[0], edges[1], edges[2],
            np.vstack([p300, p210, p120, p030]),
            np.vstack([p030, p021, p012, p003]),
            np.vstack([p003, p102, p201, p300]))

    def test_edges_property(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        surface = self._make_one(nodes)

        edge1, edge2, edge3 = surface.edges
        nodes1 = nodes[:2, :]
        nodes2 = nodes[1:, :]
        nodes3 = nodes[(2, 0), :]
        self._edges_helper(edge1, edge2, edge3,
                           nodes1, nodes2, nodes3)

    def test_edges_property_cached(self):
        import mock

        surface = self._make_one(np.zeros((3, 2)))

        # Create mock "edges" to be computed.
        sentinel1 = mock.Mock()
        sentinel2 = mock.Mock()
        sentinel3 = mock.Mock()
        expected = sentinel1, sentinel2, sentinel3
        surface._compute_edges = mock.Mock(return_value=expected)

        # Make sure the "edges" when copied just return themselves.
        sentinel1.copy.return_value = sentinel1
        sentinel2.copy.return_value = sentinel2
        sentinel3.copy.return_value = sentinel3

        # Access the property and check the mocks.
        self.assertEqual(surface.edges, expected)
        surface._compute_edges.assert_any_call()
        self.assertEqual(surface._compute_edges.call_count, 1)
        sentinel1.copy.assert_called_once_with()
        sentinel2.copy.assert_called_once_with()
        sentinel3.copy.assert_called_once_with()

        # Access again but make sure no more calls to _compute_edges().
        self.assertEqual(surface.edges, expected)
        self.assertEqual(surface._compute_edges.call_count, 1)

    def test_evaluate_barycentric_linear(self):
        lambda_vals = (0.25, 0.5, 0.25)
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.5],
            [0.0, 1.25],
        ])
        surface = self._make_one(nodes)

        expected = np.array([0.5, 0.5625])
        result = surface.evaluate_barycentric(*lambda_vals)
        self.assertTrue(np.all(expected == result))

    def test_evaluate_barycentric_quadratic(self):
        lambda_vals = (0.0, 0.25, 0.75)
        nodes = np.array([
            [0.0, 0.0],
            [0.5, 0.0],
            [1.0, 0.5],
            [0.5, 1.25],
            [0.0, 1.25],
            [0.0, 0.5],
        ])
        surface = self._make_one(nodes)

        expected = np.array([0.0625, 0.78125])
        result = surface.evaluate_barycentric(*lambda_vals)
        self.assertTrue(np.all(expected == result))

    def test_evaluate_barycentric_cubic(self):
        lambda_vals = (0.125, 0.5, 0.375)
        nodes = np.array([
            [0.0, 0.0],
            [0.25, 0.0],
            [0.75, 0.25],
            [1.0, 0.0],
            [0.0, 0.25],
            [0.375, 0.25],
            [0.5, 0.25],
            [0.0, 0.5],
            [0.25, 0.75],
            [0.0, 1.0],
        ])
        surface = self._make_one(nodes)

        expected = np.array([0.447265625, 0.37060546875])
        result = surface.evaluate_barycentric(*lambda_vals)
        self.assertTrue(np.all(expected == result))

    def test_evaluate_barycentric_negative_weights(self):
        surface = self._make_one(np.zeros((3, 2)))

        lambda_vals = (0.25, -0.5, 1.25)
        self.assertEqual(sum(lambda_vals), 1.0)

        with self.assertRaises(ValueError):
            surface.evaluate_barycentric(*lambda_vals)

    def test_evaluate_barycentric_non_unity_weights(self):
        surface = self._make_one(np.zeros((3, 2)))

        lambda_vals = (0.25, 0.25, 0.25)
        self.assertNotEqual(sum(lambda_vals), 1.0)

        with self.assertRaises(ValueError):
            surface.evaluate_barycentric(*lambda_vals)

    def test_evaluate_barycentric_unsupported(self):
        surface = self._make_one(np.zeros((15, 2)))

        lambda_vals = (1.0, 0.0, 0.0)
        self.assertEqual(sum(lambda_vals), 1.0)

        with self.assertRaises(NotImplementedError):
            surface.evaluate_barycentric(*lambda_vals)

    def test_evaluate_cartesian(self):
        s_t_vals = (0.125, 0.125)
        nodes = np.array([
            [1.0, 1.0],
            [2.0, 1.5],
            [1.0, 2.75],
        ])
        surface = self._make_one(nodes)

        expected = np.array([1.125, 1.28125])
        result = surface.evaluate_cartesian(*s_t_vals)
        self.assertTrue(np.all(expected == result))

    def test_evaluate_cartesian_calls_barycentric(self):
        import mock

        surface = self._make_one(np.zeros((3, 2)))
        eval_method = mock.Mock()
        surface.evaluate_barycentric = eval_method

        s_val = 0.25
        t_val = 0.25
        eval_method.return_value = mock.sentinel.point
        result = surface.evaluate_cartesian(s_val, t_val)
        self.assertIs(result, mock.sentinel.point)

        eval_method.assert_called_once_with(0.5, s_val, t_val)

    def test_evaluate_multi_with_barycentric(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.75],
            [2.0, 1.0],
            [-1.5, 1.0],
            [-0.5, 1.5],
            [-3.0, 2.0],
        ])
        surface = self._make_one(nodes)
        expected = np.array([
            [-1.75, 1.75],
            [0.0, 0.0],
            [0.25, 1.0625],
            [-0.625, 1.046875],
        ])

        param_vals = np.array([
            [0.25, 0.75],
            [0.0, 0.0],
            [0.5, 0.25],
            [0.25, 0.375],
        ])
        result = surface.evaluate_multi(param_vals)
        self.assertTrue(np.all(result == expected))

    def test_evaluate_multi_with_cartesian(self):
        nodes = np.array([
            [0.0, 0.0],
            [2.0, 1.0],
            [-3.0, 2.0],
        ])
        surface = self._make_one(nodes)
        expected = np.array([
            [0.0, 0.0],
            [2.0, 1.0],
            [-0.5, 1.5],
        ])

        param_vals = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.5, 0.5],
        ])
        result = surface.evaluate_multi(param_vals)
        self.assertTrue(np.all(result == expected))

    def test_evaluate_multi_wrong_dimension(self):
        surface = self._make_one(np.zeros((3, 2)))
        param_vals_1d = np.zeros((4,))
        with self.assertRaises(ValueError):
            surface.evaluate_multi(param_vals_1d)

    def test_evaluate_multi_wrong_param_cols(self):
        surface = self._make_one(np.zeros((3, 2)))
        param_vals = np.zeros((4, 4))
        with self.assertRaises(ValueError):
            surface.evaluate_multi(param_vals)

    def test__add_patch(self):
        import mock

        klass = self._get_target_class()

        ax = mock.Mock()
        color = (0.5, 0.0, 0.5)
        all_nodes = np.array([
            [0.0, 1.0],
            [1.0, 3.0],
            [2.0, 6.0],
            [3.0, 10.0],
            [4.0, 15.0],
            [5.0, 21.0],
        ])
        edge1 = all_nodes[(-1, 0, 1), :]
        edge2 = all_nodes[(1, 2, 3), :]
        edge3 = all_nodes[(3, 4, 5), :]
        self.assertIsNone(
            klass._add_patch(ax, color, edge1, edge2, edge3))

        # Check the call to ax.add_patch(). We can't
        # assert_called_once_with() since == breaks on NumPy arrays.
        self.assertEqual(ax.add_patch.call_count, 1)
        call = ax.add_patch.mock_calls[0]
        # Unpack the call as name, positional args, keyword args
        _, positional, keyword = call
        self.assertEqual(keyword, {})
        self.assertEqual(len(positional), 1)
        patch = positional[0]
        self.assertTrue(np.all(patch.get_path().vertices == all_nodes))

    def _check_plot_call(self, call, expected, **kwargs):
        # Unpack the call as name, positional args, keyword args
        _, positional, keyword = call
        self.assertEqual(keyword, kwargs)
        self.assertEqual(len(positional), 2)
        self.assertTrue(np.all(positional[0] == expected[:, 0]))
        self.assertTrue(np.all(positional[1] == expected[:, 1]))

    def _plot_helper(self, show=False):
        import matplotlib.lines
        import mock

        nodes = np.array([
            [0.0, 0.0],
            [0.0, 1.0],
            [1.0, 0.0],
        ])
        curve = self._make_one(nodes)
        plt = mock.Mock()

        figure = mock.Mock()
        plt.figure.return_value = figure
        ax = mock.Mock()
        figure.gca.return_value = ax

        color = (0.5, 0.5, 0.75)
        line = matplotlib.lines.Line2D([], [], color=color)
        ax.plot.return_value = (line,)

        with mock.patch('bezier.surface.plt', new=plt):
            if show:
                result = curve.plot(2, show=True)
            else:
                result = curve.plot(2)

        self.assertIs(result, figure)

        # Check mocks.
        plt.figure.assert_called_once_with()
        figure.gca.assert_called_once_with()

        # Check the calls to ax.plot(). We can't assert_any_call()
        # since == breaks on NumPy arrays.
        self.assertEqual(ax.plot.call_count, 4)
        calls = ax.plot.mock_calls
        self._check_plot_call(calls[0], nodes[:2, :])
        self._check_plot_call(calls[1], nodes[1:, :], color=color)
        self._check_plot_call(calls[2], nodes[(2, 0), :], color=color)
        self._check_plot_call(calls[3], nodes,
                              color='black', marker='o', linestyle='None')
        # Check the calls to ax.add_patch().
        self.assertEqual(ax.add_patch.call_count, 1)

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
            [2.0, 6.0, 9.0],
        ])
        surface = self._make_one(nodes)
        with self.assertRaises(NotImplementedError):
            surface.plot(32)

    def _subdivide_helper(self, nodes, expected_a, expected_b,
                          expected_c, expected_d):
        klass = self._get_target_class()

        surface = self._make_one(nodes)
        surface_a, surface_b, surface_c, surface_d = surface.subdivide()

        self.assertIsInstance(surface_a, klass)
        self.assertTrue(np.all(surface_a._nodes == expected_a))
        self.assertIsInstance(surface_b, klass)
        self.assertTrue(np.all(surface_b._nodes == expected_b))
        self.assertIsInstance(surface_c, klass)
        self.assertTrue(np.all(surface_c._nodes == expected_c))
        self.assertIsInstance(surface_d, klass)
        self.assertTrue(np.all(surface_d._nodes == expected_d))

    def _subdivide_points_check(self, surface, pts_exponent=5):
        # Using the exponent means that we will divide by
        # 2**exp, which can be done without roundoff (for small
        # enough exponents).
        sub_surfaces = surface.subdivide()

        ref_triangle = utils.ref_triangle_uniform_nodes(pts_exponent)
        quarter_a = 0.5 * ref_triangle
        quarters = [
            quarter_a,
            np.array([0.5, 0.5]) - quarter_a,  # B
            quarter_a + np.array([0.5, 0.0]),  # C
            quarter_a + np.array([0.0, 0.5]),  # D
        ]

        for sub_surface, quarter in zip(sub_surfaces, quarters):
            # Make sure sub_surface(ref_triangle) == surface(quarter)
            main_vals = surface.evaluate_multi(quarter)
            sub_vals = sub_surface.evaluate_multi(ref_triangle)
            self.assertTrue(np.all(main_vals == sub_vals))

    def test_subdivide_linear(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 0.0],
            [0.0, 1.0],
        ])
        expected_a = np.array([
            [0.0, 0.0],
            [0.5, 0.0],
            [0.0, 0.5],
        ])
        expected_b = np.array([
            [0.5, 0.5],
            [0.0, 0.5],
            [0.5, 0.0],
        ])
        expected_c = np.array([
            [0.5, 0.0],
            [1.0, 0.0],
            [0.5, 0.5],
        ])
        expected_d = np.array([
            [0.0, 0.5],
            [0.5, 0.5],
            [0.0, 1.0],
        ])
        self._subdivide_helper(nodes, expected_a, expected_b,
                               expected_c, expected_d)

    def test_subdivide_line_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(3, 2), seed=123987, num_bits=8)

        surface = self._make_one(nodes)
        self.assertEqual(surface.degree, 1)
        self._subdivide_points_check(surface)

    def test_subdivide_quadratic(self):
        nodes = np.array([
            [0.0, 0.0],
            [0.5, 0.25],
            [1.0, 0.0],
            [0.5, 0.75],
            [0.0, 1.0],
            [0.0, 0.5],
        ])
        expected_a = np.array([
            [0.0, 0.0],
            [0.25, 0.125],
            [0.5, 0.125],
            [0.25, 0.375],
            [0.25, 0.5],
            [0.25, 0.5],
        ])
        expected_b = np.array([
            [0.25, 0.625],
            [0.25, 0.625],
            [0.25, 0.5],
            [0.5, 0.5],
            [0.25, 0.5],
            [0.5, 0.125],
        ])
        expected_c = np.array([
            [0.5, 0.125],
            [0.75, 0.125],
            [1.0, 0.0],
            [0.5, 0.5],
            [0.5, 0.5],
            [0.25, 0.625],
        ])
        expected_d = np.array([
            [0.25, 0.5],
            [0.25, 0.625],
            [0.25, 0.625],
            [0.25, 0.625],
            [0.0, 0.75],
            [0.0, 0.5],
        ])
        self._subdivide_helper(nodes, expected_a, expected_b,
                               expected_c, expected_d)

    def test_subdivide_quadratic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(6, 2), seed=45001, num_bits=8)

        surface = self._make_one(nodes)
        self.assertEqual(surface.degree, 2)
        self._subdivide_points_check(surface)

    def test_subdivide_cubic(self):
        nodes = np.array([
            [0.0, 0.0],
            [3.25, 1.5],
            [6.5, 1.5],
            [10.0, 0.0],
            [1.5, 3.25],
            [5.0, 5.0],
            [10.0, 5.25],
            [1.5, 6.5],
            [5.25, 10.0],
            [0.0, 10.0],
        ])
        expected_a = np.array([
            [0.0, 0.0],
            [1.625, 0.75],
            [3.25, 1.125],
            [4.90625, 1.125],
            [0.75, 1.625],
            [2.4375, 2.4375],
            [4.3125, 2.875],
            [1.125, 3.25],
            [2.875, 4.3125],
            [1.125, 4.90625],
        ])
        expected_b = np.array([
            [6.96875, 6.96875],
            [4.8125, 6.65625],
            [2.875, 5.96875],
            [1.125, 4.90625],
            [6.65625, 4.8125],
            [4.75, 4.75],
            [2.875, 4.3125],
            [5.96875, 2.875],
            [4.3125, 2.875],
            [4.90625, 1.125],
        ])
        expected_c = np.array([
            [4.90625, 1.125],
            [6.5625, 1.125],
            [8.25, 0.75],
            [10.0, 0.0],
            [5.96875, 2.875],
            [7.875, 2.9375],
            [10.0, 2.625],
            [6.65625, 4.8125],
            [8.8125, 5.125],
            [6.96875, 6.96875],
        ])
        expected_d = np.array([
            [1.125, 4.90625],
            [2.875, 5.96875],
            [4.8125, 6.65625],
            [6.96875, 6.96875],
            [1.125, 6.5625],
            [2.9375, 7.875],
            [5.125, 8.8125],
            [0.75, 8.25],
            [2.625, 10.0],
            [0.0, 10.0],
        ])
        self._subdivide_helper(nodes, expected_a, expected_b,
                               expected_c, expected_d)

    def test_subdivide_cubic_check_evaluate(self):
        # Use a fixed seed so the test is deterministic and round
        # the nodes to 8 bits of precision to avoid round-off.
        nodes = utils.get_random_nodes(
            shape=(10, 2), seed=346323, num_bits=8)

        surface = self._make_one(nodes)
        self.assertEqual(surface.degree, 3)
        self._subdivide_points_check(surface)

    def test_subdivide_unsupported_degree(self):
        nodes = np.zeros((78, 3))
        surface = self._make_one(nodes)
        with self.assertRaises(NotImplementedError):
            surface.subdivide()

    def test_copy(self):
        import mock

        surface = self._make_one(np.zeros((6, 2)))
        fake_nodes = mock.Mock()
        surface._nodes = fake_nodes

        copied_nodes = np.zeros((3, 2))
        fake_nodes.copy.return_value = copied_nodes

        new_surface = surface.copy()
        self.assertIsInstance(new_surface, self._get_target_class())
        self.assertIs(new_surface._nodes, copied_nodes)

        fake_nodes.copy.assert_called_once_with()
