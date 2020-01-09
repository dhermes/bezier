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


class TestSurface(utils.NumPyTestCase):

    REF_TRIANGLE = utils.ref_triangle_uniform_nodes(5)
    REF_TRIANGLE3 = utils.ref_triangle_uniform_nodes(3)
    QUADRATIC = np.asfortranarray(
        [[0.0, 1.25, 2.0, -1.5, 0.0, -3.0], [0.0, 0.5, 1.0, 0.75, 2.0, 3.0]]
    )
    UNIT_TRIANGLE = np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
    ZEROS = np.zeros((2, 3), order="F")

    @staticmethod
    def _get_target_class():
        from bezier import surface

        return surface.Surface

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def _make_one_no_slots(self, *args, **kwargs):
        class NoSlots(self._get_target_class()):
            pass

        return NoSlots(*args, **kwargs)

    def test_constructor(self):
        nodes = np.asfortranarray([[0.0, 0.625, 1.0], [0.0, 0.5, 0.75]])
        surface = self._make_one(nodes, 1, copy=False)
        self.assertEqual(surface._degree, 1)
        self.assertEqual(surface._dimension, 2)
        self.assertIs(surface._nodes, nodes)
        self.assertIsNone(surface._edges)

    def test_constructor_wrong_dimension(self):
        nodes = np.asfortranarray([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes, 0)
        nodes = np.zeros((3, 2, 2), order="F")
        with self.assertRaises(ValueError):
            self._make_one(nodes, 1)

    def test_constructor_invalid_degree(self):
        nodes = np.empty((1, 6), order="F")
        with self.assertRaises(ValueError) as exc_info:
            self._make_one(nodes, 1)

        exc_args = exc_info.exception.args
        self.assertEqual(
            exc_args, ("A degree 1 surface should have 3 nodes, not 6.",)
        )

    def test_from_nodes_factory(self):
        nodes = np.asfortranarray(
            [
                [0.0, 1.0, 2.0, 0.0, 1.0, 2.0],
                [0.0, 0.5, 0.0, 1.0, 1.0, 3.0],
                [0.0, 0.0, 0.0, 2.0, 0.0, 0.0],
            ]
        )
        klass = self._get_target_class()
        surface = klass.from_nodes(nodes)
        self.assertIsInstance(surface, klass)
        self.assertEqual(surface._degree, 2)
        self.assertEqual(surface._dimension, 3)
        self.assertEqual(surface._nodes, nodes)
        self.assertIsNone(surface._edges)

    def test_from_nodes_factory_non_array(self):
        nodes = [[0.0, 1.0, 2.0]]
        klass = self._get_target_class()
        surface = klass.from_nodes(nodes)
        self.assertIsInstance(surface, klass)
        self.assertEqual(surface._degree, 1)
        self.assertEqual(surface._dimension, 1)
        self.assertTrue(np.all(surface._nodes == nodes))
        self.assertIsNone(surface._edges)

    def test_from_nodes_factory_invalid_degree(self):
        klass = self._get_target_class()
        messages = {
            2: "A degree 1 surface should have 3 nodes, not 2.",
            9: "A degree 3 surface should have 10 nodes, not 9.",
        }
        for num_nodes, message in messages.items():
            with self.assertRaises(ValueError) as exc_info:
                nodes = np.empty((1, num_nodes), order="F")
                klass.from_nodes(nodes)

            exc_args = exc_info.exception.args
            self.assertEqual(exc_args, (message,))

    def test___repr__(self):
        nodes = np.zeros((3, 15), order="F")
        surface = self._make_one(nodes, 4)
        expected = "<Surface (degree=4, dimension=3)>"
        self.assertEqual(repr(surface), expected)

    def test__get_degree_valid(self):
        klass = self._get_target_class()
        self.assertEqual(0, klass._get_degree(1))
        self.assertEqual(1, klass._get_degree(3))
        self.assertEqual(2, klass._get_degree(6))
        self.assertEqual(3, klass._get_degree(10))
        self.assertEqual(11, klass._get_degree(78))

    def test_area_property_wrong_dimension(self):
        nodes = np.asfortranarray(
            [[0.0, 0.0, 0.0], [1.0, 2.0, 0.0], [2.0, 3.0, 0.0]]
        )
        surface = self._make_one(nodes, 1)
        with self.assertRaises(NotImplementedError) as exc_info:
            getattr(surface, "area")

        exc_args = exc_info.exception.args
        expected_args = (
            "2D is the only supported dimension",
            "Current dimension",
            3,
        )
        self.assertEqual(exc_args, expected_args)

    def test_area_property(self):
        nodes = np.asfortranarray([[0.0, 2.0, 1.0], [0.0, 3.0, 2.0]])
        surface = self._make_one(nodes, 1)
        self.assertEqual(surface.area, 0.5)

    def _edges_helper(self, edge1, edge2, edge3, nodes1, nodes2, nodes3):
        import bezier

        self.assertIsInstance(edge1, bezier.Curve)
        self.assertEqual(edge1._nodes, nodes1)
        self.assertIsInstance(edge2, bezier.Curve)
        self.assertEqual(edge2._nodes, nodes2)
        self.assertIsInstance(edge3, bezier.Curve)
        self.assertEqual(edge3._nodes, nodes3)

    def test__compute_edges_linear(self):
        nodes = np.asfortranarray([[0.0, 2.0, -3.0], [0.0, 1.0, 3.0]])
        p100, p010, p001 = nodes.T
        surface = self._make_one(nodes, 1)
        edge1, edge2, edge3 = surface._compute_edges()
        nodes1 = np.asfortranarray(np.vstack([p100, p010]).T)
        nodes2 = np.asfortranarray(np.vstack([p010, p001]).T)
        nodes3 = np.asfortranarray(np.vstack([p001, p100]).T)
        self._edges_helper(edge1, edge2, edge3, nodes1, nodes2, nodes3)

    def test__compute_edges_quadratic(self):
        nodes = self.QUADRATIC
        p200, p110, p020, p101, p011, p002 = nodes.T
        surface = self._make_one(nodes, 2)
        edges = surface._compute_edges()
        nodes1 = np.asfortranarray(np.vstack([p200, p110, p020]).T)
        nodes2 = np.asfortranarray(np.vstack([p020, p011, p002]).T)
        nodes3 = np.asfortranarray(np.vstack([p002, p101, p200]).T)
        self._edges_helper(
            edges[0], edges[1], edges[2], nodes1, nodes2, nodes3
        )

    def test__compute_edges_cubic(self):
        nodes = np.asfortranarray(
            [
                [
                    0.0,
                    0.328125,
                    0.65625,
                    1.0,
                    0.1484375,
                    0.5,
                    1.0,
                    0.1484375,
                    0.53125,
                    0.0,
                ],
                [
                    0.0,
                    0.1484375,
                    0.1484375,
                    0.0,
                    0.328125,
                    0.5,
                    0.53125,
                    0.65625,
                    1.0,
                    1.0,
                ],
            ]
        )
        (
            p300,
            p210,
            p120,
            p030,
            p201,
            unused_p111,
            p021,
            p102,
            p012,
            p003,
        ) = nodes.T
        surface = self._make_one(nodes, 3)
        edges = surface._compute_edges()
        self._edges_helper(
            edges[0],
            edges[1],
            edges[2],
            np.asfortranarray(np.vstack([p300, p210, p120, p030]).T),
            np.asfortranarray(np.vstack([p030, p021, p012, p003]).T),
            np.asfortranarray(np.vstack([p003, p102, p201, p300]).T),
        )

    def test__get_edges(self):
        surface = self._make_one_no_slots(self.ZEROS, 1)
        compute_mock = unittest.mock.Mock(
            return_value=unittest.mock.sentinel.edges
        )
        surface._compute_edges = compute_mock
        self.assertIsNone(surface._edges)
        self.assertIs(surface._get_edges(), unittest.mock.sentinel.edges)
        self.assertIs(surface._edges, unittest.mock.sentinel.edges)
        compute_mock.assert_called_once_with()

    def test__get_edges_cached(self):
        surface = self._make_one_no_slots(self.ZEROS, 1)
        compute_mock = unittest.mock.Mock(spec=[])
        surface._compute_edges = compute_mock
        surface._edges = unittest.mock.sentinel.edges
        self.assertIs(surface._get_edges(), unittest.mock.sentinel.edges)
        compute_mock.assert_not_called()

    def test_edges_property(self):
        nodes = self.UNIT_TRIANGLE
        surface = self._make_one(nodes, 1)
        edge1, edge2, edge3 = surface.edges
        nodes1 = nodes[:, :2]
        nodes2 = nodes[:, 1:]
        nodes3 = np.asfortranarray(nodes[:, (2, 0)])
        self._edges_helper(edge1, edge2, edge3, nodes1, nodes2, nodes3)

    def test_edges_property_cached(self):
        surface = self._make_one_no_slots(self.ZEROS, 1)
        # Create mock "edges" to be computed.
        sentinel1 = unittest.mock.Mock(spec=["copy"])
        sentinel2 = unittest.mock.Mock(spec=["copy"])
        sentinel3 = unittest.mock.Mock(spec=["copy"])
        expected = sentinel1, sentinel2, sentinel3
        surface._compute_edges = unittest.mock.Mock(return_value=expected)
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

    def test__verify_barycentric(self):
        klass = self._get_target_class()
        # Valid inside.
        self.assertIsNone(klass._verify_barycentric(0.5, 0.25, 0.25))
        # Valid boundary.
        self.assertIsNone(klass._verify_barycentric(0.5, 0.0, 0.5))
        self.assertIsNone(klass._verify_barycentric(0.25, 0.75, 0.0))
        self.assertIsNone(klass._verify_barycentric(0.0, 0.0, 1.0))
        self.assertIsNone(klass._verify_barycentric(0.0, 0.5, 0.5))
        # Invalid sum.
        with self.assertRaises(ValueError):
            klass._verify_barycentric(0.5, 0.5, 0.5)
        # Invalid lamdba1
        with self.assertRaises(ValueError):
            klass._verify_barycentric(-0.5, 0.75, 0.75)
        # Invalid lamdba2.
        with self.assertRaises(ValueError):
            klass._verify_barycentric(0.75, -0.5, 0.75)
        # Invalid lamdba3.
        with self.assertRaises(ValueError):
            klass._verify_barycentric(0.875, 0.25, -0.125)

    def test_evaluate_barycentric(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1, copy=False)
        lambda_vals = (0.25, 0.0, 0.75)
        # Just make sure we call the helper.
        patch = unittest.mock.patch(
            "bezier._surface_helpers.evaluate_barycentric",
            return_value=unittest.mock.sentinel.evaluated,
        )
        with patch as mocked:
            result = surface.evaluate_barycentric(*lambda_vals)
            self.assertIs(result, unittest.mock.sentinel.evaluated)
            mocked.assert_called_once_with(self.UNIT_TRIANGLE, 1, *lambda_vals)

    def test_evaluate_barycentric_negative_weights_no_verify(self):
        lambda_vals = (0.25, -0.5, 1.25)
        nodes = np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.5, 1.25]])
        surface = self._make_one(nodes, 1)
        self.assertLess(min(lambda_vals), 0.0)
        result = surface.evaluate_barycentric(*lambda_vals, _verify=False)
        expected = np.asfortranarray([[-0.5], [1.3125]])
        self.assertEqual(result, expected)

    def test_evaluate_barycentric_non_unity_weights_no_verify(self):
        lambda_vals = (0.25, 0.25, 0.25)
        nodes = np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.5, 1.25]])
        surface = self._make_one(nodes, 1)
        self.assertNotEqual(sum(lambda_vals), 1.0)
        result = surface.evaluate_barycentric(*lambda_vals, _verify=False)
        expected = np.asfortranarray([[0.25], [0.4375]])
        self.assertEqual(result, expected)

    def test_evaluate_barycentric_multi_wrong_dimension(self):
        surface = self._make_one(self.ZEROS, 1)
        param_vals_1d = np.zeros((4,), order="F")
        with self.assertRaises(ValueError):
            surface.evaluate_barycentric_multi(param_vals_1d)

    def _eval_bary_multi_helper(self, **kwargs):
        nodes = np.asfortranarray([[0.0, 2.0, -3.0], [0.0, 1.0, 2.0]])
        surface = self._make_one(nodes, 1, copy=False)
        param_vals = np.asfortranarray([[1.0, 0.0, 0.0]])
        patch = unittest.mock.patch(
            "bezier._surface_helpers.evaluate_barycentric_multi",
            return_value=unittest.mock.sentinel.evaluated,
        )
        with patch as mocked:
            result = surface.evaluate_barycentric_multi(param_vals, **kwargs)
            self.assertEqual(result, unittest.mock.sentinel.evaluated)
            mocked.assert_called_once_with(nodes, 1, param_vals, 2)

    def test_evaluate_barycentric_multi(self):
        self._eval_bary_multi_helper()

    def test_evaluate_barycentric_multi_no_verify(self):
        self._eval_bary_multi_helper(_verify=False)

    def test__verify_cartesian(self):
        klass = self._get_target_class()
        # Valid inside.
        self.assertIsNone(klass._verify_cartesian(0.25, 0.25))
        # Valid boundary.
        self.assertIsNone(klass._verify_cartesian(0.0, 0.5))
        self.assertIsNone(klass._verify_cartesian(0.75, 0.0))
        self.assertIsNone(klass._verify_cartesian(0.0, 1.0))
        self.assertIsNone(klass._verify_cartesian(0.5, 0.5))
        # Invalid s.
        with self.assertRaises(ValueError):
            klass._verify_cartesian(-0.5, 0.75)
        # Invalid t.
        with self.assertRaises(ValueError):
            klass._verify_cartesian(0.25, -0.125)
        # Invalid (1 - s - t).
        with self.assertRaises(ValueError):
            klass._verify_cartesian(0.75, 0.75)

    def test_evaluate_cartesian(self):
        s_t_vals = (0.125, 0.125)
        nodes = np.asfortranarray([[1.0, 2.0, 1.0], [1.0, 1.5, 2.75]])
        surface = self._make_one(nodes, 1)
        expected = np.asfortranarray([[1.125], [1.28125]])
        result = surface.evaluate_cartesian(*s_t_vals)
        self.assertEqual(result, expected)

    def test_evaluate_cartesian_no_verify(self):
        s_t_vals = (0.25, 1.0)
        nodes = np.asfortranarray([[1.0, 2.0, 1.0], [1.0, 1.5, 2.75]])
        surface = self._make_one(nodes, 1)
        expected = np.asfortranarray([[1.25], [2.875]])
        result = surface.evaluate_cartesian(*s_t_vals, _verify=False)
        self.assertEqual(result, expected)

    def test_evaluate_cartesian_calls_helper(self):
        nodes = self.ZEROS
        surface = self._make_one_no_slots(nodes, 1, copy=False)
        patch = unittest.mock.patch(
            "bezier._surface_helpers.evaluate_barycentric",
            return_value=unittest.mock.sentinel.point,
        )
        s_val = 0.25
        t_val = 0.25
        with patch as mocked:
            result = surface.evaluate_cartesian(s_val, t_val)
            self.assertIs(result, unittest.mock.sentinel.point)
            mocked.assert_called_once_with(nodes, 1, 0.5, s_val, t_val)

    def test_evaluate_cartesian_multi_wrong_dimension(self):
        surface = self._make_one(self.ZEROS, 1)
        param_vals_1d = np.zeros((4,), order="F")
        with self.assertRaises(ValueError):
            surface.evaluate_cartesian_multi(param_vals_1d)

    def _eval_cartesian_multi_helper(self, **kwargs):
        nodes = np.asfortranarray([[2.0, 0.0, 3.0], [3.0, 2.0, 7.5]])
        surface = self._make_one(nodes, 1, copy=False)
        param_vals = np.asfortranarray([[1.0, 0.0]])
        patch = unittest.mock.patch(
            "bezier._surface_helpers.evaluate_cartesian_multi",
            return_value=unittest.mock.sentinel.evaluated,
        )
        with patch as mocked:
            result = surface.evaluate_cartesian_multi(param_vals, **kwargs)
            self.assertEqual(result, unittest.mock.sentinel.evaluated)
            mocked.assert_called_once_with(nodes, 1, param_vals, 2)

    def test_evaluate_cartesian_multi(self):
        self._eval_cartesian_multi_helper()

    def test_evaluate_cartesian_multi_no_verify(self):
        self._eval_cartesian_multi_helper(_verify=False)

    def test_plot_wrong_dimension(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0, 2.0], [0.0, 3.0, 6.0], [0.0, 4.0, 9.0]]
        )
        surface = self._make_one(nodes, 1, copy=False)
        with self.assertRaises(NotImplementedError):
            surface.plot(32)

    @unittest.mock.patch("bezier._plot_helpers.new_axis")
    @unittest.mock.patch("bezier._plot_helpers.add_patch")
    def test_plot_defaults(self, add_patch_mock, new_axis_mock):
        ax = unittest.mock.Mock(spec=[])
        new_axis_mock.return_value = ax
        curve = self._make_one(self.UNIT_TRIANGLE, 1, copy=False)
        pts_per_edge = 16
        result = curve.plot(pts_per_edge)
        self.assertIs(result, ax)
        # Verify mocks.
        new_axis_mock.assert_called_once_with()
        add_patch_mock.assert_called_once_with(
            ax, None, pts_per_edge, *curve._edges
        )

    @unittest.mock.patch("bezier._plot_helpers.new_axis")
    @unittest.mock.patch("bezier._plot_helpers.add_patch")
    def test_plot_explicit(self, add_patch_mock, new_axis_mock):
        ax = unittest.mock.Mock(spec=["plot"])
        color = (0.5, 0.5, 0.5)
        curve = self._make_one(self.UNIT_TRIANGLE, 1, copy=False)
        pts_per_edge = 16
        result = curve.plot(pts_per_edge, color=color, ax=ax, with_nodes=True)
        self.assertIs(result, ax)
        # Verify mocks.
        new_axis_mock.assert_not_called()
        add_patch_mock.assert_called_once_with(
            ax, color, pts_per_edge, *curve._edges
        )
        # Check the call to ax.plot(). We can't assert_any_call()
        # since == breaks on NumPy arrays.
        self.assertEqual(ax.plot.call_count, 1)
        call = ax.plot.mock_calls[0]
        utils.check_plot_call(
            self,
            call,
            self.UNIT_TRIANGLE,
            color="black",
            marker="o",
            linestyle="None",
        )

    def test_subdivide(self):
        klass = self._get_target_class()
        degree = 1
        surface = self._make_one(self.UNIT_TRIANGLE, degree, copy=False)
        surface_a, surface_b, surface_c, surface_d = surface.subdivide()
        # Check sub-surface A.
        self.assertIsInstance(surface_a, klass)
        self.assertEqual(surface_a._degree, degree)
        expected_a = np.asfortranarray([[0.0, 0.5, 0.0], [0.0, 0.0, 0.5]])
        self.assertEqual(surface_a._nodes, expected_a)
        # Check sub-surface B.
        self.assertIsInstance(surface_b, klass)
        self.assertEqual(surface_b._degree, degree)
        expected_b = np.asfortranarray([[0.5, 0.0, 0.5], [0.5, 0.5, 0.0]])
        self.assertEqual(surface_b._nodes, expected_b)
        # Check sub-surface C.
        self.assertIsInstance(surface_c, klass)
        self.assertEqual(surface_c._degree, degree)
        expected_c = np.asfortranarray([[0.5, 1.0, 0.5], [0.0, 0.0, 0.5]])
        self.assertEqual(surface_c._nodes, expected_c)
        # Check sub-surface D.
        self.assertIsInstance(surface_d, klass)
        self.assertEqual(surface_d._degree, degree)
        expected_d = np.asfortranarray([[0.0, 0.5, 0.0], [0.5, 0.5, 1.0]])
        self.assertEqual(surface_d._nodes, expected_d)

    def test__compute_valid_bad_dimension(self):
        nodes = np.zeros((3, 6), order="F")
        surface = self._make_one(nodes, 2)
        with self.assertRaises(NotImplementedError):
            surface._compute_valid()

    def test__compute_valid_linear(self):
        # Positively signed Jacobian (i.e. edges have inward normals).
        nodes = self.UNIT_TRIANGLE
        surface = self._make_one(nodes, 1)
        self.assertTrue(surface._compute_valid())
        # Change the nodes from counterclockwise to clockwise, so the
        # Jacobian becomes negatively signed.
        nodes = np.asfortranarray(nodes[:, (1, 0, 2)])
        surface = self._make_one(nodes, 1, copy=False)
        self.assertFalse(surface._compute_valid())
        # Collinear.
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 2.0, 4.0]])
        surface = self._make_one(nodes, 1, copy=False)
        self.assertFalse(surface._compute_valid())

    def test__compute_valid_quadratic(self):
        # Positively signed Jacobian (i.e. edges have inward normals).
        nodes = np.asfortranarray(
            [
                [0.0, 0.5, 1.0, 0.1875, 0.625, 0.0],
                [0.0, -0.1875, 0.0, 0.5, 0.625, 1.0],
            ]
        )
        surface = self._make_one(nodes, 2, copy=False)
        self.assertTrue(surface._compute_valid())
        # Change the nodes from counterclockwise to clockwise, so the
        # Jacobian becomes negatively signed.
        nodes = np.asfortranarray(nodes[:, (2, 1, 0, 4, 3, 5)])
        surface = self._make_one(nodes, 2, copy=False)
        self.assertFalse(surface._compute_valid())
        # Mixed sign Jacobian: B(L1, L2, L3) = [L1^2 + L2^2, L2^2 + L3^2]
        nodes = np.asfortranarray(
            [[1.0, 0.0, 1.0, 0.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]]
        )
        surface = self._make_one(nodes, 2, copy=False)
        self.assertFalse(surface._compute_valid())

    def test__compute_valid_cubic(self):
        # Positively signed Jacobian (i.e. edges have inward normals).
        nodes = np.asfortranarray(
            [
                [0.0, 1.0, 2.0, 3.0, 0.0, 1.0, 2.25, 0.0, 1.25, 0.0],
                [0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.25, 2.0, 2.25, 3.0],
            ]
        )
        surface = self._make_one(nodes, 3, copy=False)
        self.assertTrue(surface._compute_valid())
        # Change the nodes from counterclockwise to clockwise, so the
        # Jacobian becomes negatively signed.
        nodes = np.asfortranarray(nodes[:, (3, 2, 1, 0, 6, 5, 4, 8, 7, 9)])
        surface = self._make_one(nodes, 3, copy=False)
        self.assertFalse(surface._compute_valid())
        # Mixed sign Jacobian: B(L1, L2, L3) = [L1^3 + L2^3, L2^3 + L3^3]
        nodes = np.asfortranarray(
            [
                [1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0],
            ]
        )
        surface = self._make_one(nodes, 3, copy=False)
        self.assertFalse(surface._compute_valid())

    def test__compute_valid_bad_degree(self):
        from bezier import _py_helpers

        degree = 4
        num_nodes = ((degree + 1) * (degree + 2)) // 2
        nodes = np.zeros((2, num_nodes), order="F")
        surface = self._make_one(nodes, degree=degree)
        with self.assertRaises(_py_helpers.UnsupportedDegree) as exc_info:
            surface._compute_valid()
        self.assertEqual(exc_info.exception.degree, degree)
        self.assertEqual(exc_info.exception.supported, (1, 2, 3))

    def test_is_valid_property(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1)
        self.assertTrue(surface.is_valid)

    def test___dict___property(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1, copy=False)
        props_dict = surface.__dict__
        expected = {
            "_nodes": self.UNIT_TRIANGLE,
            "_dimension": 2,
            "_degree": 1,
            "_edges": None,
        }
        self.assertEqual(props_dict, expected)
        # Check that modifying ``props_dict`` won't modify ``surface``.
        expected["_dimension"] = -42
        self.assertNotEqual(surface._dimension, expected["_dimension"])

    def test_locate(self):
        surface = self._make_one(self.QUADRATIC, 2)
        point = surface.evaluate_cartesian(0.5, 0.25)
        s, t = surface.locate(point)
        self.assertEqual(s, 0.5)
        self.assertEqual(t, 0.25)

    def test_locate_no_verify(self):
        surface = self._make_one(self.QUADRATIC, 2)
        s = 0.125
        t = 0.125
        x_val, y_val = surface.evaluate_cartesian(s, t).flatten()
        point = np.asfortranarray([[x_val, np.nan], [y_val, np.nan]])
        # Make sure it fails.
        with self.assertRaises(ValueError):
            surface.locate(point)
        # Will only use the first row if _verify=False.
        computed_s, computed_t = surface.locate(point, _verify=False)
        self.assertEqual(s, computed_s)
        self.assertEqual(t, computed_t)

    def test_locate_bad_dimension(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0]])
        surface = self._make_one(nodes, 1)
        with self.assertRaises(NotImplementedError):
            surface.locate(None)

    def test_locate_bad_point(self):
        surface = self._make_one(self.QUADRATIC, 2)
        point1 = np.asfortranarray([0.0, 1.0])
        point2 = np.asfortranarray([[0.0, 1.0, 2.0]])
        with self.assertRaises(ValueError):
            surface.locate(point1)
        with self.assertRaises(ValueError):
            surface.locate(point2)

    def _basic_intersect_helper(self, **kwargs):
        import bezier

        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        # Similar triangle with overlapping square.
        nodes = np.asfortranarray([[0.5, 0.5, -0.5], [0.0, 1.0, 1.0]])
        surface2 = self._make_one(nodes, 1)
        intersections = surface1.intersect(surface2, **kwargs)
        self.assertEqual(len(intersections), 1)
        intersection = intersections[0]
        self.assertIsInstance(intersection, bezier.CurvedPolygon)
        self.assertEqual(
            intersection._metadata,  # pylint: disable=no-member
            ((5, 0.5, 1.0), (3, 0.0, 0.5), (1, 0.5, 1.0), (2, 0.0, 0.5)),
        )
        cp_edges = intersection._edges
        self.assertEqual(len(cp_edges), 4)
        # Edge 0.
        expected = np.asfortranarray([[0.0, 0.5], [0.5, 0.0]])
        self.assertEqual(cp_edges[0]._nodes, expected)
        # Edge 1.
        expected = np.asfortranarray([[0.5, 0.5], [0.0, 0.5]])
        self.assertEqual(cp_edges[1]._nodes, expected)
        # Edge 2.
        expected = np.asfortranarray([[0.5, 0.0], [0.5, 1.0]])
        self.assertEqual(cp_edges[2]._nodes, expected)
        # Edge 3.
        expected = np.asfortranarray([[0.0, 0.0], [1.0, 0.5]])
        self.assertEqual(cp_edges[3]._nodes, expected)

    def test_intersect(self):
        self._basic_intersect_helper()

    def test_intersect_no_verify(self):
        self._basic_intersect_helper(_verify=False)

    def test_intersect_bad_strategy(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1)
        strategy = unittest.mock.sentinel.bad_strategy
        with self.assertRaises(ValueError) as exc_info:
            surface.intersect(surface, strategy=strategy)
        exc_args = exc_info.exception.args
        self.assertEqual(exc_args, ("Unexpected strategy.", strategy))

    def test_intersect_algebraic(self):
        from bezier import _py_intersection_helpers

        strategy = _py_intersection_helpers.IntersectionStrategy.ALGEBRAIC
        self._basic_intersect_helper(strategy=strategy)

    def test_intersect_disjoint_bbox(self):
        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        nodes = np.asfortranarray([[4.0, 5.0, 4.0], [0.0, 0.0, 1.0]])
        surface2 = self._make_one(nodes, 1)
        intersections = surface1.intersect(surface2)
        self.assertEqual(intersections, [])

    def test_intersect_tangent_bbox(self):
        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        nodes = np.asfortranarray([[0.0, 0.0, -1.0], [0.0, 1.0, 1.0]])
        surface2 = self._make_one(nodes, 1)
        intersections = surface1.intersect(surface2)
        self.assertEqual(intersections, [])

    def test_intersect_first_contained(self):
        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        nodes = np.asfortranarray([[-1.0, 3.0, -1.0], [-1.0, -1.0, 3.0]])
        surface2 = self._make_one(nodes, 1)
        intersections = surface1.intersect(surface2)
        self.assertEqual(intersections, [surface1])

    def test_intersect_second_contained(self):
        nodes = np.asfortranarray([[-1.0, 3.0, -1.0], [-1.0, -1.0, 3.0]])
        surface1 = self._make_one(nodes, 1)
        surface2 = self._make_one(self.UNIT_TRIANGLE, 1)
        intersections = surface1.intersect(surface2)
        self.assertEqual(intersections, [surface2])

    def test_intersect_non_surface(self):
        surface = self._make_one(self.UNIT_TRIANGLE, 1)
        with self.assertRaises(TypeError):
            surface.intersect(object())

    def test_intersect_unsupported_dimension(self):
        surface1 = self._make_one(self.UNIT_TRIANGLE, 1)
        nodes2 = np.zeros((3, 3), order="F")
        surface2 = self._make_one(nodes2, 1)
        with self.assertRaises(NotImplementedError):
            surface1.intersect(surface2)
        with self.assertRaises(NotImplementedError):
            surface2.intersect(surface1)

    def test_elevate_linear(self):
        nodes = np.asfortranarray([[0.0, 2.0, -1.0], [0.0, 1.0, 2.0]])
        surface = self._make_one(nodes, 1)
        elevated = surface.elevate()
        expected = np.asfortranarray(
            [[0.0, 1.0, 2.0, -0.5, 0.5, -1.0], [0.0, 0.5, 1.0, 1.0, 1.5, 2.0]]
        )
        self.assertEqual(surface._degree, 1)
        self.assertEqual(elevated._degree, 2)
        self.assertEqual(elevated._nodes, expected)
        main_vals = surface.evaluate_cartesian_multi(self.REF_TRIANGLE3)
        sub_vals = elevated.evaluate_cartesian_multi(self.REF_TRIANGLE3)
        self.assertEqual(main_vals, sub_vals)

    def test_elevate_quadratic(self):
        klass = self._get_target_class()
        nodes = np.asfortranarray([[0.0, 6.0, 9.0, 0.0, 6.0, -3.0]])
        surface = klass.from_nodes(nodes)
        elevated = surface.elevate()
        expected = np.asfortranarray(
            [[0.0, 4.0, 7.0, 9.0, 0.0, 4.0, 7.0, -1.0, 3.0, -3.0]]
        )
        self.assertEqual(surface._degree, 2)
        self.assertEqual(elevated._degree, 3)
        self.assertEqual(elevated._nodes, expected)
        main_vals = surface.evaluate_cartesian_multi(self.REF_TRIANGLE3)
        sub_vals = elevated.evaluate_cartesian_multi(self.REF_TRIANGLE3)
        self.assertEqual(main_vals, sub_vals)


class Test__make_intersection(utils.NumPyTestCase):
    @staticmethod
    def _call_function_under_test(edge_info, all_edge_nodes):
        from bezier import surface

        return surface._make_intersection(edge_info, all_edge_nodes)

    def test_it(self):
        import bezier

        nodes1 = np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
        surface1 = bezier.Surface(nodes1, degree=1, copy=False)
        nodes2 = np.asfortranarray([[0.25, -0.75, 0.25], [0.25, 0.25, -0.75]])
        surface2 = bezier.Surface(nodes2, degree=1, copy=False)
        edge_nodes = tuple(edge._nodes for edge in surface1.edges) + tuple(
            edge._nodes for edge in surface2.edges
        )
        edge_info = (
            (0, 0.0, 0.25),
            (5, 0.75, 1.0),
            (3, 0.0, 0.25),
            (2, 0.75, 1.0),
        )
        result = self._call_function_under_test(edge_info, edge_nodes)
        self.assertIsInstance(result, bezier.CurvedPolygon)
        self.assertEqual(result._metadata, edge_info)
        self.assertEqual(result.num_sides, 4)
        # pylint: disable=unbalanced-tuple-unpacking
        edge0, edge1, edge2, edge3 = result._edges
        # pylint: enable=unbalanced-tuple-unpacking
        # First edge.
        expected = np.asfortranarray([[0.0, 0.25], [0.0, 0.0]])
        self.assertEqual(edge0._nodes, expected)
        # Second edge.
        expected = np.asfortranarray([[0.25, 0.25], [0.0, 0.25]])
        self.assertEqual(edge1._nodes, expected)
        # Third edge.
        expected = np.asfortranarray([[0.25, 0.0], [0.25, 0.25]])
        self.assertEqual(edge2._nodes, expected)
        # Fourth edge.
        expected = np.asfortranarray([[0.0, 0.0], [0.25, 0.0]])
        self.assertEqual(edge3._nodes, expected)
