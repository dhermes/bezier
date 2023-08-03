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

try:
    import sympy
except ImportError:  # pragma: NO COVER
    sympy = None

from tests.unit import test__symbolic
from tests.unit import utils


SPACING = np.spacing  # pylint: disable=no-member


class TestCurve(utils.NumPyTestCase):
    ZEROS = np.zeros((2, 2), order="F")

    @staticmethod
    def _get_target_class():
        from bezier import curve

        return curve.Curve

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        nodes = np.asfortranarray([[0.0, 0.625, 1.0], [0.0, 0.5, 0.75]])
        curve = self._make_one(nodes, 2, copy=False)
        self.assertEqual(curve._degree, 2)
        self.assertEqual(curve._dimension, 2)
        self.assertIs(curve._nodes, nodes)

    def test_constructor_invalid_num_nodes(self):
        nodes = np.empty((1, 3), order="F")
        with self.assertRaises(ValueError) as exc_info:
            self._make_one(nodes, 7, copy=False)

        exc_args = exc_info.exception.args
        self.assertEqual(
            exc_args, ("A degree 7 curve should have 8 nodes, not 3.",)
        )

    def test_constructor_wrong_dimension(self):
        nodes = np.asfortranarray([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes, None)
        nodes = np.zeros((2, 2, 2), order="F")
        with self.assertRaises(ValueError):
            self._make_one(nodes, None)

    def test_from_nodes_factory(self):
        nodes = np.asfortranarray([[1.0, 1.0], [2.0, 3.0], [0.0, 0.0]])
        klass = self._get_target_class()
        curve = klass.from_nodes(nodes)
        self.assertIsInstance(curve, klass)
        self.assertEqual(curve._degree, 1)
        self.assertEqual(curve._dimension, 3)
        self.assertEqual(curve._nodes, nodes)

    def test_from_nodes_factory_non_array(self):
        nodes = [[1.0, 1.0, 2.0], [2.0, 3.0, 4.0]]
        klass = self._get_target_class()
        curve = klass.from_nodes(nodes)
        self.assertIsInstance(curve, klass)
        self.assertEqual(curve._degree, 2)
        self.assertEqual(curve._dimension, 2)
        self.assertTrue(np.all(curve._nodes == nodes))

    def test__get_degree(self):
        klass = self._get_target_class()
        self.assertEqual(0, klass._get_degree(1))
        self.assertEqual(1, klass._get_degree(2))

    def test_length_property(self):
        nodes = np.asfortranarray([[0.0, 3.0], [0.0, 4.0]])
        curve = self._make_one(nodes, 1)
        self.assertEqual(curve.length, 5.0)

    def test___dict___property(self):
        curve = self._make_one(self.ZEROS, 1, copy=False)
        props_dict = curve.__dict__
        expected = {"_nodes": self.ZEROS, "_dimension": 2, "_degree": 1}
        self.assertEqual(props_dict, expected)
        # Check that modifying ``props_dict`` won't modify ``curve``.
        expected["_dimension"] = 47
        self.assertNotEqual(curve._dimension, expected["_dimension"])

    def test_copy(self):
        nodes = np.asfortranarray([[2.0, 3.5, 4.0], [0.0, 1.0, 0.0]])
        curve = self._make_one(nodes, 2)

        new_curve = curve.copy()
        self.assertEqual(curve._degree, new_curve._degree)
        self.assertEqual(curve._dimension, new_curve._dimension)
        self.assertTrue(np.all(curve._nodes == new_curve._nodes))
        self.assertIsNot(curve._nodes, new_curve._nodes)

    def test_evaluate(self):
        s = 0.25
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 0.5, 1.25]])
        curve = self._make_one(nodes, 2)
        expected = np.asfortranarray([[0.25], [0.265625]])
        result = curve.evaluate(s)
        self.assertEqual(expected, result)

    def test_evaluate_multi(self):
        s_vals = np.asfortranarray([0.0, 0.25, 0.5, 1.0, 1.25])
        nodes = np.asfortranarray([[0.0, 0.375, 1.0], [0.0, 0.375, 1.0]])
        curve = self._make_one(nodes, 2)
        expected = np.asfortranarray(
            [
                [0.0, 0.203125, 0.4375, 1.0, 1.328125],
                [0.0, 0.203125, 0.4375, 1.0, 1.328125],
            ]
        )
        result = curve.evaluate_multi(s_vals)
        self.assertEqual(expected, result)

    def test_evaluate_hodograph(self):
        s = 0.25
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 0.5, 1.25]])
        curve = self._make_one(nodes, 2)
        expected = np.asfortranarray([[1.0], [1.125]])
        result = curve.evaluate_hodograph(s)
        self.assertEqual(expected, result)

    def test_plot_wrong_dimension(self):
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 3.0], [0.0, 4.0]])
        curve = self._make_one(nodes, 1)
        with self.assertRaises(NotImplementedError):
            curve.plot(32)

    @unittest.mock.patch("bezier._plot_helpers.new_axis")
    def test_plot_defaults(self, new_axis_mock):
        ax = unittest.mock.Mock(spec=["plot"])
        new_axis_mock.return_value = ax
        nodes = np.asfortranarray([[0.0, 1.0], [1.0, 3.0]])
        curve = self._make_one(nodes, 1, copy=False)
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

    @unittest.mock.patch("bezier._plot_helpers.new_axis")
    def test_plot_explicit(self, new_axis_mock):
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve = self._make_one(nodes, 1, copy=False)
        num_pts = 2  # This value is crucial for the plot call.
        ax = unittest.mock.Mock(spec=["plot"])
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
        nodes = np.asfortranarray([[0.0, 4.0], [1.0, 6.0]])
        klass = self._get_target_class()
        curve = klass.from_nodes(nodes)
        # Call ``subdivide()`` and then compare.
        left, right = curve.subdivide()
        # Check the "left" sub-curve.
        self.assertEqual(left._degree, 1)
        self.assertIsInstance(left, klass)
        expected_l = np.asfortranarray([[0.0, 2.0], [1.0, 3.5]])
        self.assertEqual(left._nodes, expected_l)
        # Check the "right" sub-curve.
        self.assertIsInstance(right, klass)
        expected_r = np.asfortranarray([[2.0, 4.0], [3.5, 6.0]])
        self.assertEqual(right._nodes, expected_r)

    def test_intersect_bad_strategy(self):
        curve = self._make_one(self.ZEROS, 1)
        strategy = unittest.mock.sentinel.bad_strategy
        with self.assertRaises(ValueError) as exc_info:
            curve.intersect(curve, strategy=strategy)
        exc_args = exc_info.exception.args
        self.assertEqual(exc_args, ("Unexpected strategy.", strategy))

    def test_intersect_algebraic(self):
        from bezier.hazmat import intersection_helpers

        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = self._make_one(nodes1, 1)
        nodes2 = np.asfortranarray([[0.0, 1.0], [1.0, 0.0]])
        curve2 = self._make_one(nodes2, 1)
        strategy = intersection_helpers.IntersectionStrategy.ALGEBRAIC
        intersections = curve1.intersect(curve2, strategy=strategy)
        expected = np.asfortranarray([[0.5], [0.5]])
        self.assertEqual(intersections, expected)

    def test_intersect_empty(self):
        nodes1 = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve1 = self._make_one(nodes1, 1)
        nodes2 = np.asfortranarray([[3.0, 2.0], [0.0, 1.0]])
        curve2 = self._make_one(nodes2, 1)
        result = curve1.intersect(curve2)
        self.assertEqual(result.shape, (2, 0))

    def test_intersect_at_boundary(self):
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, -0.25, 0.0]])
        curve = self._make_one(nodes, 2)
        left, right = curve.subdivide()
        result = left.intersect(right)
        # The "right" end of ``left`` and the "left" end of ``right``.
        expected = np.asfortranarray([[1.0], [0.0]])
        self.assertEqual(result, expected)

    def _intersect_helper(self, **kwargs):
        # NOTE: ``nodes1`` is a specialization of [0, 0], [1/2, 1], [1, 1]
        #       onto the interval [1/4, 1] and ``nodes`` is a specialization
        #       of [0, 1], [1/2, 1], [1, 0] onto the interval [0, 3/4].
        #       We expect them to intersect at s = 1/3, t = 2/3.
        nodes_left = np.asfortranarray(
            [[0.25, 0.625, 1.0], [0.4375, 1.0, 1.0]]
        )
        left = self._make_one(nodes_left, 2)
        nodes_right = np.asfortranarray(
            [[0.0, 0.375, 0.75], [1.0, 1.0, 0.4375]]
        )
        right = self._make_one(nodes_right, 2)
        result = left.intersect(right, **kwargs)
        expected = np.asfortranarray([[1.0], [2.0]]) / 3.0
        self.assertTrue(
            np.allclose(result, expected, atol=0.0, rtol=0.5**52)
        )

    def test_intersect(self):
        self._intersect_helper()

    def test_intersect_no_verify(self):
        self._intersect_helper(verify=False)

    def test_intersect_non_curve(self):
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, -0.25, 0.0]])
        curve = self._make_one(nodes, 2)
        with self.assertRaises(TypeError):
            curve.intersect(object())

    def test_intersect_unsupported_dimension(self):
        nodes = np.asfortranarray(
            [[0.0, 0.5, 1.0], [0.0, -0.25, 0.0], [0.0, 0.75, 1.25]]
        )
        curve1 = self._make_one(nodes, 2)
        curve2 = self._make_one(nodes[:, :2], 1)
        with self.assertRaises(NotImplementedError):
            curve1.intersect(curve2)
        with self.assertRaises(NotImplementedError):
            curve2.intersect(curve1)

    def test_self_intersections_unsupported_dimension(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0]])
        curve = self._make_one(nodes, 2)

        with self.assertRaises(NotImplementedError) as exc_info:
            curve.self_intersections()

        expected = (
            "Self-intersection only implemented in 2D",
            "Current dimension",
            1,
        )
        self.assertEqual(expected, exc_info.exception.args)

    def test_self_intersections_unsupported_strategy(self):
        from bezier.hazmat import intersection_helpers

        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 2.0]])
        curve = self._make_one(nodes, 2)

        strategy = intersection_helpers.IntersectionStrategy.ALGEBRAIC
        with self.assertRaises(NotImplementedError) as exc_info:
            curve.self_intersections(strategy=strategy)

        expected = ("Only geometric strategy for self-intersection detection",)
        self.assertEqual(expected, exc_info.exception.args)

    def test_self_intersections_valid(self):
        nodes = np.asfortranarray(
            [
                [0.0, -1.0, 1.0, -0.75],
                [2.0, 0.0, 1.0, 1.625],
            ]
        )
        curve = self._make_one(nodes, 3)
        intersections = curve.self_intersections(verify=False)
        self.assertEqual((2, 1), intersections.shape)
        intersections_verified = curve.self_intersections(verify=True)
        self.assertEqual(intersections, intersections_verified)

        sq5 = np.sqrt(5.0)
        expected_s = 0.5 - sq5 / 6.0
        local_eps = 3 * abs(SPACING(expected_s))
        self.assertAlmostEqual(
            expected_s, intersections[0, 0], delta=local_eps
        )
        expected_t = 0.5 + sq5 / 6.0
        local_eps = abs(SPACING(expected_t))
        self.assertAlmostEqual(
            expected_t, intersections[1, 0], delta=local_eps
        )

    def test_elevate(self):
        nodes = np.asfortranarray([[0.0, 1.0, 3.0, 3.5], [0.5, 1.0, 2.0, 4.0]])
        curve = self._make_one(nodes, 3)
        self.assertEqual(curve.degree, 3)
        elevated = curve.elevate()
        self.assertEqual(elevated.degree, 4)
        s_vals = np.linspace(0.0, 1.0, 64 + 1)
        orig_vals = curve.evaluate_multi(s_vals)
        new_vals = elevated.evaluate_multi(s_vals)
        self.assertEqual(orig_vals, new_vals)

    def test_reduce_(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0, 3.0], [0.0, 3.0, 3.0, 0.0]])
        curve = self._make_one(nodes, 3)
        self.assertEqual(curve.degree, 3)
        reduced = curve.reduce_()
        expected = np.asfortranarray([[0.0, 1.5, 3.0], [0.0, 4.5, 0.0]])
        self.assertEqual(reduced.nodes, expected)
        s_vals = np.linspace(0.0, 1.0, 64 + 1)
        orig_vals = curve.evaluate_multi(s_vals)
        new_vals = reduced.evaluate_multi(s_vals)
        self.assertEqual(orig_vals, new_vals)

    def test_specialize(self):
        nodes = np.asfortranarray([[0.0, 1.0, 5.0], [0.0, 6.0, 2.0]])
        curve = self._make_one(nodes, 2)
        new_curve = curve.specialize(0.25, 0.875)
        expected = np.asfortranarray(
            [[0.6875, 1.78125, 4.046875], [2.375, 4.5625, 2.84375]]
        )
        self.assertEqual(new_curve.nodes, expected)

    def test_locate_wrong_shape(self):
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
        curve = self._make_one(nodes, 1)
        point = np.asfortranarray([0.0, 1.0, 2.0])
        with self.assertRaises(ValueError):
            curve.locate(point)

    def test_locate(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0, 2.0, 5.0], [0.0, 1.0, -1.0, 1.0]]
        )
        curve = self._make_one(nodes, 3)
        s_val = 0.75
        point = curve.evaluate(s_val)
        result = curve.locate(point)
        self.assertEqual(result, s_val)

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_to_symbolic(self):
        nodes = np.asfortranarray([[3, 3, 4, 6], [3, 3, 3, 0]])
        curve = self._make_one(nodes, 3, copy=False)
        b_polynomial = curve.to_symbolic()

        s = sympy.Symbol("s")
        expected = 3 * sympy.Matrix([[1 + s**2, 1 - s**3]]).T
        self.assertTrue(
            test__symbolic.sympy_matrix_equal(b_polynomial, expected)
        )

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_implicitize(self):
        nodes = np.asfortranarray([[3, 3, 4, 6], [3, 3, 3, 0]])
        curve = self._make_one(nodes, 3, copy=False)
        f_polynomial = curve.implicitize()

        x_sym, y_sym = sympy.symbols("x, y")
        expected = -9 * (
            x_sym**3
            - 9 * x_sym**2
            + 27 * x_sym
            - 3 * y_sym**2
            + 18 * y_sym
            - 54
        )
        self.assertTrue(test__symbolic.sympy_equal(f_polynomial, expected))

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_implicitize_bad_dimension(self):
        nodes = np.empty((1, 2), order="F")
        curve = self._make_one(nodes, 1, copy=False)
        with self.assertRaises(ValueError) as exc_info:
            curve.implicitize()

        exc_args = exc_info.exception.args
        self.assertEqual(
            exc_args,
            (
                "Only a planar (2D) curve can be implicitized",
                "Current dimension",
                1,
            ),
        )
