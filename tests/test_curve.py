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


class TestCurve(unittest.TestCase):

    @staticmethod
    def _get_target_class():
        from bezier import curve
        return curve.Curve

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        import numpy as np

        nodes = np.array([
            [0.0, 0.0],
            [0.625, 0.5],
            [1.0, 0.75],
        ])
        curve = self._make_one(nodes)
        self.assertEqual(curve._degree, 2)
        self.assertEqual(curve._dimension, 2)
        self.assertIs(curve._nodes, nodes)

    def test_constructor_bad_input(self):
        import numpy as np

        nodes = np.array([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

        nodes = np.zeros((2, 2, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test_degree_property(self):
        import numpy as np

        degree = 6
        num_nodes = degree + 1
        nodes = np.zeros((num_nodes, 2))
        curve = self._make_one(nodes)
        self.assertEqual(curve.degree, degree)

    def test_dimension_property(self):
        import numpy as np

        dimension = 4
        nodes = np.zeros((3, dimension))
        curve = self._make_one(nodes)
        self.assertEqual(curve.dimension, dimension)

    def test_evaluate(self):
        import numpy as np

        s = 0.25
        nodes = np.array([
            [0.0, 0.0],
            [0.5, 0.5],
            [1.0, 1.25],
        ])
        curve = self._make_one(nodes)
        expected = np.array([0.25,  0.265625])
        result = curve.evaluate(s)
        self.assertTrue(np.all(expected == result))

    def test_evaluate_degree_zero(self):
        import numpy as np

        nodes = np.array([
            [0.0, 0.0],
        ])
        curve = self._make_one(nodes)
        expected = nodes.flatten()
        result = curve.evaluate(0.0)
        self.assertTrue(np.all(expected == result))

    def test_evaluate_multi(self):
        import numpy as np

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
        import numpy as np

        s1 = 3.14159
        s2 = 2.817281728
        s_vals = np.array([s1, s2])
        num_pts = len(s_vals)
        curve = self._make_one(np.zeros((2, 1)))
        ret_vals = [10.0, -1.0]
        curve.evaluate = mock.Mock(side_effect=ret_vals)

        result = curve.evaluate_multi(s_vals)
        self.assertEqual(result.shape, (num_pts, 1))
        self.assertTrue(np.all(result == np.array([ret_vals]).T))

        curve.evaluate.assert_any_call(s1)
        curve.evaluate.assert_any_call(s2)
        self.assertEqual(curve.evaluate.call_count, 2)

    def _plot_helper(self, show=False):
        import mock
        import numpy as np

        nodes = np.array([
            [0.0, 0.0],
        ])
        curve = self._make_one(nodes)
        plt = mock.Mock()

        figure = mock.Mock()
        plt.figure.return_value = figure
        ax = mock.Mock()
        figure.gca.return_value = ax

        if show:
            result = curve.plot(1, plt, show=True)
        else:
            result = curve.plot(1, plt)

        # Check mocks.
        plt.figure.assert_called_once_with()
        figure.gca.assert_called_once_with()
        ax.plot.assert_called_once_with(nodes[:, 0], nodes[:, 1])
        if show:
            plt.show.assert_called_once_with()
        else:
            plt.show.assert_not_called()

    def test_plot(self):
        self._plot_helper()

    def test_plot_show(self):
        self._plot_helper(show=True)
