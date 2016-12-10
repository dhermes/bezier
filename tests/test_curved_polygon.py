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


class TestCurvedPolygon(utils.NumPyTestCase):

    NODES0 = np.array([
        [0.0, 0.0],
        [0.5, -1.0],
        [1.0, 0.0],
    ])
    NODES1 = np.array([
        [1.0, 0.0],
        [0.5, 1.0],
        [0.0, 0.0],
    ])
    COLOR = (0.125, 0.125, 0.0)

    @staticmethod
    def _get_target_class():
        from bezier import curved_polygon

        return curved_polygon.CurvedPolygon

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def _make_default(self):
        import bezier

        edge0 = bezier.Curve(self.NODES0)
        edge1 = bezier.Curve(self.NODES1)
        return self._make_one(edge0, edge1)

    def test_constructor(self):
        import bezier

        edge0 = bezier.Curve(self.NODES0)
        edge1 = bezier.Curve(self.NODES1)
        curved_poly = self._make_one(edge0, edge1)
        self.assertEqual(curved_poly._edges, (edge0, edge1))
        self.assertEqual(curved_poly._num_sides, 2)

    def test__verify_too_few(self):
        with self.assertRaises(ValueError):
            self._make_one()
        with self.assertRaises(ValueError):
            self._make_one(None)

    def test__verify_bad_dimension(self):
        import bezier

        edge0 = bezier.Curve(np.array([
            [1.0, 1.0],
            [2.0, 2.0],
        ]))
        edge1 = bezier.Curve(self.NODES1)
        with self.assertRaises(ValueError):
            self._make_one(edge0, edge1)

    def test__verify_not_aligned(self):
        import bezier

        edge0 = bezier.Curve(np.array([[0.0], [0.0]]))
        edge1 = bezier.Curve(self.NODES1)
        with self.assertRaises(ValueError):
            self._make_one(edge0, edge1)

    def test_num_sides_property(self):
        curved_poly = self._make_default()
        self.assertIs(curved_poly.num_sides, 2)

    def test___repr__(self):
        curved_poly = self._make_default()
        self.assertEqual(repr(curved_poly),
                         '<CurvedPolygon (num_sides=2)>')

    def _check_plot_calls(self, ax):
        # Check the calls to ax.plot(). We can't assert_any_call()
        # since == breaks on NumPy arrays.
        self.assertEqual(ax.plot.call_count, 2)

        calls = ax.plot.mock_calls
        utils.check_plot_call(
            self, calls[0], self.NODES0[(0, 2), :], color=None)
        utils.check_plot_call(
            self, calls[1], self.NODES1[(0, 2), :], color=self.COLOR)

    def _make_axis(self):
        import matplotlib.lines

        ax = mock.Mock()

        line = matplotlib.lines.Line2D([], [], color=self.COLOR)
        ax.plot.return_value = (line,)

        return ax

    def _plot_helper(self, show=False, ax=None):
        if ax is None:
            orig_ax = False
            ax = self._make_axis()
        else:
            orig_ax = True

        curved_poly = self._make_default()
        plt = mock.Mock()

        figure = mock.Mock()
        plt.figure.return_value = figure
        figure.gca.return_value = ax

        with mock.patch('bezier.curved_polygon.plt', new=plt):
            kwargs = {}
            if show:
                kwargs['show'] = True
            if orig_ax:
                kwargs['ax'] = ax
            result = curved_poly.plot(2, **kwargs)

        self.assertIs(result, ax)

        # Check mocks.
        if orig_ax:
            plt.figure.assert_not_called()
            figure.gca.assert_not_called()
        else:
            plt.figure.assert_called_once_with()
            figure.gca.assert_called_once_with()

        ax.add_patch.assert_called_once()

        self._check_plot_calls(ax)

        if show:
            plt.show.assert_called_once_with()
        else:
            plt.show.assert_not_called()

    def test_plot(self):
        self._plot_helper()

    def test_plot_show(self):
        self._plot_helper(show=True)

    def test_plot_existing_axis(self):
        ax = self._make_axis()
        self._plot_helper(ax=ax)
