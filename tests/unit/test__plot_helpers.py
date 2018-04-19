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
import functools
import sys
import unittest
import unittest.mock

import numpy as np

from tests.unit import utils


def run_fake_modules(modules, func):
    existing = {}
    for name, mod_obj in modules.items():
        if name in sys.modules:  # pragma: NO COVER
            existing[name] = sys.modules.pop(name)
        sys.modules[name] = mod_obj
    try:
        return func()

    finally:
        for name in modules.keys():
            sys.modules.pop(name)
            if name in existing:  # pragma: NO COVER
                sys.modules[name] = existing[name]


class Test_new_axis(unittest.TestCase):

    @staticmethod
    def _call_function_under_test():
        from bezier import _plot_helpers
        return _plot_helpers.new_axis()

    def test_it(self):
        figure = unittest.mock.Mock(spec=["gca"])
        figure.gca.return_value = unittest.mock.sentinel.ax
        plt = unittest.mock.Mock(spec=["figure"])
        plt.figure.return_value = figure
        matplotlib = unittest.mock.Mock(pyplot=plt, spec=["pyplot"])
        modules = {"matplotlib": matplotlib, "matplotlib.pyplot": plt}
        result = run_fake_modules(modules, self._call_function_under_test)
        self.assertIs(result, unittest.mock.sentinel.ax)
        # Verify mocks.
        plt.figure.assert_called_once_with()
        figure.gca.assert_called_once_with()


class Test_add_plot_boundary(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(ax, **kwargs):
        from bezier import _plot_helpers
        return _plot_helpers.add_plot_boundary(ax, **kwargs)

    def _helper(self, **kwargs):
        line = unittest.mock.Mock(spec=["get_xydata"])
        line.get_xydata.return_value = np.asfortranarray(
            [[-1.0, -1.0], [1.0, 1.0]]
        )
        ax = unittest.mock.Mock(
            lines=[line], spec=["lines", "set_xlim", "set_ylim"]
        )
        self.assertIsNone(self._call_function_under_test(ax, **kwargs))
        padding = kwargs.get("padding", 0.125)
        ax.set_xlim.assert_called_once_with(-1.0 - padding, 1.0 + padding)
        ax.set_ylim.assert_called_once_with(-1.0 - padding, 1.0 + padding)

    def test_default(self):
        self._helper()

    def test_with_padding(self):
        self._helper(padding=0.5)


class Test_add_patch(utils.NumPyTestCase):

    @staticmethod
    def _call_function_under_test(ax, color, pts_per_edge, *edges):
        from bezier import _plot_helpers
        return _plot_helpers.add_patch(ax, color, pts_per_edge, *edges)

    def _path_val(self, path, expected_transpose):
        self.assertEqual(path.Path.call_count, 1)
        call = path.Path.mock_calls[0]
        _, positional, keyword = call
        self.assertEqual(keyword, {})
        self.assertEqual(len(positional), 1)
        self.assertEqual(positional[0].T, expected_transpose)

    def _plot_check(self, ax, expected, color):
        self.assertEqual(ax.plot.call_count, 1)
        call = ax.plot.mock_calls[0]
        utils.check_plot_call(self, call, expected, color=color)

    @staticmethod
    def _get_info(nodes):
        from bezier import surface as surface_mod

        surface = surface_mod.Surface(nodes, 1, _copy=False)
        edges = surface._get_edges()
        expected = np.empty((2, 7), order="F")
        expected[:, 0] = 0.5 * (nodes[:, 0] + nodes[:, 1])
        expected[:, 1] = nodes[:, 1]
        expected[:, 2] = 0.5 * (nodes[:, 1] + nodes[:, 2])
        expected[:, 3] = nodes[:, 2]
        expected[:, 4] = 0.5 * (nodes[:, 2] + nodes[:, 0])
        expected[:, 5] = nodes[:, 0]
        expected[:, 6] = expected[:, 0]
        return expected, edges

    def test_it(self):
        # Set-up input values.
        color = (0.25, 0.5, 0.75)
        pts_per_edge = 3
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [1.0, 3.0, 6.0]])
        expected_transpose, edges = self._get_info(nodes)
        # Set-up mocks (quite a lot).
        patches = unittest.mock.Mock(spec=["PathPatch"])
        path = unittest.mock.Mock(spec=["Path"])
        matplotlib = unittest.mock.Mock(
            patches=patches, path=path, spec=["patches", "path"]
        )
        ax = unittest.mock.Mock(spec=["add_patch", "plot"])
        line = unittest.mock.Mock(spec=["get_color"])
        line.get_color.return_value = color
        ax.plot.return_value = (line,)
        # Run the code with fake modules.
        modules = {
            "matplotlib": matplotlib,
            "matplotlib.patches": patches,
            "matplotlib.path": path,
        }
        func = functools.partial(
            self._call_function_under_test, ax, color, pts_per_edge, *edges
        )
        result = run_fake_modules(modules, func)
        self.assertIsNone(result)
        # Verify mocks (quite a lot).
        self._path_val(path, expected_transpose)
        patches.PathPatch.assert_called_once_with(
            path.Path.return_value, facecolor=color, alpha=0.625
        )
        ax.add_patch.assert_called_once_with(patches.PathPatch.return_value)
        line.get_color.assert_called_once_with()
        self._plot_check(ax, expected_transpose, color)
