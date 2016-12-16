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

import mock


class Test_add_plot_boundary(unittest.TestCase):

    AX_PROPS = ('lines', 'set_xlim', 'set_ylim')

    @staticmethod
    def _call_function_under_test(ax, **kwargs):
        from bezier import _plot_helpers

        return _plot_helpers.add_plot_boundary(ax, **kwargs)

    def test_default(self):
        import matplotlib.lines

        line = matplotlib.lines.Line2D([-1.0, 1.0], [-1.0, 1.0])
        ax = mock.Mock(name='AxesSubplot', spec=self.AX_PROPS, lines=[line])

        self.assertIsNone(self._call_function_under_test(ax))
        ax.set_xlim.assert_called_once_with(-1.125, 1.125)
        ax.set_ylim.assert_called_once_with(-1.125, 1.125)

    def test_with_padding(self):
        import matplotlib.lines

        line = matplotlib.lines.Line2D([-1.0, 1.0], [-1.0, 1.0])
        ax = mock.Mock(name='AxesSubplot', spec=self.AX_PROPS, lines=[line])

        self.assertIsNone(self._call_function_under_test(ax, padding=0.5))
        ax.set_xlim.assert_called_once_with(-1.5, 1.5)
        ax.set_ylim.assert_called_once_with(-1.5, 1.5)
