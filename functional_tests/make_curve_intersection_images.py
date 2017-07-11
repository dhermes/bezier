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

from __future__ import absolute_import

import matplotlib.pyplot as plt

import bezier
from bezier import _plot_helpers

import runtime_utils


CONFIG = runtime_utils.Config()
CURVES, INTERSECTIONS = runtime_utils.get_intersections_info()
FILENAME_TEMPLATE = 'test_curves{:d}_and_{:d}'


def _get_curve(curve_id):
    curve_info = CURVES[curve_id - 1]
    assert curve_info['id'] == curve_id
    return bezier.Curve.from_nodes(
        curve_info['control_points'], _copy=False)


def make_plot(intersection_info, save_plot):
    curve_id1 = intersection_info['curve1']
    curve1 = _get_curve(curve_id1)
    curve_id2 = intersection_info['curve2']
    curve2 = _get_curve(curve_id2)
    intersection_pts = intersection_info['intersections']

    ax = curve1.plot(64)
    curve2.plot(64, ax=ax)
    ax.plot(intersection_pts[:, 0], intersection_pts[:, 1],
            marker='o', linestyle='None', color='black')
    ax.axis('scaled')
    _plot_helpers.add_plot_boundary(ax)

    filename = FILENAME_TEMPLATE.format(curve_id1, curve_id2)
    if save_plot:
        # NOTE: This is an abuse of ``current_test``, but we don't need
        #       a full-fledged ``Config``, we **just** need ``save_fig``.
        CONFIG.current_test = filename
        CONFIG.save_fig()
    else:
        plt.title(filename)
        plt.show()

    plt.close(ax.figure)


def main():
    parser = runtime_utils.get_parser()
    args = parser.parse_args()
    for intersection_info in INTERSECTIONS:
        make_plot(intersection_info, args.save_plot)


if __name__ == '__main__':
    main()
