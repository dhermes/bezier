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

from bezier import _plot_helpers

import runtime_utils


CONFIG = runtime_utils.Config()
_, INTERSECTIONS = runtime_utils.curve_intersections_info()


def make_plot(intersection_info, save_plot):
    curve1 = intersection_info.curve1
    curve2 = intersection_info.curve2
    intersection_pts = intersection_info.intersections

    ax = curve1.plot(64)
    curve2.plot(64, ax=ax)
    ax.plot(intersection_pts[:, 0], intersection_pts[:, 1],
            marker='o', linestyle='None', color='black')
    ax.axis('scaled')
    _plot_helpers.add_plot_boundary(ax)

    filename = intersection_info.img_filename
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
