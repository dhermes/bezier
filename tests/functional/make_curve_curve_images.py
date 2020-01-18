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

import matplotlib.pyplot as plt
import seaborn

from bezier import _plot_helpers
from tests.functional import utils


_, INTERSECTIONS = utils.curve_intersections_info()
# As of ``0.9.0``, this palette has (BLUE, ORANGE, GREEN, RED, PURPLE, BROWN).
_COLORS = seaborn.color_palette(palette="deep", n_colors=6)
BLUE = _COLORS[0]
GREEN = _COLORS[2]
del _COLORS


def make_plot(intersection_info, save_plot):
    curve1 = intersection_info.curve1
    curve2 = intersection_info.curve2
    intersection_pts = intersection_info.intersections
    ax = curve1.plot(64, color=BLUE)
    curve2.plot(64, ax=ax, color=GREEN)
    ax.plot(
        intersection_pts[0, :],
        intersection_pts[1, :],
        marker="o",
        linestyle="None",
        color="black",
    )
    ax.axis("scaled")
    _plot_helpers.add_plot_boundary(ax)
    filename = intersection_info.img_filename
    if save_plot:
        utils.save_fig(filename)
    else:
        plt.title(filename)
        plt.show()
    plt.close(ax.figure)


def main():
    parser = utils.get_parser()
    args = parser.parse_args()
    for intersection_info in INTERSECTIONS:
        make_plot(intersection_info, args.save_plot)


if __name__ == "__main__":
    seaborn.set()  # Required in `seaborn >= 0.8`
    main()
