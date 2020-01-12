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

import collections

import matplotlib.pyplot as plt
import seaborn

from bezier import _helpers
from tests.functional import test_segment_box
from tests.functional import utils


# As of ``0.9.0``, this palette has (BLUE, ORANGE, GREEN, RED, PURPLE, BROWN).
_COLORS = seaborn.color_palette(palette="deep", n_colors=6)
BLUE = _COLORS[0]
GREEN = _COLORS[2]
del _COLORS


def make_plot(name, segment, index, save_plot):
    figure = plt.figure()
    ax = figure.gca()
    (line,) = ax.plot([0, 1], [0, 1], alpha=0.0, color=BLUE)
    ax.fill_between([0, 1], [0, 0], [1, 1], alpha=0.5, color=line.get_color())
    (line,) = ax.plot(segment[0, :], segment[1, :], color=GREEN)
    ax.plot(
        segment[0, 0],
        segment[1, 0],
        marker="o",
        linestyle="None",
        color=line.get_color(),
    )
    left_, right_, bottom_, top_ = _helpers.bbox(segment)
    ax.fill_between(
        [left_, right_],
        [bottom_, bottom_],
        [top_, top_],
        alpha=0.5,
        color=line.get_color(),
    )
    ax.axis("scaled")
    ax.set_xlim(-1.125, 2.125)
    ax.set_ylim(-1.125, 2.125)
    if save_plot:
        utils.save_fig(f"test_{name}{index:02d}")
    else:
        plt.title(name.replace("_", r"\_"))
        plt.show()
    plt.close(figure)


def main():
    parser = utils.get_parser()
    args = parser.parse_args()

    index_by_name = collections.defaultdict(int)
    for value_triple in test_segment_box.SEGMENTS:
        name, segment, _ = value_triple
        curr_index = index_by_name[name]
        make_plot(name, segment, curr_index, args.save_plot)
        # For next occurence of ``name``.
        index_by_name[name] = curr_index + 1


if __name__ == "__main__":
    seaborn.set()  # Required in `seaborn >= 0.8`
    main()
