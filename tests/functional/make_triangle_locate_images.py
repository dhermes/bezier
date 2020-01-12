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
from tests.functional import test_triangle_locate
from tests.functional import utils


def make_plot(triangle_index, point_index, save_plot):
    triangle = test_triangle_locate.TRIANGLES[triangle_index]
    point = test_triangle_locate.POINTS[:, [point_index]]
    name = f"test_triangle{triangle_index}_and_point{point_index}"

    ax = triangle.plot(64)
    ax.plot(
        point[0, :], point[1, :], color="black", marker="o", linestyle="None"
    )
    ax.axis("scaled")
    _plot_helpers.add_plot_boundary(ax)
    if save_plot:
        utils.save_fig(name)
    else:
        plt.title(name.replace("_", r"\_"))
        plt.show()
    plt.close(ax.figure)


def main():
    parser = utils.get_parser()
    args = parser.parse_args()
    for case in test_triangle_locate.CASES:
        triangle_index, point_index, _, _ = case
        make_plot(triangle_index, point_index, args.save_plot)


if __name__ == "__main__":
    seaborn.set()  # Required in `seaborn >= 0.8`
    main()
