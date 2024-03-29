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

import bezier
from bezier import _plot_helpers
from tests.functional import utils


_, INTERSECTIONS = utils.triangle_intersections_info()
# As of ``0.9.0``, this palette has (BLUE, ORANGE, GREEN, RED, PURPLE, BROWN).
_COLORS = seaborn.color_palette(palette="deep", n_colors=6)
BLUE = _COLORS[0]
GREEN = _COLORS[2]
RED = _COLORS[3]
del _COLORS


def make_curved_polygon(triangle1, triangle2, curved_polygon_info):
    if isinstance(curved_polygon_info, utils.CurvedPolygonInfo):
        base_edges = (triangle1._get_edges(), triangle2._get_edges())
        edge_list = curved_polygon_info.edge_list
        start_params = curved_polygon_info.start_params
        end_params = curved_polygon_info.end_params
        info = zip(edge_list, start_params, end_params)
        edges = []
        for edge_index, start_param, end_param in info:
            if edge_index < 3:
                base_edge = base_edges[0][edge_index]
            else:
                base_edge = base_edges[1][edge_index - 3]
            edges.append(base_edge.specialize(start_param, end_param))
        return bezier.CurvedPolygon(*edges)

    else:
        assert isinstance(curved_polygon_info, utils.TriangleIntersectionInfo)
        if curved_polygon_info.first:
            return triangle1

        else:
            return triangle2


def make_plot(intersection_info, save_plot):
    triangle1 = intersection_info.triangle1
    triangle2 = intersection_info.triangle2
    ax = triangle1.plot(64, color=BLUE)
    triangle2.plot(64, ax=ax, color=GREEN)
    color = RED
    for curved_polygon_info in intersection_info.intersections:
        curved_polygon = make_curved_polygon(
            triangle1, triangle2, curved_polygon_info
        )
        curved_polygon.plot(64, color=color, ax=ax)
        # Color is (R,G,B,A) but we just want (R,G,B).
        color = ax.patches[-1].get_facecolor()[:3]
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
    seaborn.set()
    main()
