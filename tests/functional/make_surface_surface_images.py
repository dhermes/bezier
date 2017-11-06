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

from __future__ import absolute_import

import matplotlib.pyplot as plt
import six

import bezier
from bezier import _plot_helpers
from tests.functional import utils


CONFIG = utils.Config()
_, INTERSECTIONS = utils.surface_intersections_info()


def make_curved_polygon(surface1, surface2, curved_polygon_info):
    base_edges = (
        surface1._get_edges(),
        surface2._get_edges(),
    )

    edge_pairs = curved_polygon_info.edge_pairs
    start_params = curved_polygon_info.start_params
    end_params = curved_polygon_info.end_params
    info = six.moves.zip(edge_pairs, start_params, end_params)
    edges = []
    for edge_pair, start_param, end_param in info:
        surf_index, edge_index = edge_pair
        base_edge = base_edges[surf_index][edge_index]
        edges.append(base_edge.specialize(start_param, end_param))

    return bezier.CurvedPolygon(*edges)


def make_plot(intersection_info, save_plot):
    surface1 = intersection_info.surface1
    surface2 = intersection_info.surface2

    ax = surface1.plot(64)
    surface2.plot(64, ax=ax)

    color = None
    for curved_polygon_info in intersection_info.intersections:
        curved_polygon = make_curved_polygon(
            surface1, surface2, curved_polygon_info)
        curved_polygon.plot(64, color=color, ax=ax)
        # Color is (R,G,B,A) but we just want (R,G,B).
        color = ax.patches[-1].get_facecolor()[:3]

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
    parser = utils.get_parser()
    args = parser.parse_args()
    for intersection_info in INTERSECTIONS:
        make_plot(intersection_info, args.save_plot)


if __name__ == '__main__':
    main()
