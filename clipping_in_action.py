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

from matplotlib import patches
from matplotlib import path as _path_mod
import matplotlib.pyplot as plt
import numpy as np
import seaborn

from bezier import _curve_helpers
from bezier.hazmat import clipping


S_VALS = np.linspace(0.0, 1.0, 1024)
ALPHA = 0.5


def reduce_and_plot(
    nodes1, nodes2, original_nodes1=None, original_nodes2=None
):
    start, end = clipping.clip_range(nodes1, nodes2)
    new_nodes2 = _curve_helpers.specialize_curve(nodes2, start, end)

    figure = plt.figure()
    ax = figure.gca()

    if original_nodes1 is None:
        color1 = None
        marker1 = None
    else:
        points = _curve_helpers.evaluate_multi(original_nodes1, S_VALS)
        (line,) = ax.plot(points[0, :], points[1, :], alpha=ALPHA)
        color1 = line.get_color()
        marker1 = "o"

    if original_nodes2 is None:
        color2 = None
        marker2 = None
    else:
        points = _curve_helpers.evaluate_multi(original_nodes2, S_VALS)
        (line,) = ax.plot(points[0, :], points[1, :], alpha=ALPHA)
        color2 = line.get_color()
        marker2 = "o"

    a, b, _, d_min, d_max = clipping.compute_fat_line(nodes1)
    ax.plot(
        nodes1[0, (0, -1)],
        nodes1[1, (0, -1)],
        color="black",
        linestyle="dashed",
    )

    polygon = np.asfortranarray(
        [
            [nodes1[0, 0] + d_min * a, nodes1[1, 0] + d_min * b],
            [nodes1[0, -1] + d_min * a, nodes1[1, -1] + d_min * b],
            [nodes1[0, -1] + d_max * a, nodes1[1, -1] + d_max * b],
            [nodes1[0, 0] + d_max * a, nodes1[1, 0] + d_max * b],
            [nodes1[0, 0] + d_min * a, nodes1[1, 0] + d_min * b],
        ]
    )
    patch = patches.PathPatch(
        _path_mod.Path(polygon), facecolor="gray", alpha=ALPHA
    )
    ax.add_patch(patch)

    to_plot = (
        (nodes1, color1, marker1),
        (nodes2, color2, marker2),
        (new_nodes2, "black", None),
    )
    for nodes, color, marker in to_plot:
        points = _curve_helpers.evaluate_multi(nodes, S_VALS)
        ax.plot(points[0, :], points[1, :], color=color)
        if marker is not None:
            ax.plot(
                points[0, (0, -1)],
                points[1, (0, -1)],
                color=color,
                marker="o",
                linestyle="None",
            )

    ax.axis("scaled")
    plt.show()

    return start, end, new_nodes2


def main():
    original_nodes1 = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 2.0, 0.0]])
    original_nodes2 = np.asfortranarray([[0.0, 2.0, 0.0], [-1.0, 1.0, 3.0]])

    s1, t1, s2, t2 = 0.0, 1.0, 0.0, 1.0
    nodes1 = original_nodes1.copy(order="F")
    nodes2 = original_nodes2.copy(order="F")

    # Clip B2(t) on B1(s).
    start, end, nodes2 = reduce_and_plot(nodes1, nodes2)
    s2, t2 = s2 + start * (t2 - s2), s2 + end * (t2 - s2)
    print((s1, t1, s2, t2))

    # Clip B1(s) on (the 1-stage) B2(t).
    start, end, nodes1 = reduce_and_plot(
        nodes2, nodes1, original_nodes1=original_nodes2
    )
    s1, t1 = s1 + start * (t1 - s1), s1 + end * (t1 - s1)
    print((s1, t1, s2, t2))

    # Clip (the 1-stage) B2(t) on (the 1-stage) B1(s).
    start, end, nodes2 = reduce_and_plot(
        nodes1,
        nodes2,
        original_nodes1=original_nodes1,
        original_nodes2=original_nodes2,
    )
    s2, t2 = s2 + start * (t2 - s2), s2 + end * (t2 - s2)
    print((s1, t1, s2, t2))

    # Clip (the 1-stage) B1(s) on (the 2-stage) B2(t).
    start, end, nodes1 = reduce_and_plot(
        nodes2,
        nodes1,
        original_nodes1=original_nodes2,
        original_nodes2=original_nodes1,
    )
    s1, t1 = s1 + start * (t1 - s1), s1 + end * (t1 - s1)
    print((s1, t1, s2, t2))

    # Clip (the 2-stage) B2(t) on (the 2-stage) B1(s).
    start, end, nodes2 = reduce_and_plot(
        nodes1,
        nodes2,
        original_nodes1=original_nodes1,
        original_nodes2=original_nodes2,
    )
    s2, t2 = s2 + start * (t2 - s2), s2 + end * (t2 - s2)
    print((s1, t1, s2, t2))

    # Clip (the 2-stage) B1(s) on (the 3-stage) B2(t).
    start, end, nodes1 = reduce_and_plot(
        nodes2,
        nodes1,
        original_nodes1=original_nodes2,
        original_nodes2=original_nodes1,
    )
    s1, t1 = s1 + start * (t1 - s1), s1 + end * (t1 - s1)
    print((s1, t1, s2, t2))


if __name__ == "__main__":
    seaborn.set()
    main()
