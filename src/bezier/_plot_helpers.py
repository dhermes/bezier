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
"""Plotting utilities."""

import numpy as np

from bezier import _helpers


def new_axis():
    """Get a new matplotlib axis.

    Returns:
        matplotlib.artist.Artist: A newly created axis.
    """
    # NOTE: We import the plotting library at runtime to
    #       avoid the cost for users that only want to compute.
    #       The ``matplotlib`` import is a tad expensive.
    import matplotlib.pyplot as plt

    figure = plt.figure()
    return figure.gca()


def add_plot_boundary(ax, padding=0.125):
    """Add a buffer of empty space around a plot boundary.

    .. note::

       This only uses ``line`` data from the axis. It **could**
       use ``patch`` data, but doesn't at this time.

    Args:
        ax (matplotlib.artist.Artist): A matplotlib axis.
        padding (Optional[float]): Amount (as a fraction of width and height)
            of padding to add around data. Defaults to ``0.125``.
    """
    nodes = np.asfortranarray(
        np.vstack([line.get_xydata() for line in ax.lines]).T
    )
    left, right, bottom, top = _helpers.bbox(nodes)
    center_x = 0.5 * (right + left)
    delta_x = right - left
    center_y = 0.5 * (top + bottom)
    delta_y = top - bottom
    multiplier = (1.0 + padding) * 0.5
    ax.set_xlim(
        center_x - multiplier * delta_x, center_x + multiplier * delta_x
    )
    ax.set_ylim(
        center_y - multiplier * delta_y, center_y + multiplier * delta_y
    )


def add_patch(ax, color, pts_per_edge, *edges):
    """Add a polygonal surface patch to a plot.

    Args:
        ax (matplotlib.artist.Artist): A matplotlib axis.
        color (Tuple[float, float, float]): Color as RGB profile.
        pts_per_edge (int): Number of points to use in polygonal
            approximation of edge.
        edges (Tuple[~bezier.curve.Curve, ...]): Curved edges defining
            a boundary.
    """
    from matplotlib import patches
    from matplotlib import path as _path_mod

    s_vals = np.linspace(0.0, 1.0, pts_per_edge)
    # Evaluate points on each edge.
    all_points = []
    for edge in edges:
        points = edge.evaluate_multi(s_vals)
        # We assume the edges overlap and leave out the first point
        # in each.
        all_points.append(points[:, 1:])
    # Add first point as last point (polygon is closed).
    first_edge = all_points[0]
    all_points.append(first_edge[:, [0]])
    # Add boundary first.
    polygon = np.asfortranarray(np.hstack(all_points))
    line, = ax.plot(polygon[0, :], polygon[1, :], color=color)
    # Reset ``color`` in case it was ``None`` and set from color wheel.
    color = line.get_color()
    # ``polygon`` is stored Fortran-contiguous with ``x-y`` points in each
    # column but ``Path()`` wants ``x-y`` points in each row.
    path = _path_mod.Path(polygon.T)
    patch = patches.PathPatch(path, facecolor=color, alpha=0.625)
    ax.add_patch(patch)
