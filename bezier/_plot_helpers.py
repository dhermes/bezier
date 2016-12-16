# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

"""Plotting utilities."""


import numpy as np

from bezier import _helpers


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
    nodes = np.vstack([line.get_xydata() for line in ax.lines])
    left, right, bottom, top = _helpers.bbox(nodes)
    center_x = 0.5 * (right + left)
    delta_x = right - left
    center_y = 0.5 * (top + bottom)
    delta_y = top - bottom

    multiplier = (1.0 + padding) * 0.5
    ax.set_xlim(center_x - multiplier * delta_x,
                center_x + multiplier * delta_x)
    ax.set_ylim(center_y - multiplier * delta_y,
                center_y + multiplier * delta_y)
