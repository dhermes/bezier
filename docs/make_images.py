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

"""Helper to make images that are intended for docs."""


import os

import matplotlib.pyplot as plt
import numpy as np
try:
    import seaborn  # pylint: disable=unused-import
except ImportError:
    pass

import bezier


_DOCS_DIR = os.path.abspath(os.path.dirname(__file__))
IMAGES_DIR = os.path.join(_DOCS_DIR, 'images')


def save_image(figure, filename):
    """Save an image to the docs images directory.

    Args:
        filename (str): The name of the file (not containing
            directory info).
    """
    path = os.path.join(IMAGES_DIR, filename)
    figure.savefig(path, bbox_inches='tight')
    print('Saved {}'.format(filename))
    plt.close(figure)


def linearization_error():
    """Image for :func:`.linearization_error` docstring."""
    curve = bezier.Curve(np.array([
        [0.0, 0.0],
        [3.0, 1.0],
        [9.0, -2.0],
    ]))
    line = bezier.Curve(curve._nodes[(0, -1), :])

    midpoints = np.vstack([
        curve.evaluate(0.5),
        line.evaluate(0.5),
    ])

    ax = curve.plot(256)
    line.plot(256, ax=ax)
    ax.plot(midpoints[:, 0], midpoints[:, 1],
            color='black', linestyle='dashed')

    ax.axis('scaled')
    save_image(ax.figure, 'linearization_error.png')


def main():
    linearization_error()


if __name__ == '__main__':
    main()
