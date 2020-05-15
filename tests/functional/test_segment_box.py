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

import numpy as np
import pytest

from bezier.hazmat import geometric_intersection


# Always gives us the unit square.
UNIT_SQUARE = np.asfortranarray([[0.0, 1.0], [0.0, 1.0]])
INTERSECTION = geometric_intersection.BoxIntersectionType.INTERSECTION
DISJOINT = geometric_intersection.BoxIntersectionType.DISJOINT
SEGMENTS = (
    ("outside", np.asfortranarray([[0.25, -0.75], [-0.75, 0.25]]), DISJOINT),
    ("outside", np.asfortranarray([[0.75, 1.75], [-0.75, 0.25]]), DISJOINT),
    ("outside", np.asfortranarray([[0.75, 1.75], [1.75, 0.75]]), DISJOINT),
    ("outside", np.asfortranarray([[-0.75, 0.25], [0.75, 1.75]]), DISJOINT),
    (
        "start_in_box",
        np.asfortranarray([[0.5, 0.5], [0.0, -0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[0.5, 0.5], [0.0, 0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[1.0, 1.5], [0.5, 0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[1.0, 0.5], [0.5, 0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[0.5, 0.5], [1.0, 1.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[0.5, 0.5], [1.0, 0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[0.0, -0.5], [0.5, 0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[0.0, 0.5], [0.5, 0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[0.5, 1.5], [0.5, 0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[0.5, 0.5], [0.5, -0.5]]),
        INTERSECTION,
    ),
    (
        "start_in_box",
        np.asfortranarray([[0.5, -0.5], [0.5, 1.5]]),
        INTERSECTION,
    ),
    ("end_in_box", np.asfortranarray([[0.5, 0.5], [-1.0, 0.0]]), INTERSECTION),
    (
        "end_in_box",
        np.asfortranarray([[0.25, 0.5], [-1.0, 0.0]]),
        INTERSECTION,
    ),
    ("end_in_box", np.asfortranarray([[1.5, 1.0], [0.5, 0.5]]), INTERSECTION),
    ("end_in_box", np.asfortranarray([[1.0, 1.0], [-0.5, 0.5]]), INTERSECTION),
    ("end_in_box", np.asfortranarray([[0.5, 0.5], [1.5, 1.0]]), INTERSECTION),
    ("end_in_box", np.asfortranarray([[0.5, 0.5], [-0.5, 1.0]]), INTERSECTION),
    ("end_in_box", np.asfortranarray([[-0.5, 0.0], [0.5, 0.5]]), INTERSECTION),
    ("end_in_box", np.asfortranarray([[0.5, 0.5], [1.5, 0.5]]), INTERSECTION),
    ("end_in_box", np.asfortranarray([[-0.5, 0.5], [0.5, 0.5]]), INTERSECTION),
    ("end_in_box", np.asfortranarray([[1.5, 0.5], [1.5, 0.5]]), INTERSECTION),
    (
        "goes_through_box",
        np.asfortranarray([[0.5, 1.25], [-0.25, 0.5]]),
        INTERSECTION,
    ),
    (
        "goes_through_box",
        np.asfortranarray([[-0.25, 0.5], [0.5, -0.25]]),
        INTERSECTION,
    ),
    (
        "goes_through_box",
        np.asfortranarray([[1.25, 0.5], [0.5, 1.25]]),
        INTERSECTION,
    ),
    (
        "goes_through_box",
        np.asfortranarray([[-0.5, 1.5], [0.5, 0.5]]),
        INTERSECTION,
    ),
    (
        "goes_through_box",
        np.asfortranarray([[-0.25, 1.25], [0.0, 0.0]]),
        INTERSECTION,
    ),
    (
        "goes_through_box",
        np.asfortranarray([[1.0, 1.0], [-0.25, 1.25]]),
        INTERSECTION,
    ),
    (
        "goes_through_box",
        np.asfortranarray([[1.25, -0.25], [1.0, 1.0]]),
        INTERSECTION,
    ),
    (
        "goes_through_box",
        np.asfortranarray([[0.0, 0.0], [1.25, -0.25]]),
        INTERSECTION,
    ),
    (
        "goes_through_box",
        np.asfortranarray([[0.5, -0.25], [1.25, 0.5]]),
        INTERSECTION,
    ),
)


def id_func(value_triple):
    name, segment, expected = value_triple
    segment_flat = tuple(segment.flatten())
    return f"{name}-{segment_flat}-{expected}"


@pytest.mark.parametrize("value_triple", SEGMENTS, ids=id_func)
def test_intersect(value_triple):
    _, segment, expected = value_triple
    result = geometric_intersection.bbox_line_intersect(
        UNIT_SQUARE, segment[:, 0], segment[:, 1]
    )
    assert result == expected
