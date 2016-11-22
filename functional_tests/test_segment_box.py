# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

import numpy as np
import pytest

import bezier
from bezier import _intersection_helpers


# Always gives us the unit square.
UNIT_SQUARE = np.array([
    [0.0, 0.0],
    [1.0, 1.0],
])


def run_it(segment, expected=None):
    if expected is None:
        expected = _intersection_helpers.BoxIntersectionType.intersection

    result = _intersection_helpers.bbox_line_intersect(
        UNIT_SQUARE, segment[[0], :], segment[[1], :])
    assert result is expected


def test_outside():
    segment_bottom_left = np.array([
        [0.25, -0.75],
        [-0.75, 0.25],
    ])
    segment_bottom_right = np.array([
        [0.75, -0.75],
        [1.75, 0.25],
    ])
    segment_top_right = np.array([
        [0.75, 1.75],
        [1.75, 0.75],
    ])
    segment_top_left = np.array([
        [-0.75, 0.75],
        [0.25, 1.75],
    ])

    segments = (
        segment_bottom_left,
        segment_bottom_right,
        segment_top_right,
        segment_top_left,
    )
    expected = _intersection_helpers.BoxIntersectionType.disjoint
    for segment in segments:
        run_it(segment, expected)


def test_unsupported():
    segment = np.array([
        [-0.5, 0.5],
        [1.5, 0.5],
    ])
    with pytest.raises(NotImplementedError):
        run_it(segment)


def test_start_in_box():
    segments = (
        np.array([
            [0.5, 0.0],
            [0.5, -0.5],
        ]),
        np.array([
            [0.5, 0.0],
            [0.5, 0.5],
        ]),
        np.array([
            [1.0, 0.5],
            [1.5, 0.5],
        ]),
        np.array([
            [1.0, 0.5],
            [0.5, 0.5],
        ]),
        np.array([
            [0.5, 1.0],
            [0.5, 1.5],
        ]),
        np.array([
            [0.5, 1.0],
            [0.5, 0.5],
        ]),
        np.array([
            [0.0, 0.5],
            [-0.5, 0.5],
        ]),
        np.array([
            [0.0, 0.5],
            [0.5, 0.5],
        ]),
        np.array([
            [0.5, 0.5],
            [1.5, 0.5],
        ]),
        np.array([
            [0.5, 0.5],
            [0.5, -0.5],
        ]),
        np.array([
            [0.5, 0.5],
            [-0.5, 1.5],
        ]),
    )

    for segment in segments:
        run_it(segment)


def test_end_in_box():
    segments = (
        np.array([
            [0.5, -1.0],
            [0.5, 0.0],
        ]),
        np.array([
            [0.25, -1.0],
            [0.5, 0.0],
        ]),
        np.array([
            [1.5, 0.5],
            [1.0, 0.5],
        ]),
        np.array([
            [1.0, -0.5],
            [1.0, 0.5],
        ]),
        np.array([
            [0.5, 1.5],
            [0.5, 1.0],
        ]),
        np.array([
            [0.5, -0.5],
            [0.5, 1.0],
        ]),
        np.array([
            [-0.5, 0.5],
            [0.0, 0.5],
        ]),
        np.array([
            [0.5, 1.5],
            [0.5, 0.5],
        ]),
        np.array([
            [-0.5, 0.5],
            [0.5, 0.5],
        ]),
        np.array([
            [1.5, 1.5],
            [0.5, 0.5],
        ]),
    )

    for segment in segments:
        run_it(segment)


def test_goes_through_box():
    segments = (
        np.array([
            [0.5, -0.25],
            [1.25, 0.5],
        ]),
        np.array([
            [-0.25, 0.5],
            [0.5, -0.25],
        ]),
        np.array([
            [0.0, -0.25],
            [0.0, 1.25],
        ]),
        np.array([
            [1.25, 0.5],
            [0.5, 1.25],
        ]),
        np.array([
            [0.5, 1.25],
            [-0.25, 0.5],
        ]),
    )

    for segment in segments:
        run_it(segment)
