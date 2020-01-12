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

import bezier
from tests.functional import utils


CONFIG = utils.Config()
# F1 = sympy.Matrix([[s, t]])
SURFACE1 = bezier.Surface.from_nodes(
    np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]), copy=False
)
# F2 = sympy.Matrix([[
#     (-t**2 + 2 * s + t) / 2, (s**2 + 2 * s * t - s + 2 * t) / 2]])
SURFACE2 = bezier.Surface.from_nodes(
    np.asfortranarray(
        [[0.0, 0.5, 1.0, 0.25, 0.75, 0.0], [0.0, -0.25, 0.0, 0.5, 0.75, 1.0]]
    ),
    copy=False,
)
# F3 = sympy.Matrix([[
#     -(2 * s * t - 4 * s - t) / 4, (s**2 - s * t + 4 * t) / 4]])
SURFACE3 = bezier.Surface.from_nodes(
    np.asfortranarray(
        [
            [0.0, 0.5, 1.0, 0.125, 0.375, 0.25],
            [0.0, 0.0, 0.25, 0.5, 0.375, 1.0],
        ]
    ),
    copy=False,
)
# F4 = sympy.Matrix([[2 * (s + 2 * t) * (1 - t), 2 * t * (s + 1)]])
SURFACE4 = bezier.Surface.from_nodes(
    np.asfortranarray(
        [[0.0, 1.0, 2.0, 2.0, 2.0, 0.0], [0.0, 0.0, 0.0, 1.0, 2.0, 2.0]]
    ),
    copy=False,
)
SURFACES = {
    1: SURFACE1,
    2: SURFACE2,
    3: SURFACE3,
    4: SURFACE4,
}
POINTS = np.asfortranarray(
    [[0.0, 0.25, 0.59375, 0.265625, 1.25], [0.0, 0.25, 0.25, 0.73046875, 1.25]]
)
S_VAL1, _ = utils.real_roots([4, -32, -56, -24, 5])
_, T_VAL1 = utils.real_roots([4, 8, -16, 44, -11])
(S_VAL2,) = utils.real_roots([2, -5, 15, -3])
(T_VAL2,) = utils.real_roots([14, -61, 74, -15])
(S_VAL3,) = utils.real_roots([64, 101, 34, -5])
(T_VAL3,) = utils.real_roots([128, -192, 91, -8])
CASES = (
    (1, 0, 0.0, 0.0),
    (1, 1, 0.25, 0.25),
    (2, 1, S_VAL1, T_VAL1),
    (2, 2, 0.5, 0.25),
    (2, 4, None, None),
    (3, 1, S_VAL2, T_VAL2),
    (3, 3, 0.125, 0.75),
    (4, 2, S_VAL3, T_VAL3),
    (4, 4, 0.25, 0.5),
)


def id_func(case):
    surface_index, point_index, _, _ = case
    return f"surface:{surface_index}-point:{point_index}"


@pytest.mark.parametrize("case", CASES, ids=id_func)
def test_locate(case):
    surface_index, point_index, expected_s, expected_t = case
    surface = SURFACES[surface_index]

    point = POINTS[:, [point_index]]
    if expected_s is None and expected_t is None:
        assert surface.locate(point) is None
    else:
        s, t = surface.locate(point)
        CONFIG.assert_close(s, expected_s)
        CONFIG.assert_close(t, expected_t)
