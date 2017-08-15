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

import pytest

import runtime_utils


FAILURES = (11, 20, 24, 42)


def get_bounds():
    return 35.0 / 1024.0, 42.0 / 1024.0


def intersect_all(intersections):
    for intersection in intersections:
        try:
            intersection.curve1.intersect(intersection.curve2)
        except NotImplementedError:
            assert intersection.id_ in FAILURES


@pytest.mark.benchmark(
    group='curve-intersection',
    disable_gc=True,
    warmup=False
)
def test_intersections(benchmark):
    _, intersections = runtime_utils.curve_intersections_info()
    result = benchmark(intersect_all, intersections)
    assert result is None

    stats = benchmark.stats.stats
    min_time, max_time = get_bounds()
    assert min_time <= stats.median <= max_time
