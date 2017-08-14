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

import gc

import memory_profiler

import runtime_utils


FAILURES = (11, 20, 24, 42)


def intersect(curve1, curve2):
    return curve1.intersect(curve2)


def intersect_all():
    _, intersections = runtime_utils.curve_intersections_info()
    for intersection in intersections:
        try:
            intersection.curve1.intersect(intersection.curve2)
        except NotImplementedError:
            assert intersection.id_ in FAILURES


def test_main():
    # This should be the **only** test function here.
    gc.disable()
    intersect_all()
    kb_used = memory_profiler.memory_usage(max_usage=True)
    assert 28 <= kb_used <= 32
    gc.enable()


if __name__ == '__main__':
    test_main()
