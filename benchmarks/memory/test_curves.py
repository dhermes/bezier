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

from __future__ import print_function

import gc
import os
import sys

import memory_profiler

import runtime_utils


FAILURES = (11, 20, 24, 42)
ERR_TEMPLATE = 'Memory usage {:g} outside of expected range {}-{}KB.'
SUCCESS_TEMPLATE = 'Memory usage: {:g}KB.'


def get_bounds():
    # NOTE: These bounds assume **just** the interpeter is running this code.
    #       When using a test runner like `py.test`, usage goes up by 4-8 KB.
    if os.getenv('CIRCLECI') == 'true':
        return 22, 30
    else:
        return 28, 32


def intersect_all():
    _, intersections = runtime_utils.curve_intersections_info()
    for intersection in intersections:
        try:
            intersection.curve1.intersect(intersection.curve2)
        except NotImplementedError:
            assert intersection.id_ in FAILURES


def test_main():
    min_kb, max_kb = get_bounds()

    # This should be the **only** test function here.
    gc.disable()
    intersect_all()
    kb_used_after = memory_profiler.memory_usage(max_usage=True)
    if min_kb <= kb_used_after <= max_kb:
        status = 0
        msg = SUCCESS_TEMPLATE.format(kb_used_after)
        print(msg)
    else:
        status = 1
        msg = ERR_TEMPLATE.format(kb_used_after, min_kb, max_kb)
        print(msg, file=sys.stderr)

    gc.enable()
    sys.exit(status)


if __name__ == '__main__':
    test_main()
