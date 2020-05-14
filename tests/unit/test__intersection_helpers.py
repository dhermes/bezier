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

from tests.unit import utils
from tests.unit.hazmat import test_intersection_helpers


@utils.needs_speedup
class Test_speedup_newton_refine(test_intersection_helpers.Test_newton_refine):
    @staticmethod
    def _call_function_under_test(s, nodes1, t, nodes2):
        from bezier import _speedup

        return _speedup.newton_refine_curve_intersect(s, nodes1, t, nodes2)
