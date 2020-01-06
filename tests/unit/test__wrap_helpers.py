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

from tests.unit import test__py_helpers
from tests.unit import utils


@utils.needs_speedup
class Test_speedup_vector_close(test__py_helpers.Test__vector_close):
    @staticmethod
    def _call_function_under_test(vec1, vec2, **kwargs):
        from bezier import _speedup

        return _speedup.vector_close(vec1, vec2, **kwargs)


@utils.needs_speedup
class Test_speedup_in_interval(test__py_helpers.Test__in_interval):
    @staticmethod
    def _call_function_under_test(value, start, end):
        from bezier import _speedup

        return _speedup.in_interval(value, start, end)


@utils.needs_speedup
class Test_speedup_bbox(test__py_helpers.Test__bbox):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _speedup

        return _speedup.bbox(nodes)


@utils.needs_speedup
class Test_speedup_contains_nd(test__py_helpers.Test__contains_nd):
    @staticmethod
    def _call_function_under_test(nodes, point):
        from bezier import _speedup

        return _speedup.contains_nd(nodes, point)


@utils.needs_speedup
class Test_speedup_cross_product(test__py_helpers.Test__cross_product):
    @staticmethod
    def _call_function_under_test(vec0, vec1):
        from bezier import _speedup

        return _speedup.cross_product(vec0, vec1)


@utils.needs_speedup
class Test_speedup_wiggle_interval(test__py_helpers.Test__wiggle_interval):
    def _call_function_under_test(self, value, **kwargs):
        from bezier import _speedup

        self.assertEqual(kwargs, {})
        return _speedup.wiggle_interval(value, **kwargs)

    def test_custom_wiggle(self):
        # Fortran implementation doesn't support optional wiggle. This
        # isn't because Fortran **can't** (just use "optional"), it's just
        # to allow the compiler to pre-compute 1 + wiggle / 1 - wiggle
        # rather than having to deal with it at run-time.
        pass


@utils.needs_speedup
class Test_speedup_simple_convex_hull(
    test__py_helpers.Test__simple_convex_hull
):
    @staticmethod
    def _call_function_under_test(points):
        from bezier import _speedup

        return _speedup.simple_convex_hull(points)


@utils.needs_speedup
class Test_speedup_polygon_collide(test__py_helpers.Test__polygon_collide):
    @staticmethod
    def _call_function_under_test(polygon1, polygon2):
        from bezier import _speedup

        return _speedup.polygon_collide(polygon1, polygon2)
