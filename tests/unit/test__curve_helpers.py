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

from tests.unit import test__py_curve_helpers
from tests.unit import utils


@utils.needs_speedup
class Test_speedup_subdivide_nodes(
    test__py_curve_helpers.Test_subdivide_nodes
):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _speedup

        return _speedup.subdivide_nodes_curve(nodes)


@utils.needs_speedup
class Test_speedup_evaluate_multi_barycentric(
    test__py_curve_helpers.Test_evaluate_multi_barycentric
):
    @staticmethod
    def _call_function_under_test(nodes, lambda1, lambda2):
        from bezier import _speedup

        return _speedup.evaluate_multi_barycentric(nodes, lambda1, lambda2)


@utils.needs_speedup
class Test_speedup_evaluate_multi(test__py_curve_helpers.Test_evaluate_multi):
    @staticmethod
    def _call_function_under_test(nodes, s_vals):
        from bezier import _speedup

        return _speedup.evaluate_multi(nodes, s_vals)


@utils.needs_speedup
class Test_speedup_compute_length(test__py_curve_helpers.Test_compute_length):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _speedup

        return _speedup.compute_length(nodes)

    def _scipy_skip(self):
        # Fortran implementation directly includes QUADPACK, so the presence
        # or absence of SciPy is irrelevant.
        pass

    def test_without_scipy(self):
        # Fortran implementation directly includes QUADPACK, so the presence
        # or absence of SciPy is irrelevant.
        pass


@utils.needs_speedup
class Test_speedup_elevate_nodes(test__py_curve_helpers.Test_elevate_nodes):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _speedup

        return _speedup.elevate_nodes(nodes)


@utils.needs_speedup
class Test_speedup_specialize_curve(
    test__py_curve_helpers.Test_specialize_curve
):
    @staticmethod
    def _call_function_under_test(nodes, start, end):
        from bezier import _speedup

        return _speedup.specialize_curve(nodes, start, end)


@utils.needs_speedup
class Test_speedup_evaluate_hodograph(
    test__py_curve_helpers.Test_evaluate_hodograph
):
    @staticmethod
    def _call_function_under_test(s, nodes):
        from bezier import _speedup

        return _speedup.evaluate_hodograph(s, nodes)


@utils.needs_speedup
class Test_speedup_get_curvature(test__py_curve_helpers.Test_get_curvature):
    @staticmethod
    def _call_function_under_test(nodes, tangent_vec, s):
        from bezier import _speedup

        return _speedup.get_curvature(nodes, tangent_vec, s)


@utils.needs_speedup
class Test_speedup_newton_refine(test__py_curve_helpers.Test_newton_refine):
    @staticmethod
    def _call_function_under_test(nodes, point, s):
        from bezier import _speedup

        return _speedup.newton_refine_curve(nodes, point, s)


@utils.needs_speedup
class Test_speedup_locate_point(test__py_curve_helpers.Test_locate_point):
    @staticmethod
    def _call_function_under_test(nodes, point):
        from bezier import _speedup

        return _speedup.locate_point_curve(nodes, point)


@utils.needs_speedup
class Test_speedup_reduce_pseudo_inverse(
    test__py_curve_helpers.Test_reduce_pseudo_inverse
):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _speedup

        return _speedup.reduce_pseudo_inverse(nodes)


@utils.needs_speedup
class Test_speedup_full_reduce(test__py_curve_helpers.Test_full_reduce):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _speedup

        return _speedup.full_reduce(nodes)
