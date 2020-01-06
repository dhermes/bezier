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

from tests.unit import test__py_surface_helpers
from tests.unit import utils


@utils.needs_speedup
class Test_speedup_de_casteljau_one_round(
    test__py_surface_helpers.Test_de_casteljau_one_round
):
    @staticmethod
    def _call_function_under_test(nodes, degree, lambda1, lambda2, lambda3):
        from bezier import _speedup

        return _speedup.de_casteljau_one_round(
            nodes, degree, lambda1, lambda2, lambda3
        )


@utils.needs_speedup
class Test_speedup_specialize_surface(
    test__py_surface_helpers.Test_specialize_surface
):
    @staticmethod
    def _call_function_under_test(
        nodes, degree, weights_a, weights_b, weights_c
    ):
        from bezier import _speedup

        return _speedup.specialize_surface(
            nodes, degree, weights_a, weights_b, weights_c
        )


@utils.needs_speedup
class Test_speedup_subdivide_nodes(
    test__py_surface_helpers.Test_subdivide_nodes
):
    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _speedup

        return _speedup.subdivide_nodes_surface(nodes, degree)


@utils.needs_speedup
class Test_speedup_jacobian_both(test__py_surface_helpers.Test_jacobian_both):
    @staticmethod
    def _call_function_under_test(nodes, degree, dimension):
        from bezier import _speedup

        return _speedup.jacobian_both(nodes, degree, dimension)


@utils.needs_speedup
class Test_speedup_jacobian_det(test__py_surface_helpers.Test_jacobian_det):
    @staticmethod
    def _call_function_under_test(nodes, degree, st_vals):
        from bezier import _speedup

        return _speedup.jacobian_det(nodes, degree, st_vals)


@utils.needs_speedup
class Test_speedup_evaluate_barycentric(
    test__py_surface_helpers.Test_evaluate_barycentric
):
    @staticmethod
    def _call_function_under_test(nodes, degree, lambda1, lambda2, lambda3):
        from bezier import _speedup

        return _speedup.evaluate_barycentric(
            nodes, degree, lambda1, lambda2, lambda3
        )


@utils.needs_speedup
class Test_speedup_evaluate_barycentric_multi(
    test__py_surface_helpers.Test_evaluate_barycentric_multi
):
    @staticmethod
    def _call_function_under_test(nodes, degree, param_vals, dimension):
        from bezier import _speedup

        return _speedup.evaluate_barycentric_multi(
            nodes, degree, param_vals, dimension
        )


@utils.needs_speedup
class Test_speedup_evaluate_cartesian_multi(
    test__py_surface_helpers.Test_evaluate_cartesian_multi
):
    @staticmethod
    def _call_function_under_test(nodes, degree, param_vals, dimension):
        from bezier import _speedup

        return _speedup.evaluate_cartesian_multi(
            nodes, degree, param_vals, dimension
        )


@utils.needs_speedup
class Test_speedup_compute_edge_nodes(
    test__py_surface_helpers.Test_compute_edge_nodes
):
    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _speedup

        return _speedup.compute_edge_nodes(nodes, degree)


@utils.needs_speedup
class Test_speedup_compute_area(test__py_surface_helpers.Test_compute_area):
    @staticmethod
    def _call_function_under_test(edges):
        from bezier import _speedup

        return _speedup.compute_area(edges)
