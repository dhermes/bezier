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

from __future__ import absolute_import

import os

from tests.functional import utils


SCRIPTS_DIR = os.path.abspath(os.path.dirname(__file__))
GENERATED_FILENAME = os.path.abspath(os.path.join(
    SCRIPTS_DIR, '..', 'tests', 'fortran', 'functional_curve.f90'))
COL_TEMPLATE = '    nodes{:d}(:, {:d}) = [{}_dp, {}_dp]'
PARAM_TEMPLATE_COMPACT = '    expected_params({:d}, :) = [{}]'
PARAM_TEMPLATE = '    expected_params({:d}, :) = [ &\n         {}]'
PARAM_JOIN_ON = ', &\n         '
CASE_TEMPLATE = """\
  subroutine case{case_id:d}()

    real(c_double) :: nodes{id1}(2, {num_nodes1:d})
    real(c_double) :: nodes{id2}(2, {num_nodes2:d})
    real(c_double) :: expected_params(2, {expected_n:d})

{nodes1_vals}

{nodes2_vals}

{expected_params}    call intersect_and_check( &
         {case_id:d}, &
         {num_nodes1:d}, nodes{id1}, &
         {num_nodes2:d}, nodes{id2}, &
         {expected_n:d}, expected_params)

  end subroutine case{case_id:d}
"""
FILE_TEMPLATE = """\
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     https://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module functional_curve

  use, intrinsic :: iso_c_binding, only: c_double
  use types, only: dp
  use functional_test_helpers, only: intersect_and_check
  implicit none
  private &
{cases}
  public all_cases

contains

{case_impls}
  subroutine all_cases()

{case_calls}

  end subroutine all_cases

end module functional_curve
"""


def params_string(index, curve_params):
    inside = ['{}_dp'.format(param)
              for param in curve_params]
    result_compact = PARAM_TEMPLATE_COMPACT.format(index, ', '.join(inside))
    if len(result_compact) <= 80:
        return result_compact

    return PARAM_TEMPLATE.format(index, PARAM_JOIN_ON.join(inside))


def make_test_case(intersection):
    id1 = int(intersection.curve1_info.id_)
    nodes1 = intersection.nodes1
    _, num_nodes1 = nodes1.shape

    id2 = int(intersection.curve2_info.id_)
    nodes2 = intersection.nodes2
    _, num_nodes2 = nodes2.shape

    num_params = intersection.num_params

    parts = []
    for col in range(num_nodes1):
        x_val, y_val = nodes1[:, col]
        parts.append(COL_TEMPLATE.format(id1, col + 1, x_val, y_val))
    nodes1_vals = '\n'.join(parts)

    parts = []
    for col in range(num_nodes2):
        x_val, y_val = nodes2[:, col]
        parts.append(COL_TEMPLATE.format(id2, col + 1, x_val, y_val))
    nodes2_vals = '\n'.join(parts)

    if num_params == 0:
        expected_params = ''
    else:
        part1 = params_string(1, intersection.curve1_params)
        part2 = params_string(2, intersection.curve2_params)
        expected_params = part1 + '\n' + part2 + '\n\n'

    return CASE_TEMPLATE.format(
        case_id=intersection.id_,
        id1=id1,
        num_nodes1=num_nodes1,
        id2=id2,
        num_nodes2=num_nodes2,
        expected_n=num_params,
        nodes1_vals=nodes1_vals,
        nodes2_vals=nodes2_vals,
        expected_params=expected_params,
    )


def write_content(intersections, file_obj):
    parts = []
    for intersection in intersections[:-1]:
        parts.append('       case{:d}, &'.format(intersection.id_))
    parts.append('       case{:d}'.format(intersections[-1].id_))
    cases = '\n'.join(parts)

    parts = []
    for intersection in intersections:
        parts.append(make_test_case(intersection))
    case_impls = '\n'.join(parts)

    parts = []
    for intersection in intersections:
        parts.append('    call case{:d}()'.format(intersection.id_))
    case_calls = '\n'.join(parts)

    content = FILE_TEMPLATE.format(
        cases=cases,
        case_impls=case_impls,
        case_calls=case_calls,
    )
    file_obj.write(content)


def main():
    _, intersections = utils.curve_intersections_info()

    with open(GENERATED_FILENAME, 'w') as file_obj:
        write_content(intersections, file_obj)


if __name__ == '__main__':
    main()
