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

module test_curve

  use iso_c_binding, only: c_bool, c_double
  use curve, only: evaluate_curve_barycentric, evaluate_multi
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private test_evaluate_curve_barycentric, test_evaluate_multi
  public curve_all_tests

contains

  subroutine curve_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_evaluate_curve_barycentric(success)
    call test_evaluate_multi(success)

  end subroutine curve_all_tests

  subroutine test_evaluate_curve_barycentric(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes(4, 3)
    real(c_double) :: lambda1(3), lambda2(3)
    real(c_double) :: evaluated(3, 3), expected(3, 3)
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "evaluate_curve_barycentric"

    ! CASE 1: Slightly complex example.
    nodes(1, :) = [0.0_dp, 0.0_dp, 0.0_dp]
    nodes(2, :) = [0.5_dp, 3.0_dp, 1.0_dp]
    nodes(3, :) = [1.5_dp, 4.0_dp, 1.0_dp]
    nodes(4, :) = [2.0_dp, 8.0_dp, 1.0_dp]
    lambda1 = [0.25_dp, 0.5_dp, 0.75_dp]
    lambda2 = [0.25_dp, 0.125_dp, -0.75_dp]
    expected(1, :) = [0.125_dp, 0.453125_dp, 0.109375_dp]
    expected(2, :) = [0.0859375_dp, 0.390625_dp, 0.119140625_dp]
    expected(3, :) = [0.421875_dp, -2.109375_dp, -0.421875_dp]
    call evaluate_curve_barycentric( &
         3, 3, nodes, 3, lambda1, lambda2, evaluated)

    if (all(evaluated == expected)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_evaluate_curve_barycentric

  subroutine test_evaluate_multi(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes1(2, 3), nodes2(3, 2)
    real(c_double) :: s_vals1(129), s_vals2(65)
    real(c_double) :: evaluated1(129, 3), expected1(129, 3)
    real(c_double) :: evaluated2(65, 2), expected2(65, 2)
    integer :: i
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "evaluate_multi"

    ! CASE 1: Linear curve.
    do i = 1, 129
       s_vals1(i) = (i - 1) / 128.0_dp
    end do
    nodes1(1, :) = [1.0_dp, 1.0_dp, -7.0_dp]
    nodes1(2, :) = [2.0_dp, -1.0_dp, -4.0_dp]
    ! B(s) = [s + 1, 1 - 2 s, 3 s - 7]
    expected1(:, 1) = 1.0_dp + s_vals1
    expected1(:, 2) = 1.0_dp - 2.0_dp * s_vals1
    expected1(:, 3) = -7.0_dp + 3.0_dp * s_vals1
    call evaluate_multi( &
         1, 3, nodes1, 129, s_vals1, evaluated1)

    if (all(evaluated1 == expected1)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Quadratic curve.
    do i = 1, 65
       s_vals2(i) = (i - 1) / 64.0_dp
    end do
    nodes2(1, :) = [0.0_dp, 0.0_dp]
    nodes2(2, :) = [2.0_dp, -1.0_dp]
    nodes2(3, :) = [3.0_dp, 2.0_dp]
    ! B(s) = [s(4 - s), 2s(2s - 1)]
    expected2(:, 1) = s_vals2 * (4.0_dp - s_vals2)
    expected2(:, 2) = 2.0_dp * s_vals2 * (2.0_dp * s_vals2 - 1.0_dp)
    call evaluate_multi( &
         2, 2, nodes2, 65, s_vals2, evaluated2)

    if (all(evaluated2 == expected2)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_evaluate_multi

end module test_curve
