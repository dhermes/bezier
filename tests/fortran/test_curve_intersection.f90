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

module test_curve_intersection

  use iso_c_binding, only: c_bool, c_double
  use curve, only: subdivide_nodes
  use curve_intersection, only: linearization_error
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private test_linearization_error
  public curve_intersection_all_tests

contains

  subroutine curve_intersection_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_linearization_error(success)

  end subroutine curve_intersection_all_tests

  subroutine test_linearization_error(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes1(2, 2), nodes2(3, 2), nodes3(5, 2)
    real(c_double) :: nodes4(3, 3), nodes5(4, 2), nodes6(6, 2)
    real(c_double) :: left_nodes2(3, 2), right_nodes2(3, 2)
    real(c_double) :: error1, error2, expected
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "linearization_error"

    ! CASE 1: Linear curve (i.e. no error).
    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 2.0_dp]
    call linearization_error( &
         1, 2, nodes1, error1)
    if (error1 == 0.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Degree-elevated line (i.e. no error) as quadratic.
    nodes2(1, :) = 0
    nodes2(2, :) = [0.5_dp, 1.0_dp]
    nodes2(3, :) = [1.0_dp, 2.0_dp]
    call linearization_error( &
         2, 2, nodes2, error1)
    if (error1 == 0.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: Degree-elevated line (i.e. no error) as quartic.
    nodes3(1, :) = 0
    nodes3(2, :) = [0.25_dp, 0.5_dp]
    nodes3(3, :) = [0.5_dp, 1.0_dp]
    nodes3(4, :) = [0.75_dp, 1.5_dp]
    nodes3(5, :) = [1.0_dp, 2.0_dp]
    call linearization_error( &
         4, 2, nodes3, error1)
    if (error1 == 0.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 4: Line with bad parameterization.
    ! NOTE: This is the line 3 y = 4 x, but with the parameterization
    !       x(s) = 3 s (4 - 3 s).
    nodes2(1, :) = 0
    nodes2(2, :) = [6.0_dp, 8.0_dp]
    nodes2(3, :) = [3.0_dp, 4.0_dp]
    call linearization_error( &
         2, 2, nodes2, error1)
    ! D^2 v = [-9, -12]
    expected = 0.125_dp * 2 * 1 * 15.0_dp
    if (error1 == expected) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 5: Quadratic curve.
    nodes2(1, :) = 0
    nodes2(2, :) = [1.0_dp, 1.0_dp]
    nodes2(3, :) = [5.0_dp, 6.0_dp]
    call linearization_error( &
         2, 2, nodes2, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [1, 1], [4, 5]
    !           d^2 Nodes = [3, 4]
    !       so that sqrt(3^2 + 4^2) = 5.0
    expected = 0.125_dp * 2 * 1 * 5.0_dp
    if (error1 == expected) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 6: Subdivided curves (left and right) from CASE 5.
    call subdivide_nodes( &
         3, 2, nodes2, left_nodes2, right_nodes2)
    call linearization_error( &
         2, 2, left_nodes2, error1)
    call linearization_error( &
         2, 2, right_nodes2, error2)
    ! For a degree two curve, the 2nd derivative is constant
    ! so by subdividing, our error should drop by a factor
    ! of (1/2)^2 = 4.
    expected = 0.25_dp * expected
    if (error1 == expected .AND. error2 == expected) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 7: Quadratic curve in 3D.
    nodes4(1, :) = [1.5_dp, 0.0_dp, 6.25_dp]
    nodes4(2, :) = [3.5_dp, -5.0_dp, 10.25_dp]
    nodes4(3, :) = [8.5_dp, 2.0_dp, 10.25_dp]
    call linearization_error( &
         2, 3, nodes4, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [2, -5, 4], [5, 7, 0]
    !           d^2 Nodes = [3, 12, -4]
    !       so that sqrt(3^2 + 12^2 + 4^2) = 13.0
    expected = 0.125_dp * 2 * 1 * 13.0_dp
    if (error1 == expected) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 8: Quadratic with bad parameterization.
    ! NOTE: This is the quadratic y = 1 + x^2 / 4, but with the
    !       parameterization x(s) = (3 s - 1)^2.
    nodes3(1, :) = [1.0_dp, 1.25_dp]
    nodes3(2, :) = [-0.5_dp, 0.5_dp]
    nodes3(3, :) = [-0.5_dp, 2.0_dp]
    nodes3(4, :) = [1.0_dp, -1.0_dp]
    nodes3(5, :) = [4.0_dp, 5.0_dp]
    call linearization_error( &
         4, 2, nodes3, error1)
    ! D^2 v = [1.5, 2.25], [1.5, -4.5], [1.5, 9]
    expected = 0.125_dp * 4 * 3 * sqrt(1.5_dp**2 + 9.0_dp**2)
    if (abs(error1 - expected) <= spacing(expected)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 9: Cubic curve.
    nodes5(1, :) = 0
    nodes5(2, :) = 1
    nodes5(3, :) = [5.0_dp, 6.0_dp]
    nodes5(4, :) = [6.0_dp, 7.0_dp]
    call linearization_error( &
         3, 2, nodes5, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [1, 1], [4, 5], [1, 1]
    !           d^2 Nodes = [3, 4], [-3, -4]
    !       so that sqrt(3^2 + 4^2) = 5.0
    expected = 0.125_dp * 3 * 2 * 5.0_dp
    if (error1 == expected) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 10: Quartic curve.
    nodes3(1, :) = 0
    nodes3(2, :) = 1
    nodes3(3, :) = [5.0_dp, 6.0_dp]
    nodes3(4, :) = [6.0_dp, 7.0_dp]
    nodes3(5, :) = [4.0_dp, 7.0_dp]
    call linearization_error( &
         4, 2, nodes3, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [1, 1], [4, 5], [1, 1], [-2, 0]
    !           d^2 Nodes = [3, 4], [-3, -4], [-3, -1]
    !       so that sqrt(3^2 + 4^2) = 5.0
    expected = 0.125_dp * 4 * 3 * 5.0_dp
    if (error1 == expected) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 11: Quintic curve (i.e. degree 5).
    nodes6(1, :) = [0.0_dp, 0.0_dp]
    nodes6(2, :) = [1.0_dp, 1.0_dp]
    nodes6(3, :) = [7.0_dp, 3.0_dp]
    nodes6(4, :) = [11.0_dp, 8.0_dp]
    nodes6(5, :) = [15.0_dp, 1.0_dp]
    nodes6(6, :) = [16.0_dp, -3.0_dp]
    call linearization_error( &
         5, 2, nodes6, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [1, 1], [6, 2], [4, 5], [4, -7], [1, -4]
    !           d^2 Nodes = [5, 1], [-2, 3], [0, -12], [-3, 3]
    !       so that sqrt(5^2 + 12^2) = 13.0
    expected = 0.125_dp * 5 * 4 * 13.0_dp
    if (error1 == expected) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_linearization_error

end module test_curve_intersection
