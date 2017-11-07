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

module test_surface_intersection

  use, intrinsic :: iso_c_binding, only: c_bool, c_double, c_int
  use curve, only: LOCATE_MISS
  use surface_intersection, only: newton_refine, locate_point
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private test_newton_refine, test_locate_point
  public surface_intersection_all_tests

contains

  subroutine surface_intersection_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_newton_refine(success)
    call test_locate_point(success)

  end subroutine surface_intersection_all_tests

  subroutine test_newton_refine(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(:), allocatable :: name
    real(c_double) :: nodes(6, 2)
    real(c_double) :: updated_s, updated_t

    case_id = 1
    name = "newton_refine (Surface)"

    ! CASE 1: Quadratic surface.
    ! NOTE: This surface is given by
    !     [(4 s - t^2) / 4, (4 s^2 + 4 s t - t^2 - 4 s + 8 t) / 8]
    nodes(1, :) = [0.0_dp, 0.0_dp]
    nodes(2, :) = [0.5_dp, -0.25_dp]
    nodes(3, :) = [1.0_dp, 0.0_dp]
    nodes(4, :) = [0.0_dp, 0.5_dp]
    nodes(5, :) = [0.5_dp, 0.5_dp]
    nodes(6, :) = [-0.25_dp, 0.875_dp]
    ! At our point (s, t) = (1/4, 1/2), the Jacobian is
    !     [1, -1/4]
    !     [0,  1  ]
    ! hence there will be no round-off when applying the inverse.
    call newton_refine( &
         6, nodes, 2, 0.484375_dp, 0.1796875_dp, &
         0.25_dp, 0.5_dp, updated_s, updated_t)
    case_success = ( &
         updated_s == 247.0_dp / 512 .AND. &
         updated_t == 31.0_dp / 128)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Solution is **exact**.
    ! NOTE: This surface is given by [s, t]
    nodes(1, :) = [0.0_dp, 0.0_dp]
    nodes(2, :) = [0.5_dp, 0.0_dp]
    nodes(3, :) = [1.0_dp, 0.0_dp]
    nodes(4, :) = [0.0_dp, 0.5_dp]
    nodes(5, :) = [0.5_dp, 0.5_dp]
    nodes(6, :) = [0.0_dp, 1.0_dp]
    ! Since x(s) = s and y(t) = t, we simply use the same x/y and s/t.
    call newton_refine( &
         6, nodes, 2, 0.375_dp, 0.75_dp, &
         0.375_dp, 0.75_dp, updated_s, updated_t)
    case_success = ( &
         updated_s == 0.375_dp .AND. &
         updated_t == 0.75_dp)
    call print_status(name, case_id, case_success, success)

  end subroutine test_newton_refine

  subroutine test_locate_point(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    real(c_double), allocatable :: nodes(:, :)
    integer(c_int) :: degree
    real(c_double) :: x_val, y_val, s_val, t_val
    real(c_double) :: expected_s, expected_t
    character(:), allocatable :: name

    case_id = 1
    name = "newton_refine (Surface)"

    ! CASE 1: B(s, t) = (s, t), match.
    allocate(nodes(3, 2))
    nodes(1, :) = 0
    nodes(2, :) = [1.0_dp, 0.0_dp]
    nodes(3, :) = [0.0_dp, 1.0_dp]
    x_val = 0.25_dp
    y_val = 0.625_dp
    call locate_point( &
         3, nodes, 1, x_val, y_val, s_val, t_val)
    case_success = (x_val == s_val .AND. y_val == t_val)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: B(s, t) = (s, t), no match.
    x_val = -0.125_dp
    y_val = 0.25_dp
    call locate_point( &
         3, nodes, 1, x_val, y_val, s_val, t_val)
    case_success = (s_val == LOCATE_MISS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: B(s, t) = (-2 (s + 2 t) (t - 1), 2 (s + 1) t), one extra
    !         Newton refinement required (verified by hand).
    deallocate(nodes)
    allocate(nodes(6, 2))
    nodes(1, :) = 0
    nodes(2, :) = [1.0_dp, 0.0_dp]
    nodes(3, :) = [2.0_dp, 0.0_dp]
    nodes(4, :) = [2.0_dp, 1.0_dp]
    nodes(5, :) = [2.0_dp, 2.0_dp]
    nodes(6, :) = [0.0_dp, 2.0_dp]
    x_val = 0.59375_dp
    y_val = 0.25_dp
    call locate_point( &
         6, nodes, 2, x_val, y_val, s_val, t_val)
    ! NOTE: We can use the resultant to find that the **true** answers
    !       are roots of the following polynomials.
    !           64 s^3 + 101 s^2 + 34 s - 5 = 0
    !           128 t^3 - 192 t^2 + 91 t - 8 = 0
    !       Using extended precision, we can find these values to more
    !       digits than what is supported by IEEE-754.
    expected_s = 0.109190958136897160638_dp
    expected_t = 0.11269475204698919699_dp
    case_success = ( &
         abs(s_val - expected_s) <= 5 * spacing(expected_s) .AND. &
         abs(t_val - expected_t) <= spacing(expected_t))
    call print_status(name, case_id, case_success, success)

  end subroutine test_locate_point

end module test_surface_intersection
