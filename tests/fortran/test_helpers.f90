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

module test_helpers

  use, intrinsic :: iso_c_binding, only: c_double, c_bool
  use helpers, only: &
       WIGGLE, cross_product, bbox, wiggle_interval, contains_nd, &
       vector_close, in_interval, ulps_away
  use types, only: dp
  use unit_test_helpers, only: MACHINE_EPS, print_status
  implicit none
  private &
       test_cross_product, test_bbox, test_wiggle_interval, &
       test_contains_nd, test_vector_close, test_in_interval, test_ulps_away
  public helpers_all_tests

contains

  subroutine helpers_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_cross_product(success)
    call test_bbox(success)
    call test_wiggle_interval(success)
    call test_contains_nd(success)
    call test_vector_close(success)
    call test_in_interval(success)
    call test_ulps_away(success)

  end subroutine helpers_all_tests

  subroutine test_cross_product(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: vec0(1, 2)
    real(c_double) :: vec1(1, 2)
    real(c_double) :: result_
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "cross_product"

    ! CASE 1: Identical vector.
    vec0(1, :) = [1.0_dp, 7.0_dp] / 8
    vec1(1, :) = [-11.0_dp, 24.0_dp] / 32
    call cross_product(vec0, vec1, result_)
    if (result_ == 101.0_dp / 256) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_cross_product

  subroutine test_bbox(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes(6, 2)
    real(c_double) :: left, right, bottom, top
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "bbox"

    ! CASE 1: Simple case (just two nodes).
    nodes(1, :) = [0.0_dp, 5.0_dp]
    nodes(2, :) = [1.0_dp, 3.0_dp]
    call bbox(2, nodes(:2, :), left, right, bottom, top)
    if (left == 0.0_dp .AND. right == 1.0_dp &
         .AND. bottom == 3.0_dp .AND. top == 5.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: "Many" nodes.
    nodes(1, :) = [1.0_dp, 0.0_dp]
    nodes(2, :) = [2.0_dp, 1.0_dp]
    nodes(3, :) = [-1.0_dp, 2.0_dp]
    nodes(4, :) = [5.0_dp, -3.0_dp]
    nodes(5, :) = [4.0_dp, 4.0_dp]
    nodes(6, :) = [0.0_dp, 0.0_dp]
    call bbox(6, nodes, left, right, bottom, top)
    if (left == -1.0_dp .AND. right == 5.0_dp &
         .AND. bottom == -3.0_dp .AND. top == 4.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_bbox

  subroutine test_wiggle_interval(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: value1, value2, value3, value4, value5
    real(c_double) :: result1, result2, result3, result4, result5
    logical(c_bool) :: &
         wiggle_success1, wiggle_success2, wiggle_success3, &
         wiggle_success4, wiggle_success5
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "wiggle_interval"

    ! CASE 1: **Exactly** at endpoint.
    value1 = 0.0_dp
    call wiggle_interval(value1, result1, wiggle_success1)
    if (wiggle_success1 .AND. result1 == 0.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Near endpoint.
    value1 = 1.0_dp + 0.5_dp**20
    call wiggle_interval(value1, result1, wiggle_success1)
    if (.NOT. wiggle_success1) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: "Below" the interval.
    value1 = -0.25_dp
    call wiggle_interval(value1, result1, wiggle_success1)
    if (.NOT. wiggle_success1) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 4: "Above" the interval.
    value1 = 1.5_dp
    call wiggle_interval(value1, result1, wiggle_success1)
    if (.NOT. wiggle_success1) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 5: Inside the interval.
    value1 = 0.25_dp
    call wiggle_interval(value1, result1, wiggle_success1)
    if (wiggle_success1 .AND. result1 == value1) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 6: Wiggle "below" the interval.
    value1 = -0.5_dp**60
    call wiggle_interval(value1, result1, wiggle_success1)
    if (wiggle_success1 .AND. result1 == 0.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 7: Wiggle "above" the interval.
    value1 = 1.0_dp + MACHINE_EPS
    call wiggle_interval(value1, result1, wiggle_success1)
    if (wiggle_success1 .AND. result1 == 1.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 8: Test "lower outer" boundary, i.e. values in (-epsilon, 0).
    value1 = -WIGGLE + MACHINE_EPS * WIGGLE
    value2 = -WIGGLE + MACHINE_EPS * WIGGLE / 2
    value3 = -WIGGLE
    value4 = -WIGGLE - MACHINE_EPS * WIGGLE
    call wiggle_interval(value1, result1, wiggle_success1)
    call wiggle_interval(value2, result2, wiggle_success2)
    call wiggle_interval(value3, result3, wiggle_success3)
    call wiggle_interval(value4, result4, wiggle_success4)
    if (wiggle_success1 .AND. result1 == 0.0_dp .AND. &
         wiggle_success2 .AND. result2 == 0.0_dp .AND. &
         .NOT. wiggle_success3 .AND. .NOT. wiggle_success4) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 9: Test "upper outer" boundary, i.e. values in (1, 1 + epsilon).
    value1 = 1.0_dp + WIGGLE - 2 * MACHINE_EPS
    value2 = 1.0_dp + WIGGLE - MACHINE_EPS
    value3 = 1.0_dp + WIGGLE
    value4 = 1.0_dp + WIGGLE + MACHINE_EPS
    call wiggle_interval(value1, result1, wiggle_success1)
    call wiggle_interval(value2, result2, wiggle_success2)
    call wiggle_interval(value3, result3, wiggle_success3)
    call wiggle_interval(value4, result4, wiggle_success4)
    if (wiggle_success1 .AND. result1 == 1.0_dp .AND. &
         wiggle_success2 .AND. result2 == 1.0_dp .AND. &
         .NOT. wiggle_success3 .AND. .NOT. wiggle_success4) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 10: Test "lower inner" boundary, i.e. values in (0, epsilon).
    value1 = WIGGLE - WIGGLE * MACHINE_EPS
    value2 = WIGGLE - WIGGLE * MACHINE_EPS / 2
    value3 = WIGGLE
    value4 = WIGGLE + WIGGLE * MACHINE_EPS
    call wiggle_interval(value1, result1, wiggle_success1)
    call wiggle_interval(value2, result2, wiggle_success2)
    call wiggle_interval(value3, result3, wiggle_success3)
    call wiggle_interval(value4, result4, wiggle_success4)
    if (wiggle_success1 .AND. result1 == 0.0_dp .AND. &
         wiggle_success2 .AND. result2 == 0.0_dp .AND. &
         wiggle_success3 .AND. result3 == value3 .AND. &
         wiggle_success4 .AND. result4 == value4) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 11: Test "uper inner" boundary, i.e. values in (1 - epsilon, 1).
    value1 = 1.0_dp - WIGGLE - MACHINE_EPS
    value2 = 1.0_dp - WIGGLE - MACHINE_EPS / 2
    value3 = 1.0_dp - WIGGLE
    value4 = 1.0_dp - WIGGLE + MACHINE_EPS / 2
    value5 = 1.0_dp - WIGGLE + MACHINE_EPS
    call wiggle_interval(value1, result1, wiggle_success1)
    call wiggle_interval(value2, result2, wiggle_success2)
    call wiggle_interval(value3, result3, wiggle_success3)
    call wiggle_interval(value4, result4, wiggle_success4)
    call wiggle_interval(value5, result5, wiggle_success5)
    if (wiggle_success1 .AND. result1 == value1 .AND. &
         wiggle_success2 .AND. result2 == value2 .AND. &
         wiggle_success3 .AND. result3 == value3 .AND. &
         wiggle_success4 .AND. result4 == 1.0_dp .AND. &
         wiggle_success5 .AND. result5 == 1.0_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_wiggle_interval

  subroutine test_contains_nd(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes1(3, 2), nodes2(2, 3), nodes3(2, 4)
    real(c_double) :: point1(2), point2(3), point3(4)
    logical(c_bool) :: is_contained
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "contains_nd"

    ! CASE 1: Below bounding box.
    nodes1(1, :) = [0.0_dp, 1.0_dp]
    nodes1(2, :) = [0.5_dp, 0.0_dp]
    nodes1(3, :) = [1.0_dp, 2.0_dp]
    point1 = [-0.5_dp, 1.0_dp]
    call contains_nd(3, 2, nodes1, point1, is_contained)
    if (.NOT. is_contained) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Above bounding box.
    nodes2(1, :) = [0.0_dp, -4.0_dp, 2.0_dp]
    nodes2(2, :) = [-1.0_dp, 1.0_dp, 3.0_dp]
    point2 = [-0.5_dp, 2.0_dp, 2.5_dp]
    call contains_nd(2, 3, nodes2, point2, is_contained)
    if (.NOT. is_contained) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: Inside bounding box.
    nodes3(1, :) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    nodes3(2, :) = [1.0_dp, -2.0_dp, -4.0_dp, 1.0_dp]
    point3 = [0.5_dp, 0.0_dp, 0.0_dp, 2.0_dp]
    call contains_nd(2, 4, nodes3, point3, is_contained)
    if (is_contained) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_contains_nd

  subroutine test_vector_close(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: eps
    real(c_double) :: vec1(1, 2)
    real(c_double) :: vec2(1, 2)
    logical(c_bool) :: is_close
    integer :: case_id
    character(:), allocatable :: name

    eps = 0.5_dp**40
    case_id = 1
    name = "vector_close"

    ! CASE 1: Identical vector.
    vec1(1, :) = [0.5_dp, 4.0_dp]
    is_close = vector_close(2, vec1, vec1, eps)
    if (is_close) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Far apart vectors.
    vec1(1, :) = [0.0_dp, 6.0_dp]
    vec2(1, :) = [1.0_dp, -4.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (.NOT. is_close) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: Close but different.
    vec1(1, :) = [2.25_dp, -3.5_dp]
    vec2(1, :) = vec1(1, :) + 0.5_dp**43 * [-5.0_dp, 12.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (is_close) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 4: Custom epsilon.
    vec1(1, :) = [3.0_dp, 4.0_dp]
    vec2(1, :) = [2.0_dp, 5.0_dp]
    is_close = vector_close(2, vec1, vec2, 0.5_dp)
    if (is_close .AND. .NOT. vector_close(2, vec1, vec2, eps)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 5: Near zero.
    vec1(1, :) = [0.0_dp, 0.0_dp]
    vec2(1, :) = 0.5_dp**45 * [3.0_dp, 4.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (is_close) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 6: Near zero failure (i.e. not near enough).
    vec1(1, :) = 0.5_dp**20 * [1.0_dp, 0.0_dp]
    vec2(1, :) = [0.0_dp, 0.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (.NOT. is_close) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_vector_close

  subroutine test_in_interval(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical(c_bool) :: is_inside
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "in_interval"

    ! CASE 1: Interior value.
    is_inside = in_interval(1.5_dp, 1.0_dp, 2.0_dp)
    if (is_inside) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Barely inside.
    is_inside = in_interval(1.0_dp + MACHINE_EPS, 1.0_dp, 2.0_dp)
    if (is_inside) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: Barely outside.
    is_inside = in_interval(1.0_dp - MACHINE_EPS / 2, 1.0_dp, 2.0_dp)
    if (.NOT. is_inside) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 4: Exterior value.
    is_inside = in_interval(-1.0_dp, 1.0_dp, 2.0_dp)
    if (.NOT. is_inside) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_in_interval

  subroutine test_ulps_away(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: eps, value1, value2
    logical(c_bool) :: is_near
    integer :: case_id
    character(:), allocatable :: name

    eps = 0.5_dp**40
    case_id = 1
    name = "ulps_away"

    ! CASE 1: First value is zero.
    is_near = ulps_away(0.0_dp, 0.0_dp, 1, eps)
    if (is_near) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Second value is zero.
    is_near = ulps_away(1.0_dp, 0.0_dp, 1, eps)
    if (.NOT. is_near) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: First value positive **AND** negative.
    is_near = ulps_away(1.0_dp, 1.0_dp + MACHINE_EPS, 1, eps)
    if (is_near .AND. ulps_away(-1.0_dp, -1.0_dp - MACHINE_EPS, 1, eps)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 4: Boundaries where a single bit works.
    if ( &
         ulps_away(1.0_dp, 1.0_dp, 1, eps) .AND. &
         ulps_away(1.0_dp, 1.0_dp + MACHINE_EPS, 1, eps) .AND. &
         ulps_away(1.0_dp + MACHINE_EPS, 1.0_dp, 1, eps) .AND. &
         ulps_away(1.0_dp, 1.0_dp - MACHINE_EPS / 2, 1, eps) .AND. &
         ulps_away(1.0_dp - MACHINE_EPS / 2, 1.0_dp, 1, eps)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 5: Non-default ``num_bits``.
    value1 = 1.5_dp
    value2 = value1 + 0.5_dp**43
    if ( &
         .NOT. ulps_away(value1, value2, 1, eps) .AND. &
         ulps_away(value1, value2, 1000, eps)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 6: Very close, but not close enough.
    value1 = 0.25_dp - 5.625_dp * MACHINE_EPS
    value2 = 0.25_dp - 6.25_dp * MACHINE_EPS
    if ( &
         .NOT. ulps_away(value1, value2, 1, eps) .AND. &
         .NOT. ulps_away(value1, value2, 4, eps) .AND. &
         ulps_away(value1, value2, 5, eps)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_ulps_away

end module test_helpers
