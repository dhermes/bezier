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

  use iso_c_binding, only: c_bool, c_double, c_int
  use curve, only: evaluate_multi, subdivide_nodes
  use curve_intersection, only: &
       BoxIntersectionType_INTERSECTION, BoxIntersectionType_TANGENT, &
       BoxIntersectionType_DISJOINT, linearization_error, &
       segment_intersection, newton_refine_intersect, &
       bbox_intersect
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private &
       test_linearization_error, test_segment_intersection, &
       test_newton_refine_intersect, test_bbox_intersect
  public curve_intersection_all_tests

contains

  subroutine curve_intersection_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_linearization_error(success)
    call test_segment_intersection(success)
    call test_newton_refine_intersect(success)
    call test_bbox_intersect(success)

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

  subroutine test_segment_intersection(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: start0(1, 2), end0(1, 2)
    real(c_double) :: start1(1, 2), end1(1, 2)
    real(c_double) :: s, t
    logical(c_bool) :: si_success
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "segment_intersection"

    ! CASE 1: Segments that intersect.
    start0(1, :) = [1.75_dp, 2.125_dp]
    end0(1, :) = [-1.25_dp, 1.625_dp]
    start1(1, :) = [-0.25_dp, 2.625_dp]
    end1(1, :) = [1.75_dp, 1.625_dp]
    ! D0 x D1 == 4.0, so there will be no round-off in answer.
    call segment_intersection( &
         start0, end0, start1, end1, s, t, si_success)
    if (si_success .AND. s == 0.25_dp .AND. t == 0.625_dp) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Parallel segments.
    start0(1, :) = [0.0_dp, 0.5_dp]
    end0(1, :) = [0.0_dp, -0.5_dp]
    start1(1, :) = [0.0_dp, 1.0_dp]
    end1(1, :) = [0.0_dp, -1.0_dp]
    call segment_intersection( &
         start0, end0, start1, end1, s, t, si_success)
    if (.NOT. si_success) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_segment_intersection

  subroutine test_newton_refine_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes1(2, 2), nodes2(2, 2)
    real(c_double) :: nodes3(3, 2), nodes4(3, 2), nodes5(5, 2)
    real(c_double) :: known_s, known_t
    real(c_double) :: wrong_s, wrong_t
    real(c_double) :: new_s, new_t
    logical(c_bool) :: local_success
    integer :: i
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "newton_refine_intersect"

    ! CASE 1: Intersection of two lines.
    ! Newton's method is exact on linear problems so will
    ! always converge after one step.
    nodes1(1, :) = 0
    nodes1(2, :) = 1
    known_s = 0.75_dp
    nodes2(1, :) = [1.0_dp, 0.0_dp]
    nodes2(2, :) = [0.0_dp, 3.0_dp]
    known_t = 0.25_dp
    ! NOTE: By construction, the Jacobian matrix will be
    !           [1, 1], [1, -3]
    !       which has determinant -4.0, hence there will
    !       be no round-off when solving.
    wrong_s = known_s - 0.125_dp
    wrong_t = known_t + 0.125_dp
    call newton_refine_intersect( &
         wrong_s, 1, nodes1, wrong_t, 1, nodes2, new_s, new_t)
    if (new_s == known_s .AND. new_t == known_t .AND. &
         curves_intersect(nodes1, known_s, nodes2, known_t)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Mixed degree (line and quadratic).
    nodes3(1, :) = 0
    nodes3(2, :) = [0.5_dp, 1.0_dp]
    nodes3(3, :) = [1.0_dp, 0.0_dp]
    known_s = 0.5_dp
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [0.0_dp, 1.0_dp]
    known_t = 0.5_dp
    ! NOTE: By construction, the Jacobian matrix will be
    !           [1, 1], [1, -1]
    !       which has determinant -2.0, hence there will
    !       be no round-off when solving.
    wrong_s = 0.25_dp
    wrong_t = 0.25_dp
    call newton_refine_intersect( &
         wrong_s, 2, nodes3, wrong_t, 1, nodes1, new_s, new_t)
    if (new_s == 0.4375_dp .AND. &
         abs(known_s - new_s) < abs(known_s - wrong_s) .AND. &
         new_t == 0.5625_dp .AND. &
         abs(known_t - new_t) < abs(known_t - wrong_t) .AND. &
         curves_intersect(nodes3, known_s, nodes1, known_t)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: Early exit (i.e. already an intersection).
    nodes3(1, :) = 0
    nodes3(2, :) = [0.5_dp, 1.0_dp]
    nodes3(3, :) = [1.0_dp, 0.0_dp]
    known_s = 0.25_dp
    nodes4(1, :) = [1.0_dp, 0.75_dp]
    nodes4(2, :) = [0.5_dp, -0.25_dp]
    nodes4(3, :) = [0.0_dp, 0.75_dp]
    known_t = 0.75_dp
    call newton_refine_intersect( &
         known_s, 2, nodes3, known_t, 2, nodes4, new_s, new_t)
    if (new_s == known_s .AND. new_t == known_t .AND. &
         curves_intersect(nodes3, known_s, nodes4, known_t)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 4: Intersect quadratics.
    nodes3(1, :) = 0
    nodes3(2, :) = [0.5_dp, 1.0_dp]
    nodes3(3, :) = [1.0_dp, 0.0_dp]
    known_s = 0.25_dp
    nodes4(1, :) = [1.0_dp, 0.75_dp]
    nodes4(2, :) = [0.5_dp, -0.25_dp]
    nodes4(3, :) = [0.0_dp, 0.75_dp]
    known_t = 0.75_dp
    ! NOTE: By construction, the Jacobian matrix will be
    !           [1, 3/4], [1, -5/4]
    !       which has determinant -2.0, hence there will
    !       be no round-off when solving.
    wrong_s = known_s + 0.0625_dp
    wrong_t = known_t + 0.0625_dp
    call newton_refine_intersect( &
         wrong_s, 2, nodes3, wrong_t, 2, nodes4, new_s, new_t)
    if (new_s == 0.2421875_dp .AND. &
         abs(known_s - new_s) < abs(known_s - wrong_s) .AND. &
         new_t == 0.7578125_dp .AND. &
         abs(known_t - new_t) < abs(known_t - wrong_t) .AND. &
         curves_intersect(nodes3, known_s, nodes4, known_t)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 5: Convergence test.
    nodes5(1, :) = [0.0_dp, 0.0_dp]
    nodes5(2, :) = [0.25_dp, 1.0_dp]
    nodes5(3, :) = [0.5_dp, -0.75_dp]
    nodes5(4, :) = [0.75_dp, 1.0_dp]
    nodes5(5, :) = [1.0_dp, 0.0_dp]
    known_s = 0.5_dp
    nodes1(1, :) = [0.5_dp, 0.0_dp]
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    known_t = 0.21875_dp
    ! We start with a very wrong guess and update it three times.
    wrong_s = 0.0_dp
    wrong_t = 0.0_dp
    local_success = .TRUE.
    do i = 1, 3
       call newton_refine_intersect( &
            wrong_s, 4, nodes5, wrong_t, 1, nodes1, new_s, new_t)
       wrong_s = new_s
       wrong_t = new_t
       if (i == 1) then
          if (new_s /= known_s .OR. new_t /= 2.0_dp) then
             local_success = .FALSE.
          end if
       else
          if (new_s /= known_s .OR. new_t /= known_t) then
             local_success = .FALSE.
          end if
       end if
    end do

    if (local_success .AND. &
         curves_intersect(nodes5, known_s, nodes1, known_t)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_newton_refine_intersect

  function curves_intersect(nodes1, s, nodes2, t) result(predicate)
    real(c_double), intent(in) :: nodes1(:, :), nodes2(:, :)
    real(c_double), intent(in) :: s, t
    logical(c_bool) :: predicate
    ! Variables outside of signature.
    real(c_double) :: point1(1, 2), point2(1, 2)

    call evaluate_multi( &
         size(nodes1, 1) - 1, size(nodes1, 2), nodes1, 1, [s], point1)
    call evaluate_multi( &
         size(nodes2, 1) - 1, size(nodes2, 2), nodes2, 1, [t], point2)
    predicate = all(point1 == point2)

  end function curves_intersect

  subroutine test_bbox_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes1(4, 2), nodes2(4, 2), delta
    integer(c_int) :: enum_
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "bbox_intersect"

    ! CASE 1: Intersecting bbox-es.
    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    nodes1(3, :) = [1.0_dp, 1.0_dp]
    nodes1(4, :) = [0.0_dp, 1.0_dp]
    nodes2(1, :) = [0.5_dp, 0.5_dp]
    nodes2(2, :) = [1.5_dp, 0.5_dp]
    nodes2(3, :) = [1.5_dp, 1.5_dp]
    nodes2(4, :) = [0.5_dp, 1.5_dp]
    call bbox_intersect( &
         4, nodes1, 4, nodes2, enum_)
    if (enum_ == BoxIntersectionType_INTERSECTION) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Far apart bbox-es.
    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    nodes1(3, :) = [1.0_dp, 1.0_dp]
    nodes1(4, :) = [0.0_dp, 1.0_dp]
    nodes2(1, :) = [100.0_dp, 100.0_dp]
    nodes2(2, :) = [101.0_dp, 100.0_dp]
    nodes2(3, :) = [101.0_dp, 101.0_dp]
    nodes2(4, :) = [100.0_dp, 101.0_dp]
    call bbox_intersect( &
         4, nodes1, 4, nodes2, enum_)
    if (enum_ == BoxIntersectionType_DISJOINT) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: Disjoint bbox-es that have an "aligned" edge.
    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    nodes1(3, :) = [1.0_dp, 1.0_dp]
    nodes1(4, :) = [0.0_dp, 1.0_dp]
    nodes2(1, :) = [1.0_dp, 2.0_dp]
    nodes2(2, :) = [2.0_dp, 2.0_dp]
    nodes2(3, :) = [2.0_dp, 3.0_dp]
    nodes2(4, :) = [1.0_dp, 3.0_dp]
    call bbox_intersect( &
         4, nodes1, 4, nodes2, enum_)
    if (enum_ == BoxIntersectionType_DISJOINT) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 4: Tangent bbox-es.
    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    nodes1(3, :) = [1.0_dp, 1.0_dp]
    nodes1(4, :) = [0.0_dp, 1.0_dp]
    nodes2(1, :) = [1.0_dp, 0.0_dp]
    nodes2(2, :) = [2.0_dp, 0.0_dp]
    nodes2(3, :) = [2.0_dp, 1.0_dp]
    nodes2(4, :) = [1.0_dp, 1.0_dp]
    call bbox_intersect( &
         4, nodes1, 4, nodes2, enum_)
    if (enum_ == BoxIntersectionType_TANGENT) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 5: Almost tangent bbox-es.
    delta = 1.0_dp + spacing(1.0_dp)
    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    nodes1(3, :) = [1.0_dp, 1.0_dp]
    nodes1(4, :) = [0.0_dp, 1.0_dp]
    nodes2(1, :) = [delta, 0.0_dp]
    nodes2(2, :) = [1.0_dp + delta, 0.0_dp]
    nodes2(3, :) = [1.0_dp + delta, 1.0_dp]
    nodes2(4, :) = [delta, 1.0_dp]
    call bbox_intersect( &
         4, nodes1, 4, nodes2, enum_)
    if (enum_ == BoxIntersectionType_DISJOINT) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_bbox_intersect

end module test_curve_intersection
