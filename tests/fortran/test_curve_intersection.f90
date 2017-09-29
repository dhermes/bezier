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

  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use, intrinsic :: iso_c_binding, only: c_bool, c_double, c_int
  use curve, only: evaluate_multi, subdivide_nodes
  use curve_intersection, only: &
       BoxIntersectionType_INTERSECTION, BoxIntersectionType_TANGENT, &
       BoxIntersectionType_DISJOINT, FROM_LINEARIZED_SUCCESS, &
       FROM_LINEARIZED_PARALLEL, linearization_error, segment_intersection, &
       newton_refine_intersect, bbox_intersect, parallel_different, &
       from_linearized, bbox_line_intersect
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private &
       test_linearization_error, test_segment_intersection, &
       test_newton_refine_intersect, test_bbox_intersect, &
       test_parallel_different, test_from_linearized, test_bbox_line_intersect
  public curve_intersection_all_tests

contains

  subroutine curve_intersection_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_linearization_error(success)
    call test_segment_intersection(success)
    call test_newton_refine_intersect(success)
    call test_bbox_intersect(success)
    call test_parallel_different(success)
    call test_from_linearized(success)
    call test_bbox_line_intersect(success)

  end subroutine curve_intersection_all_tests

  subroutine test_linearization_error(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
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
         2, 2, nodes1, error1)
    case_success = (error1 == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Degree-elevated line (i.e. no error) as quadratic.
    nodes2(1, :) = 0
    nodes2(2, :) = [0.5_dp, 1.0_dp]
    nodes2(3, :) = [1.0_dp, 2.0_dp]
    call linearization_error( &
         3, 2, nodes2, error1)
    case_success = (error1 == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Degree-elevated line (i.e. no error) as quartic.
    nodes3(1, :) = 0
    nodes3(2, :) = [0.25_dp, 0.5_dp]
    nodes3(3, :) = [0.5_dp, 1.0_dp]
    nodes3(4, :) = [0.75_dp, 1.5_dp]
    nodes3(5, :) = [1.0_dp, 2.0_dp]
    call linearization_error( &
         5, 2, nodes3, error1)
    case_success = (error1 == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Line with bad parameterization.
    ! NOTE: This is the line 3 y = 4 x, but with the parameterization
    !       x(s) = 3 s (4 - 3 s).
    nodes2(1, :) = 0
    nodes2(2, :) = [6.0_dp, 8.0_dp]
    nodes2(3, :) = [3.0_dp, 4.0_dp]
    call linearization_error( &
         3, 2, nodes2, error1)
    ! D^2 v = [-9, -12]
    expected = 0.125_dp * 2 * 1 * 15.0_dp
    case_success = (error1 == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Quadratic curve.
    nodes2(1, :) = 0
    nodes2(2, :) = [1.0_dp, 1.0_dp]
    nodes2(3, :) = [5.0_dp, 6.0_dp]
    call linearization_error( &
         3, 2, nodes2, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [1, 1], [4, 5]
    !           d^2 Nodes = [3, 4]
    !       so that sqrt(3^2 + 4^2) = 5.0
    expected = 0.125_dp * 2 * 1 * 5.0_dp
    case_success = (error1 == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Subdivided curves (left and right) from CASE 5.
    call subdivide_nodes( &
         3, 2, nodes2, left_nodes2, right_nodes2)
    call linearization_error( &
         3, 2, left_nodes2, error1)
    call linearization_error( &
         3, 2, right_nodes2, error2)
    ! For a degree two curve, the 2nd derivative is constant
    ! so by subdividing, our error should drop by a factor
    ! of (1/2)^2 = 4.
    expected = 0.25_dp * expected
    case_success = (error1 == expected .AND. error2 == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Quadratic curve in 3D.
    nodes4(1, :) = [1.5_dp, 0.0_dp, 6.25_dp]
    nodes4(2, :) = [3.5_dp, -5.0_dp, 10.25_dp]
    nodes4(3, :) = [8.5_dp, 2.0_dp, 10.25_dp]
    call linearization_error( &
         3, 3, nodes4, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [2, -5, 4], [5, 7, 0]
    !           d^2 Nodes = [3, 12, -4]
    !       so that sqrt(3^2 + 12^2 + 4^2) = 13.0
    expected = 0.125_dp * 2 * 1 * 13.0_dp
    case_success = (error1 == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Quadratic with bad parameterization.
    ! NOTE: This is the quadratic y = 1 + x^2 / 4, but with the
    !       parameterization x(s) = (3 s - 1)^2.
    nodes3(1, :) = [1.0_dp, 1.25_dp]
    nodes3(2, :) = [-0.5_dp, 0.5_dp]
    nodes3(3, :) = [-0.5_dp, 2.0_dp]
    nodes3(4, :) = [1.0_dp, -1.0_dp]
    nodes3(5, :) = [4.0_dp, 5.0_dp]
    call linearization_error( &
         5, 2, nodes3, error1)
    ! D^2 v = [1.5, 2.25], [1.5, -4.5], [1.5, 9]
    expected = 0.125_dp * 4 * 3 * sqrt(1.5_dp**2 + 9.0_dp**2)
    case_success = (abs(error1 - expected) <= spacing(expected))
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Cubic curve.
    nodes5(1, :) = 0
    nodes5(2, :) = 1
    nodes5(3, :) = [5.0_dp, 6.0_dp]
    nodes5(4, :) = [6.0_dp, 7.0_dp]
    call linearization_error( &
         4, 2, nodes5, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [1, 1], [4, 5], [1, 1]
    !           d^2 Nodes = [3, 4], [-3, -4]
    !       so that sqrt(3^2 + 4^2) = 5.0
    expected = 0.125_dp * 3 * 2 * 5.0_dp
    case_success = (error1 == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: Quartic curve.
    nodes3(1, :) = 0
    nodes3(2, :) = 1
    nodes3(3, :) = [5.0_dp, 6.0_dp]
    nodes3(4, :) = [6.0_dp, 7.0_dp]
    nodes3(5, :) = [4.0_dp, 7.0_dp]
    call linearization_error( &
         5, 2, nodes3, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [1, 1], [4, 5], [1, 1], [-2, 0]
    !           d^2 Nodes = [3, 4], [-3, -4], [-3, -1]
    !       so that sqrt(3^2 + 4^2) = 5.0
    expected = 0.125_dp * 4 * 3 * 5.0_dp
    case_success = (error1 == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 11: Quintic curve (i.e. degree 5).
    nodes6(1, :) = [0.0_dp, 0.0_dp]
    nodes6(2, :) = [1.0_dp, 1.0_dp]
    nodes6(3, :) = [7.0_dp, 3.0_dp]
    nodes6(4, :) = [11.0_dp, 8.0_dp]
    nodes6(5, :) = [15.0_dp, 1.0_dp]
    nodes6(6, :) = [16.0_dp, -3.0_dp]
    call linearization_error( &
         6, 2, nodes6, error1)
    ! NOTE: This is hand picked so that
    !             d Nodes = [1, 1], [6, 2], [4, 5], [4, -7], [1, -4]
    !           d^2 Nodes = [5, 1], [-2, 3], [0, -12], [-3, 3]
    !       so that sqrt(5^2 + 12^2) = 13.0
    expected = 0.125_dp * 5 * 4 * 13.0_dp
    case_success = (error1 == expected)
    call print_status(name, case_id, case_success, success)

  end subroutine test_linearization_error

  subroutine test_segment_intersection(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
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
    case_success = (si_success .AND. s == 0.25_dp .AND. t == 0.625_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Parallel segments.
    start0(1, :) = [0.0_dp, 0.5_dp]
    end0(1, :) = [0.0_dp, -0.5_dp]
    start1(1, :) = [0.0_dp, 1.0_dp]
    end1(1, :) = [0.0_dp, -1.0_dp]
    call segment_intersection( &
         start0, end0, start1, end1, s, t, si_success)
    case_success = (.NOT. si_success)
    call print_status(name, case_id, case_success, success)

  end subroutine test_segment_intersection

  subroutine test_newton_refine_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 2), nodes2(2, 2)
    real(c_double) :: nodes3(3, 2), nodes4(3, 2), nodes5(5, 2)
    real(c_double) :: known_s, known_t
    real(c_double) :: wrong_s, wrong_t
    real(c_double) :: new_s, new_t
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
         wrong_s, 2, nodes1, wrong_t, 2, nodes2, new_s, new_t)
    case_success = ( &
         new_s == known_s .AND. new_t == known_t .AND. &
         curves_intersect(nodes1, known_s, nodes2, known_t))
    call print_status(name, case_id, case_success, success)

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
         wrong_s, 3, nodes3, wrong_t, 2, nodes1, new_s, new_t)
    case_success = ( &
         new_s == 0.4375_dp .AND. &
         abs(known_s - new_s) < abs(known_s - wrong_s) .AND. &
         new_t == 0.5625_dp .AND. &
         abs(known_t - new_t) < abs(known_t - wrong_t) .AND. &
         curves_intersect(nodes3, known_s, nodes1, known_t))
    call print_status(name, case_id, case_success, success)

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
         known_s, 3, nodes3, known_t, 3, nodes4, new_s, new_t)
    case_success = ( &
         new_s == known_s .AND. new_t == known_t .AND. &
         curves_intersect(nodes3, known_s, nodes4, known_t))
    call print_status(name, case_id, case_success, success)

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
         wrong_s, 3, nodes3, wrong_t, 3, nodes4, new_s, new_t)
    case_success = ( &
         new_s == 0.2421875_dp .AND. &
         abs(known_s - new_s) < abs(known_s - wrong_s) .AND. &
         new_t == 0.7578125_dp .AND. &
         abs(known_t - new_t) < abs(known_t - wrong_t) .AND. &
         curves_intersect(nodes3, known_s, nodes4, known_t))
    call print_status(name, case_id, case_success, success)

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
    case_success = .TRUE.
    do i = 1, 3
       call newton_refine_intersect( &
            wrong_s, 5, nodes5, wrong_t, 2, nodes1, new_s, new_t)
       wrong_s = new_s
       wrong_t = new_t
       if (i == 1) then
          case_success = ( &
               case_success .AND. &
               new_s == known_s .AND. &
               new_t == 2.0_dp)
       else
          case_success = ( &
               case_success .AND. &
               new_s == known_s .AND. &
               new_t == known_t)
       end if
    end do

    case_success = ( &
         case_success .AND. &
         curves_intersect(nodes5, known_s, nodes1, known_t))
    call print_status(name, case_id, case_success, success)

  end subroutine test_newton_refine_intersect

  function curves_intersect(nodes1, s, nodes2, t) result(predicate)
    real(c_double), intent(in) :: nodes1(:, :), nodes2(:, :)
    real(c_double), intent(in) :: s, t
    logical(c_bool) :: predicate
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: point1(1, 2), point2(1, 2)

    call evaluate_multi( &
         size(nodes1, 1), size(nodes1, 2), nodes1, 1, [s], point1)
    call evaluate_multi( &
         size(nodes2, 1), size(nodes2, 2), nodes2, 1, [t], point2)
    predicate = all(point1 == point2)

  end function curves_intersect

  subroutine test_bbox_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: unit_square(4, 2), other(4, 2), delta
    integer(c_int) :: enum_, i
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "bbox_intersect"
    unit_square(1, :) = [0.0_dp, 0.0_dp]
    unit_square(2, :) = [1.0_dp, 0.0_dp]
    unit_square(3, :) = [1.0_dp, 1.0_dp]
    unit_square(4, :) = [0.0_dp, 1.0_dp]

    ! CASE 1: Intersecting bbox-es.
    forall (i = 1:4)
       other(i, :) = unit_square(i, :) + [0.5_dp, 0.5_dp]
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Far apart bbox-es.
    forall (i = 1:4)
       other(i, :) = unit_square(i, :) + [100.0_dp, 100.0_dp]
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_DISJOINT)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Disjoint bbox-es that have an "aligned" edge.
    forall (i = 1:4)
       other(i, :) = unit_square(i, :) + [1.0_dp, 2.0_dp]
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_DISJOINT)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Tangent bbox-es.
    forall (i = 1:4)
       other(i, :) = unit_square(i, :) + [1.0_dp, 0.0_dp]
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_TANGENT)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Almost tangent bbox-es.
    delta = 1.0_dp + spacing(1.0_dp)
    forall (i = 1:4)
       other(i, :) = unit_square(i, :) + [delta, 0.0_dp]
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_DISJOINT)
    call print_status(name, case_id, case_success, success)

  end subroutine test_bbox_intersect

  subroutine test_parallel_different(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    logical(c_bool) :: result_
    real(c_double) :: start0(1, 2), end0(1, 2)
    real(c_double) :: start1(1, 2), end1(1, 2)
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "parallel_different"

    ! CASE 1: Same line, but segments don't overlap.
    start0(1, :) = 0
    end0(1, :) = [3.0_dp, 4.0_dp]
    start1(1, :) = [6.0_dp, 8.0_dp]
    end1(1, :) = [9.0_dp, 12.0_dp]
    call parallel_different( &
         start0, end0, start1, end1, result_)
    case_success = (result_)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Same line, ``start1`` contained in segment "0".
    ! NOTE: All of segment "1" is contained in segment "0", but the
    !       computation stops after ``start1`` is found on the segment.
    start0(1, :) = [6.0_dp, -3.0_dp]
    end0(1, :) = [-7.0_dp, 1.0_dp]
    start1(1, :) = [1.125_dp, -1.5_dp]
    end1(1, :) = [-5.375_dp, 0.5_dp]
    call parallel_different( &
         start0, end0, start1, end1, result_)
    case_success = (.NOT. result_)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Same line, ``end1`` contained in segment "0".
    start0(1, :) = [1.0_dp, 2.0_dp]
    end0(1, :) = [3.0_dp, 5.0_dp]
    start1(1, :) = [-0.5_dp, -0.25_dp]
    end1(1, :) = [2.0_dp, 3.5_dp]
    call parallel_different( &
         start0, end0, start1, end1, result_)
    case_success = (.NOT. result_)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Same line, segment "0" fully contained in segment "1".
    start0(1, :) = [-9.0_dp, 0.0_dp]
    end0(1, :) = [4.0_dp, 5.0_dp]
    start1(1, :) = [23.5_dp, 12.5_dp]
    end1(1, :) = [-25.25_dp, -6.25_dp]
    call parallel_different( &
         start0, end0, start1, end1, result_)
    case_success = (.NOT. result_)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Parallel, but different lines.
    start0(1, :) = [3.0_dp, 2.0_dp]
    end0(1, :) = [3.0_dp, 0.75_dp]
    start1(1, :) = 0
    end1(1, :) = [0.0_dp, 2.0_dp]
    call parallel_different( &
         start0, end0, start1, end1, result_)
    case_success = (result_)
    call print_status(name, case_id, case_success, success)

  end subroutine test_parallel_different

  subroutine test_from_linearized(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nan_val
    real(c_double) :: start_node1(1, 2), end_node1(1, 2), error1
    real(c_double) :: start_node2(1, 2), end_node2(1, 2), error2
    real(c_double) :: nodes1(3, 2), nodes2(3, 2)
    real(c_double) :: nodes3(2, 2), nodes4(2, 2)
    real(c_double) :: refined_s, refined_t
    logical(c_bool) :: does_intersect
    integer(c_int) :: py_exc
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "from_linearized"
    nan_val = ieee_value(nan_val, ieee_quiet_nan)

    ! CASE 1: Basic test of not-very-linearized quadratics.
    nodes1(1, :) = 0
    nodes1(2, :) = [0.5_dp, 1.0_dp]
    nodes1(3, :) = 1
    start_node1 = nodes1(:1, :)
    end_node1 = nodes1(3:, :)
    ! NOTE: This curve isn't close to linear, but that's OK.
    error1 = nan_val

    nodes2(1, :) = [0.0_dp, 1.0_dp]
    nodes2(2, :) = [0.5_dp, 1.0_dp]
    nodes2(3, :) = [1.0_dp, 0.0_dp]
    start_node2 = nodes2(:1, :)
    end_node2 = nodes2(3:, :)
    ! NOTE: This curve isn't close to linear, but that's OK.
    error2 = nan_val

    call from_linearized( &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 3, nodes1, &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 3, nodes2, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = ( &
         does_intersect .AND. py_exc == FROM_LINEARIZED_SUCCESS .AND. &
         refined_s == 0.5_dp .AND. refined_t == 0.5_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Bounding boxes intersect but the lines do not.
    nodes3(1, :) = 0
    nodes3(2, :) = 1
    start_node1 = nodes3(:1, :)
    end_node1 = nodes3(2:, :)
    error1 = 0.0_dp

    nodes4(1, :) = [1.75_dp, -0.75_dp]
    nodes4(2, :) = [0.75_dp, 0.25_dp]
    start_node2 = nodes4(:1, :)
    end_node2 = nodes4(2:, :)
    error2 = 0.0_dp

    call from_linearized( &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 2, nodes3, &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 2, nodes4, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = ( &
         .NOT. does_intersect .AND. py_exc == FROM_LINEARIZED_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Same as CASE 2, but swap the inputs.
    call from_linearized( &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 2, nodes4, &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 2, nodes3, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = ( &
         .NOT. does_intersect .AND. py_exc == FROM_LINEARIZED_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Bounding boxes intersect but the quadratics do not.
    nodes1(1, :) = 0
    nodes1(2, :) = [0.5_dp, 0.0_dp]
    nodes1(3, :) = 1
    start_node1 = nodes1(:1, :)
    end_node1 = nodes1(3:, :)
    error1 = 0.25_dp

    nodes2(1, :) = [1.75_dp, -0.75_dp]
    nodes2(2, :) = [1.25_dp, -0.75_dp]
    nodes2(3, :) = [0.75_dp, 0.25_dp]
    start_node2 = nodes2(:1, :)
    end_node2 = nodes2(3:, :)
    error2 = 0.25_dp

    call from_linearized( &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 3, nodes1, &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 3, nodes2, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = ( &
         .NOT. does_intersect .AND. py_exc == FROM_LINEARIZED_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Same as CASE 4, but swap the inputs.
    call from_linearized( &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 3, nodes2, &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 3, nodes1, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = ( &
         .NOT. does_intersect .AND. py_exc == FROM_LINEARIZED_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Parallel lines that do not intersect.
    nodes3(1, :) = 0
    nodes3(2, :) = 1
    start_node1 = nodes3(:1, :)
    end_node1 = nodes3(2:, :)
    error1 = 0.0_dp

    nodes4(1, :) = [0.0_dp, 1.0_dp]
    nodes4(2, :) = [1.0_dp, 2.0_dp]
    start_node2 = nodes4(:1, :)
    end_node2 = nodes4(2:, :)
    error2 = 0.0_dp

    call from_linearized( &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 2, nodes3, &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 2, nodes4, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = ( &
         .NOT. does_intersect .AND. py_exc == FROM_LINEARIZED_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Parallel lines that **do** intersect.
    nodes3(1, :) = 0
    nodes3(2, :) = 1
    start_node1 = nodes3(:1, :)
    end_node1 = nodes3(2:, :)
    error1 = 0.0_dp

    nodes4(1, :) = [0.5_dp, 0.5_dp]
    nodes4(2, :) = [3.0_dp, 3.0_dp]
    start_node2 = nodes4(:1, :)
    end_node2 = nodes4(2:, :)
    error2 = 0.0_dp

    call from_linearized( &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 2, nodes3, &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 2, nodes4, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = (py_exc == FROM_LINEARIZED_PARALLEL)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Linearized parts are same line but disjoint segments.
    nodes3(1, :) = 0
    nodes3(2, :) = 1
    start_node1 = nodes3(:1, :)
    end_node1 = nodes3(2:, :)
    error1 = 0.0_dp

    nodes2(1, :) = 2
    nodes2(2, :) = 2.5009765625_dp
    nodes2(3, :) = 3
    start_node2 = nodes2(:1, :)
    end_node2 = nodes2(3:, :)
    error2 = nan_val

    call from_linearized( &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 2, nodes3, &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 3, nodes2, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = ( &
         .NOT. does_intersect .AND. py_exc == FROM_LINEARIZED_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Linearized parts parallel / diff. lines / bbox-es overlap.
    nodes3(1, :) = 0
    nodes3(2, :) = 1
    start_node1 = nodes3(:1, :)
    end_node1 = nodes3(2:, :)
    error1 = 0.0_dp

    nodes2(1, :) = [0.5_dp, 0.75_dp]
    nodes2(2, :) = [1.0009765625_dp, 1.2509765625_dp]
    nodes2(3, :) = [1.5_dp, 1.75_dp]
    start_node2 = nodes2(:1, :)
    end_node2 = nodes2(3:, :)
    error2 = nan_val

    call from_linearized( &
         error1, 0.0_dp, 1.0_dp, start_node1, end_node1, 2, nodes3, &
         error2, 0.0_dp, 1.0_dp, start_node2, end_node2, 3, nodes2, &
         refined_s, refined_t, does_intersect, py_exc)
    case_success = (py_exc == FROM_LINEARIZED_PARALLEL)
    call print_status(name, case_id, case_success, success)

  end subroutine test_from_linearized

  subroutine test_bbox_line_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: unit_square(4, 2)
    real(c_double) :: line_start(1, 2), line_end(1, 2)
    integer(c_int) :: enum_
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "bbox_line_intersect"
    unit_square(1, :) = [0.0_dp, 0.0_dp]
    unit_square(2, :) = [1.0_dp, 0.0_dp]
    unit_square(3, :) = [1.0_dp, 1.0_dp]
    unit_square(4, :) = [0.0_dp, 1.0_dp]

    ! CASE 1: Line starts inside bounding box.
    line_start(1, :) = 0.5_dp
    line_end(1, :) = [0.5_dp, 1.5_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Line ends (but does not start) in bounding box.
    line_start(1, :) = [-1.0_dp, 0.5_dp]
    line_end(1, :) = 0.5_dp
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Line segment "pierces" bbox from the bottom.
    line_start(1, :) = [0.5_dp, -0.5_dp]
    line_end(1, :) = [0.5_dp, 1.5_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Line segment "pierces" bbox from the right.
    line_start(1, :) = [-0.5_dp, 0.5_dp]
    line_end(1, :) = [1.5_dp, 0.5_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Line segment "pierces" bbox from the top.
    line_start(1, :) = [-0.25_dp, 0.5_dp]
    line_end(1, :) = [0.5_dp, 1.25_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Line segment is disjoint from bbox.
    line_start(1, :) = 2
    line_end(1, :) = [2.0_dp, 5.0_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_DISJOINT)
    call print_status(name, case_id, case_success, success)

  end subroutine test_bbox_line_intersect

end module test_curve_intersection
