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

  use, intrinsic :: iso_c_binding, only: c_bool, c_double, c_int
  use status, only: &
       Status_SUCCESS, Status_PARALLEL, Status_WIGGLE_FAIL, &
       Status_NO_CONVERGE, Status_INSUFFICIENT_SPACE
  use curve, only: &
       CurveData, evaluate_multi, specialize_curve, subdivide_nodes, &
       curves_equal, subdivide_curve
  use curve_intersection, only: &
       BoxIntersectionType_INTERSECTION, BoxIntersectionType_TANGENT, &
       BoxIntersectionType_DISJOINT, &
       Subdivide_FIRST, Subdivide_SECOND, &
       Subdivide_BOTH, Subdivide_NEITHER, &
       linearization_error, segment_intersection, &
       newton_refine_intersect, bbox_intersect, parallel_different, &
       from_linearized, bbox_line_intersect, add_intersection, &
       add_from_linearized, endpoint_check, tangent_bbox_intersection, &
       add_candidates, intersect_one_round, all_intersections, &
       all_intersections_abi, set_max_candidates, get_max_candidates, &
       set_similar_ulps, get_similar_ulps
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private &
       test_linearization_error, test_segment_intersection, &
       test_newton_refine_intersect, test_bbox_intersect, &
       test_parallel_different, test_from_linearized, &
       test_bbox_line_intersect, test_add_intersection, &
       test_add_from_linearized, test_endpoint_check, &
       test_tangent_bbox_intersection, test_add_candidates, &
       test_intersect_one_round, test_all_intersections, &
       test_all_intersections_abi
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
    call test_add_intersection(success)
    call test_add_from_linearized(success)
    call test_endpoint_check(success)
    call test_tangent_bbox_intersection(success)
    call test_add_candidates(success)
    call test_intersect_one_round(success)
    call test_all_intersections(success)
    call test_all_intersections_abi(success)
    call test_set_max_candidates(success)
    call test_set_similar_ulps(success)

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
    character(19) :: name

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
    character(20) :: name

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
    character(23) :: name

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
    character(14) :: name

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
    character(18) :: name

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
    type(CurveData) :: curve1, curve2, curve3
    real(c_double) :: cubic(4, 2)
    real(c_double) :: error1, error2
    real(c_double) :: refined_s, refined_t
    logical(c_bool) :: does_intersect
    integer(c_int) :: status
    integer :: case_id
    character(15) :: name

    case_id = 1
    name = "from_linearized"

    ! Will hold throughout.
    curve1%start = 0.0_dp
    curve1%end_ = 1.0_dp
    curve2%start = 0.0_dp
    curve2%end_ = 1.0_dp

    ! CASE 1: Basic test of not-very-linearized quadratics.
    allocate(curve1%nodes(3, 2))
    curve1%nodes(1, :) = 0
    curve1%nodes(2, :) = [0.5_dp, 1.0_dp]
    curve1%nodes(3, :) = 1
    ! NOTE: This curve isn't close to linear, but that's OK.
    call linearization_error( &
         3, 2, curve1%nodes, error1)

    allocate(curve2%nodes(3, 2))
    curve2%nodes(1, :) = [0.0_dp, 1.0_dp]
    curve2%nodes(2, :) = [0.5_dp, 1.0_dp]
    curve2%nodes(3, :) = [1.0_dp, 0.0_dp]
    ! NOTE: This curve isn't close to linear, but that's OK.
    call linearization_error( &
         3, 2, curve2%nodes, error2)

    call from_linearized( &
         error1, curve1, 3, curve1%nodes, &
         error2, curve2, 3, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         does_intersect .AND. status == Status_SUCCESS .AND. &
         refined_s == 0.5_dp .AND. refined_t == 0.5_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Bounding boxes intersect but the quadratics do not.
    curve1%nodes(1, :) = 0
    curve1%nodes(2, :) = [0.5_dp, 0.0_dp]
    curve1%nodes(3, :) = 1
    error1 = 0.25_dp

    deallocate(curve2%nodes)
    allocate(curve2%nodes(3, 2))
    curve2%nodes(1, :) = [1.75_dp, -0.75_dp]
    curve2%nodes(2, :) = [1.25_dp, -0.75_dp]
    curve2%nodes(3, :) = [0.75_dp, 0.25_dp]
    error2 = 0.25_dp

    call from_linearized( &
         error1, curve1, 3, curve1%nodes, &
         error2, curve2, 3, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         .NOT. does_intersect .AND. status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Same as CASE 2, but swap the inputs.
    call from_linearized( &
         error2, curve2, 3, curve2%nodes, &
         error1, curve1, 3, curve1%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         .NOT. does_intersect .AND. status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Linearized parts are same line but disjoint segments.
    deallocate(curve1%nodes)
    allocate(curve1%nodes(2, 2))
    curve1%nodes(1, :) = 0
    curve1%nodes(2, :) = 1
    error1 = 0.0_dp

    curve2%nodes(1, :) = 2
    curve2%nodes(2, :) = 2.5009765625_dp
    curve2%nodes(3, :) = 3
    call linearization_error( &
         3, 2, curve2%nodes, error2)

    call from_linearized( &
         error1, curve1, 2, curve1%nodes, &
         error2, curve2, 3, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         .NOT. does_intersect .AND. status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Linearized parts parallel / diff. lines / bbox-es overlap.
    curve1%nodes(1, :) = 0
    curve1%nodes(2, :) = 1
    error1 = 0.0_dp

    curve2%nodes(1, :) = [0.5_dp, 0.75_dp]
    curve2%nodes(2, :) = [1.0009765625_dp, 1.2509765625_dp]
    curve2%nodes(3, :) = [1.5_dp, 1.75_dp]
    call linearization_error( &
         3, 2, curve2%nodes, error2)

    call from_linearized( &
         error1, curve1, 2, curve1%nodes, &
         error2, curve2, 3, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = (status == Status_PARALLEL)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Bounding boxes intersect but the lines do not.
    curve1%nodes(1, :) = 0
    curve1%nodes(2, :) = 1
    error1 = 0.0_dp

    deallocate(curve2%nodes)
    allocate(curve2%nodes(2, 2))
    curve2%nodes(1, :) = [1.75_dp, -0.75_dp]
    curve2%nodes(2, :) = [0.75_dp, 0.25_dp]
    error2 = 0.0_dp

    call from_linearized( &
         error1, curve1, 2, curve1%nodes, &
         error2, curve2, 2, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         .NOT. does_intersect .AND. status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Same as CASE 6, but swap the inputs.
    call from_linearized( &
         error2, curve2, 2, curve2%nodes, &
         error1, curve1, 2, curve1%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         .NOT. does_intersect .AND. status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Parallel lines that do not intersect.
    curve1%nodes(1, :) = 0
    curve1%nodes(2, :) = 1
    error1 = 0.0_dp

    curve2%nodes(1, :) = [0.0_dp, 1.0_dp]
    curve2%nodes(2, :) = [1.0_dp, 2.0_dp]
    error2 = 0.0_dp

    call from_linearized( &
         error1, curve1, 2, curve1%nodes, &
         error2, curve2, 2, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         .NOT. does_intersect .AND. status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Parallel lines that **do** intersect.
    curve1%nodes(1, :) = 0
    curve1%nodes(2, :) = 1
    error1 = 0.0_dp

    curve2%nodes(1, :) = [0.5_dp, 0.5_dp]
    curve2%nodes(2, :) = [3.0_dp, 3.0_dp]
    error2 = 0.0_dp

    call from_linearized( &
         error1, curve1, 2, curve1%nodes, &
         error2, curve2, 2, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = (status == Status_PARALLEL)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: A "wiggle" failure caused by almost touching segments
    !          that produce local parameters within the [-2^{-16}, 1 + 2^{-16}]
    !          interval but produce a "global" parameter outside of the more
    !          narrow [-2^{-45}, 1 + 2^{-45}] (where ``WIGGLE == 2^{-45}`` in
    !          ``helpers.f90``).
    deallocate(curve1%nodes)
    allocate(curve1%nodes(4, 2))
    curve1%nodes(1, :) = [-0.7993236103108717_dp, -0.21683567278362156_dp]
    curve1%nodes(2, :) = [-0.8072986524226636_dp, -0.21898490744674426_dp]
    curve1%nodes(3, :) = [-0.8152736945344552_dp, -0.2211341421098668_dp]
    curve1%nodes(4, :) = [-0.8232487366462472_dp, -0.2232833767729893_dp]
    call linearization_error( &
         4, 2, curve1%nodes, error1)

    curve3%start = 0.99609375_dp
    curve3%end_ = 1.0_dp
    allocate(curve3%nodes(4, 2))
    cubic(1, :) = [-0.7838204403623438_dp, -0.25519640597397464_dp]
    cubic(2, :) = [-0.7894577677825452_dp, -0.24259531488131633_dp]
    cubic(3, :) = [-0.7946421067207265_dp, -0.22976394420044136_dp]
    cubic(4, :) = [-0.799367666650849_dp, -0.21671303774854855_dp]
    call specialize_curve( &
         4, 2, cubic, curve3%start, curve3%end_, curve3%nodes)
    call linearization_error( &
         4, 2, curve3%nodes, error2)

    call from_linearized( &
         error1, curve1, 4, curve1%nodes, &
         error2, curve3, 4, cubic, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         does_intersect .AND. &
         status == Status_WIGGLE_FAIL)
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
    character(19) :: name

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

  subroutine test_add_intersection(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    integer :: case_id
    character(16) :: name

    case_id = 1
    name = "add_intersection"
    num_intersections = 0

    ! CASE 1: Intersections is **not allocated**.
    case_success = .NOT. allocated(intersections)
    call add_intersection( &
         0.25_dp, 0.5_dp, num_intersections, intersections)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         intersections(1, 1) == 0.25_dp .AND. &
         intersections(2, 1) == 0.5_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Intersections is allocated, but not big enough.
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1)
    call add_intersection( &
         0.375_dp, 0.125_dp, num_intersections, intersections)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 2 .AND. &
         intersections(1, 2) == 0.375_dp .AND. &
         intersections(2, 2) == 0.125_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Intersections is allocated and has enough space for
    !         another intersection.
    deallocate(intersections)
    allocate(intersections(2, 2))
    num_intersections = 0
    call add_intersection( &
         0.75_dp, 1.0_dp, &
         num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 1 .AND. &
         intersections(1, 1) == 0.75_dp .AND. &
         intersections(2, 1) == 1.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Proposed intersection is a duplicate, matches first.
    intersections(1, 1) = 0.125_dp
    intersections(2, 1) = 0.75_dp
    intersections(1, 2) = 1.0_dp
    intersections(2, 2) = 0.25_dp
    num_intersections = 2

    call add_intersection( &
         0.125_dp, 0.75_dp, num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 2)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Proposed intersection is a duplicate, matches second.
    call add_intersection( &
         1.0_dp, 0.25_dp, num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 2)
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_intersection

  subroutine test_add_from_linearized(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer(c_int) :: status
    real(c_double) :: nodes1(2, 2), nodes2(3, 2)
    real(c_double) :: root_nodes1(3, 2), root_nodes2(3, 2)
    type(CurveData) :: first, second
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    real(c_double) :: linearization_error1, linearization_error2
    integer :: case_id
    character(19) :: name

    case_id = 1
    name = "add_from_linearized"
    num_intersections = 0

    ! CASE 1: Lines that intersect. Since lines, there are 0 subdivisions.
    first%start = 0.0_dp
    first%end_ = 1.0_dp
    nodes1(1, :) = 0
    nodes1(2, :) = 1
    first%nodes = nodes1

    second%start = 0.0_dp
    second%end_ = 1.0_dp
    nodes1(1, :) = [0.0_dp, 1.0_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    second%nodes = nodes1

    call add_from_linearized( &
         first, first%nodes, 0.0_dp, &
         second, second%nodes, 0.0_dp, &
         num_intersections, intersections, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         status == Status_SUCCESS .AND. &
         intersections(1, 1) == 0.5_dp .AND. &
         intersections(2, 1) == 0.5_dp)
    call print_status(name, case_id, case_success, success)
    num_intersections = 0

    ! CASE 2: Lines that **do not** intersect.
    first%start = 0.0_dp
    first%end_ = 1.0_dp
    nodes1(1, :) = 0
    nodes1(2, :) = 1
    first%nodes = nodes1

    second%start = 0.0_dp
    second%end_ = 1.0_dp
    nodes1(1, :) = [3.0_dp, 0.0_dp]
    nodes1(2, :) = [2.0_dp, 1.0_dp]
    second%nodes = nodes1

    call add_from_linearized( &
         first, first%nodes, 0.0_dp, &
         second, second%nodes, 0.0_dp, &
         num_intersections, intersections, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Quadratic curves that intersect after many (12) subdivisions.
    root_nodes1(1, :) = [0.25_dp, 0.4375_dp]
    root_nodes1(2, :) = [0.625_dp, 1.0_dp]
    root_nodes1(3, :) = [1.0_dp, 1.0_dp]
    first%start = 1365.0_dp / 4096.0_dp
    first%end_ = 1366.0_dp / 4096.0_dp
    ! NOTE: This is the result of
    !       call specialize_curve( &
    !            3, 2, root_nodes1, first%start, first%end_, ...)
    nodes2(1, :) = [134201344.0_dp, 201310207.0_dp]
    nodes2(2, :) = [134225920.0_dp, 201334786.0_dp]
    nodes2(3, :) = [134250496.0_dp, 201359356.0_dp]
    first%nodes = 0.5_dp**28 * nodes2

    root_nodes2(1, :) = [0.0_dp, 1.0_dp]
    root_nodes2(2, :) = [0.375_dp, 1.0_dp]
    root_nodes2(3, :) = [0.75_dp, 0.4375_dp]
    second%start = 2730.0_dp / 4096.0_dp
    second%end_ = 2731.0_dp / 4096.0_dp
    ! NOTE: This is the result of
    !       call specialize_curve( &
    !            3, 2, root_nodes2, second%start, second%end_, ...)
    nodes2(1, :) = [134184960.0_dp, 201359356.0_dp]
    nodes2(2, :) = [134209536.0_dp, 201334786.0_dp]
    nodes2(3, :) = [134234112.0_dp, 201310207.0_dp]
    second%nodes = 0.5_dp**28 * nodes2

    call linearization_error( &
         3, 2, first%nodes, linearization_error1)
    call linearization_error( &
         3, 2, second%nodes, linearization_error2)
    call add_from_linearized( &
         first, root_nodes1, linearization_error1, &
         second, root_nodes2, linearization_error2, &
         num_intersections, intersections, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         status == Status_SUCCESS .AND. &
         3.0_dp * intersections(1, 1) == 1.0_dp .AND. &
         3.0_dp * intersections(2, 1) == 2.0_dp)
    call print_status(name, case_id, case_success, success)
    num_intersections = 0

    ! CASE 4: Parallel lines that **do** intersect.
    first%start = 0.0_dp
    first%end_ = 1.0_dp
    nodes1(1, :) = 0
    nodes1(2, :) = 1
    first%nodes = nodes1

    second%start = 0.0_dp
    second%end_ = 1.0_dp
    nodes1(1, :) = [0.5_dp, 0.5_dp]
    nodes1(2, :) = [3.0_dp, 3.0_dp]
    second%nodes = nodes1

    call add_from_linearized( &
         first, first%nodes, 0.0_dp, &
         second, second%nodes, 0.0_dp, &
         num_intersections, intersections, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         status == Status_PARALLEL)
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_from_linearized

  subroutine test_endpoint_check(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    type(CurveData) :: first, second
    real(c_double) :: node_first(1, 2), node_second(1, 2)
    real(c_double) :: s, t
    real(c_double) :: nodes1(2, 2), nodes2(3, 2)
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    integer :: case_id
    character(14) :: name

    case_id = 1
    name = "endpoint_check"
    num_intersections = 0

    ! CASE 1: Endpoints that are not close to one another.
    node_first = 0
    s = 0.0_dp

    node_second = 1
    t = 0.0_dp

    ! NOTE: We use ``first`` and ``second`` without any actual info,
    !       but it won't matter because ``node_first`` and ``node_second``
    !       are not close.
    call endpoint_check( &
         first, node_first, 0.0_dp, &
         second, node_second, 0.0_dp, &
         num_intersections, intersections)
    case_success = ( &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: The endpoints are actually the same.
    nodes1(1, :) = 0
    nodes1(2, :) = 1
    first%nodes = nodes1
    s = 1.0_dp
    node_first(1, :) = first%nodes(2, :)

    nodes1(1, :) = 1
    nodes1(2, :) = [2.0_dp, 1.0_dp]
    second%nodes = nodes1
    t = 0.0_dp
    node_second(1, :) = second%nodes(1, :)

    call endpoint_check( &
         first, node_first, s, &
         second, node_second, t, &
         num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         intersections(1, 1) == s .AND. &
         intersections(2, 1) == t)
    call print_status(name, case_id, case_success, success)
    num_intersections = 0

    ! CASE 3: An intersection after one subdivision.
    first%start = 0.0_dp
    first%end_ = 0.5_dp
    nodes2(1, :) = 0
    nodes2(2, :) = [0.25_dp, 0.5_dp]
    nodes2(3, :) = 0.5_dp
    first%nodes = nodes2
    node_first(1, :) = first%nodes(3, :)

    second%start = 0.5_dp
    second%end_ = 1.0_dp
    nodes2(1, :) = 0.5_dp
    nodes2(2, :) = [0.5_dp, 0.0_dp]
    nodes2(3, :) = [1.0_dp, -0.5_dp]
    second%nodes = nodes2
    node_second(1, :) = second%nodes(1, :)

    call endpoint_check( &
         first, node_first, s, &
         second, node_second, t, &
         num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         intersections(1, 1) == 0.5_dp .AND. &
         intersections(2, 1) == 0.5_dp)
    call print_status(name, case_id, case_success, success)

  end subroutine test_endpoint_check

  subroutine test_tangent_bbox_intersection(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    type(CurveData) :: first, second
    real(c_double) :: nodes1(3, 2), nodes2(2, 2)
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    integer :: case_id
    character(25) :: name

    case_id = 1
    name = "tangent_bbox_intersection"
    num_intersections = 0

    ! CASE 1: Standard case of two quadratics that touch at a single endpoint.
    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 2.0_dp]
    nodes1(3, :) = [2.0_dp, 0.0_dp]
    first%nodes = nodes1

    nodes1(1, :) = [2.0_dp, 0.0_dp]
    nodes1(2, :) = [3.0_dp, 2.0_dp]
    nodes1(3, :) = [4.0_dp, 0.0_dp]
    second%nodes = nodes1

    call tangent_bbox_intersection( &
         first, second, num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         intersections(1, 1) == 1.0_dp .AND. &
         intersections(2, 1) == 0.0_dp)
    call print_status(name, case_id, case_success, success)
    num_intersections = 0

    ! CASE 2: Two quadratics that touch at both endpoints.
    nodes1(1, :) = 0
    nodes1(2, :) = [-1.0_dp, 0.5_dp]
    nodes1(3, :) = [0.0_dp, 1.0_dp]
    first%nodes = nodes1

    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.5_dp]
    nodes1(3, :) = [0.0_dp, 1.0_dp]
    second%nodes = nodes1

    call tangent_bbox_intersection( &
         first, second, num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 2 .AND. &
         intersections(1, 1) == 0.0_dp .AND. &
         intersections(2, 1) == 0.0_dp .AND. &
         intersections(2, 2) == 1.0_dp .AND. &
         intersections(2, 2) == 1.0_dp)
    call print_status(name, case_id, case_success, success)
    num_intersections = 0

    ! CASE 3: Two lines that don't touch at endpoints, but have
    !         tangent bounding boxes.
    nodes2(1, :) = 0
    nodes2(2, :) = [2.0_dp, 1.0_dp]
    first%nodes = nodes2

    nodes2(1, :) = [0.5_dp, 1.0_dp]
    nodes2(2, :) = [2.5_dp, 2.0_dp]
    second%nodes = nodes2

    call tangent_bbox_intersection( &
         first, second, num_intersections, intersections)
    case_success = (num_intersections == 0)
    call print_status(name, case_id, case_success, success)

  end subroutine test_tangent_bbox_intersection

  subroutine test_add_candidates(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    type(CurveData), allocatable :: candidates(:, :)
    type(CurveData) :: curve1, curve2, curve3
    type(CurveData) :: right1, left3, right3
    integer(c_int) :: num_candidates
    integer :: case_id
    character(14) :: name

    case_id = 1
    name = "add_candidates"

    ! Pre-populate the curves.
    allocate(curve1%nodes(3, 2))
    curve1%nodes(1, :) = [0.0_dp, 1.0_dp]
    curve1%nodes(2, :) = [1.0_dp, 2.0_dp]
    curve1%nodes(3, :) = [2.0_dp, 1.0_dp]
    allocate(curve2%nodes(2, 2))
    curve2%nodes(1, :) = [0.0_dp, 1.25_dp]
    curve2%nodes(2, :) = [2.0_dp, 1.25_dp]
    call subdivide_curve(curve1, curve3, right1)
    call subdivide_curve(curve3, left3, right3)

    ! CASE 1: Add when ``candidates`` is un-allocated. (Also ``enum_`` is
    !         ``Subdivide_FIRST``.)
    num_candidates = 0
    case_success = .NOT. allocated(candidates)
    call add_candidates( &
         candidates, num_candidates, curve1, curve2, Subdivide_FIRST)
    case_success = ( &
         case_success .AND. &
         num_candidates == 2 .AND. &
         all(shape(candidates) == [2, 2]) .AND. &
         curves_equal(candidates(1, 1), curve3) .AND. &
         curves_equal(candidates(2, 1), curve2) .AND. &
         curves_equal(candidates(1, 2), right1) .AND. &
         curves_equal(candidates(2, 2), curve2))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Add when ``candidates`` is allocated but too small. (Also
    !         ``enum_`` is ``Subdivide_BOTH``.)
    case_success = ( &
         num_candidates == 2 .AND. &
         allocated(candidates) .AND. &
         all(shape(candidates) == [2, 2]))
    call add_candidates( &
         candidates, num_candidates, curve3, curve1, Subdivide_BOTH)
    case_success = ( &
         case_success .AND. &
         num_candidates == 6 .AND. &
         all(shape(candidates) == [2, 6]) .AND. &
         curves_equal(candidates(1, 3), left3) .AND. &
         curves_equal(candidates(2, 3), curve3) .AND. &
         curves_equal(candidates(1, 4), left3) .AND. &
         curves_equal(candidates(2, 4), right1) .AND. &
         curves_equal(candidates(1, 5), right3) .AND. &
         curves_equal(candidates(2, 5), curve3) .AND. &
         curves_equal(candidates(1, 6), right3) .AND. &
         curves_equal(candidates(2, 6), right1))
    call print_status(name, case_id, case_success, success)
    deallocate(candidates)

    ! CASE 3: Add when ``candidates`` is bigger than necessary. (Also
    !         ``enum_`` is ``Subdivide_SECOND``.)
    num_candidates = 2
    allocate(candidates(2, 5))

    case_success = ( &
         num_candidates == 2 .AND. &
         allocated(candidates) .AND. &
         all(shape(candidates) == [2, 5]))
    call add_candidates( &
         candidates, num_candidates, curve2, curve3, Subdivide_SECOND)
    case_success = ( &
         case_success .AND. &
         num_candidates == 4 .AND. &
         all(shape(candidates) == [2, 5]) .AND. &
         curves_equal(candidates(1, 3), curve2) .AND. &
         curves_equal(candidates(2, 3), left3) .AND. &
         curves_equal(candidates(1, 4), curve2) .AND. &
         curves_equal(candidates(2, 4), right3))
    call print_status(name, case_id, case_success, success)
    deallocate(candidates)

    ! CASE 4: ``enum_`` does not require subdivision, ``candidates`` is
    !         **not** allocated.
    num_candidates = -6  ! This is a nonsense value.
    case_success = .NOT. allocated(candidates)
    call add_candidates( &
         candidates, num_candidates, curve1, curve2, Subdivide_NEITHER)
    case_success = ( &
         case_success .AND. &
         num_candidates == -6 .AND. &
         .NOT. allocated(candidates))
    call print_status(name, case_id, case_success, success)

    ! CASE 5: ``enum_`` does not require subdivision, ``candidates`` is
    !         allocated.
    num_candidates = 13  ! This does not reflect the size.
    allocate(candidates(2, 1))

    call add_candidates( &
         candidates, num_candidates, curve1, curve2, Subdivide_NEITHER)
    case_success = ( &
         case_success .AND. &
         num_candidates == 13 .AND. &
         all(shape(candidates) == [2, 1]))
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_candidates

  subroutine test_intersect_one_round(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: fixed_quadratic1(3, 2), fixed_quadratic2(3, 2)
    real(c_double) :: fixed_line1(2, 2), fixed_line2(2, 2)
    integer(c_int) :: num_candidates
    type(CurveData), allocatable :: candidates(:, :)
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    type(CurveData), allocatable :: next_candidates(:, :)
    integer(c_int) :: num_next_candidates
    integer(c_int) :: status
    type(CurveData) :: left1, right1, left2, right2
    integer :: case_id
    character(19) :: name

    case_id = 1
    name = "intersect_one_round"
    num_intersections = 0

    ! NOTE: ``fixed_quadratic1`` is a specialization of
    !       [0, 0], [1/2, 1], [1, 1] onto the interval [1/4, 1].
    fixed_quadratic1(1, :) = [0.25_dp, 0.4375_dp]
    fixed_quadratic1(2, :) = [0.625_dp, 1.0_dp]
    fixed_quadratic1(3, :) = [1.0_dp, 1.0_dp]
    ! NOTE: ``fixed_quadratic2`` is a specialization of
    !       [0, 1], [1/2, 1], [1, 0] onto the interval [0, 3/4].
    fixed_quadratic2(1, :) = [0.0_dp, 1.0_dp]
    fixed_quadratic2(2, :) = [0.375_dp, 1.0_dp]
    fixed_quadratic2(3, :) = [0.75_dp, 0.4375_dp]

    fixed_line1(1, :) = [0.0_dp, 0.0_dp]
    fixed_line1(2, :) = [1.0_dp, 1.0_dp]

    fixed_line2(1, :) = [0.0_dp, 1.0_dp]
    fixed_line2(2, :) = [1.0_dp, 0.0_dp]

    ! CASE 1: Simple test, non-linearized quadratics.
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve.
    candidates(1, 1)%nodes = fixed_quadratic1
    call subdivide_curve(candidates(1, 1), left1, right1)
    ! Populate the "second" curve.
    candidates(2, 1)%nodes = fixed_quadratic2
    call subdivide_curve(candidates(2, 1), left2, right2)

    call intersect_one_round( &
         fixed_quadratic1, fixed_quadratic2, num_candidates, candidates, &
         num_intersections, intersections, &
         next_candidates, num_next_candidates, status)
    case_success = ( &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0 .AND. &
         num_next_candidates == 4 .AND. &
         curves_equal(next_candidates(1, 1), left1) .AND. &
         curves_equal(next_candidates(2, 1), left2) .AND. &
         curves_equal(next_candidates(1, 2), left1) .AND. &
         curves_equal(next_candidates(2, 2), right2) .AND. &
         curves_equal(next_candidates(1, 3), right1) .AND. &
         curves_equal(next_candidates(2, 3), left2) .AND. &
         curves_equal(next_candidates(1, 4), right1) .AND. &
         curves_equal(next_candidates(2, 4), right2) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)
    deallocate(candidates)
    deallocate(next_candidates)

    ! CASE 2: First curve is linearized, second is not.
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve with a line.
    candidates(1, 1)%nodes = fixed_line1
    ! Populate the "second" curve.
    candidates(2, 1)%nodes = fixed_quadratic2
    call subdivide_curve(candidates(2, 1), left2, right2)

    call intersect_one_round( &
         fixed_line1, fixed_quadratic2, num_candidates, candidates, &
         num_intersections, intersections, &
         next_candidates, num_next_candidates, status)
    case_success = ( &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0 .AND. &
         num_next_candidates == 2 .AND. &
         curves_equal(next_candidates(1, 1), candidates(1, 1)) .AND. &
         curves_equal(next_candidates(2, 1), left2) .AND. &
         curves_equal(next_candidates(1, 2), candidates(1, 1)) .AND. &
         curves_equal(next_candidates(2, 2), right2) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)
    deallocate(candidates)
    deallocate(next_candidates)

    ! CASE 3: Second curve is linearized, first is not.
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve.
    candidates(1, 1)%nodes = fixed_quadratic1
    call subdivide_curve(candidates(1, 1), left1, right1)
    ! Populate the "second" curve with a line.
    candidates(2, 1)%nodes = fixed_line2

    call intersect_one_round( &
         fixed_quadratic1, fixed_line2, num_candidates, candidates, &
         num_intersections, intersections, &
         next_candidates, num_next_candidates, status)
    case_success = ( &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0 .AND. &
         num_next_candidates == 2 .AND. &
         curves_equal(next_candidates(1, 1), left1) .AND. &
         curves_equal(next_candidates(2, 1), candidates(2, 1)) .AND. &
         curves_equal(next_candidates(1, 2), right1) .AND. &
         curves_equal(next_candidates(2, 2), candidates(2, 1)) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)
    deallocate(candidates)
    deallocate(next_candidates)

    ! CASE 4: Both curves are linearized.
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve with a line.
    candidates(1, 1)%nodes = fixed_line1
    ! Populate the "second" curve with a line.
    candidates(2, 1)%nodes = fixed_line2

    call intersect_one_round( &
         fixed_line1, fixed_line2, num_candidates, candidates, &
         num_intersections, intersections, &
         next_candidates, num_next_candidates, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         intersections(1, 1) == 0.5_dp .AND. &
         intersections(2, 1) == 0.5_dp .AND. &
         num_next_candidates == 0 .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)
    deallocate(candidates)
    num_intersections = 0

    ! CASE 5: Failure caused by parallel lines that **do** intersect.
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve with a line.
    candidates(1, 1)%nodes = fixed_line1
    ! Populate the "second" curve with a line.
    allocate(candidates(2, 1)%nodes(2, 2))
    candidates(2, 1)%nodes(1, :) = [0.5_dp, 0.5_dp]
    candidates(2, 1)%nodes(2, :) = [3.0_dp, 3.0_dp]

    call intersect_one_round( &
         fixed_line1, candidates(2, 1)%nodes, num_candidates, candidates, &
         num_intersections, intersections, &
         next_candidates, num_next_candidates, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         num_next_candidates == 0 .AND. &
         status == Status_PARALLEL)
    call print_status(name, case_id, case_success, success)
    deallocate(candidates)

    ! CASE 6: Disjoint bounding boxes.
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve with a line.
    candidates(1, 1)%nodes = fixed_quadratic1
    ! Populate the "second" curve with a line.
    allocate(candidates(2, 1)%nodes(2, 2))
    candidates(2, 1)%nodes(1, :) = [1.0_dp, 1.25_dp]
    candidates(2, 1)%nodes(2, :) = [0.0_dp, 2.0_dp]

    call intersect_one_round( &
         fixed_quadratic1, candidates(2, 1)%nodes, &
         num_candidates, candidates, &
         num_intersections, intersections, &
         next_candidates, num_next_candidates, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         num_next_candidates == 0 .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)
    deallocate(candidates)

    ! CASE 7: Tangent bounding boxes (**with** an intersection), noting
    !         that tangency is only allowed for a pair with both curves
    !         can't be linearized (since tangency can be resolved in the
    !         linearized case by checking endpoints).
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve.
    allocate(candidates(1, 1)%nodes(3, 2))
    candidates(1, 1)%nodes(1, :) = [0.0_dp, 0.0_dp]
    candidates(1, 1)%nodes(2, :) = [0.5_dp, 1.0_dp]
    candidates(1, 1)%nodes(3, :) = [1.0_dp, 0.0_dp]
    ! [0 <= x <= 1.0], [0.0 <= y <= 1.0]
    ! Populate the "second" curve.
    allocate(candidates(2, 1)%nodes(3, 2))
    candidates(2, 1)%nodes(1, :) = [1.0_dp, 0.0_dp]
    candidates(2, 1)%nodes(2, :) = [1.5_dp, 0.5_dp]
    candidates(2, 1)%nodes(3, :) = [2.0_dp, -0.25_dp]
    ! [1.0 <= x <= 2.0], [-0.25 <= y <= 0.5]

    call intersect_one_round( &
         candidates(1, 1)%nodes, candidates(2, 1)%nodes, &
         num_candidates, candidates, &
         num_intersections, intersections, &
         next_candidates, num_next_candidates, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         intersections(1, 1) == 1.0_dp .AND. &
         intersections(2, 1) == 0.0_dp .AND. &
         num_next_candidates == 0 .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_intersect_one_round

  subroutine test_all_intersections(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double), allocatable :: intersections(:, :)
    real(c_double) :: linear1(2, 2), linear2(2, 2)
    real(c_double) :: quadratic1(3, 2), quadratic2(3, 2)
    real(c_double) :: cubic1(4, 2), cubic2(4, 2)
    integer(c_int) :: num_intersections
    integer(c_int) :: status
    integer :: case_id
    character(17) :: name

    case_id = 1
    name = "all_intersections"

    ! CASE 1: No intersections.
    linear1(1, :) = 0
    linear1(2, :) = 1
    linear2(1, :) = [3.0_dp, 3.0_dp]
    linear2(2, :) = [4.0_dp, 3.0_dp]

    call all_intersections( &
         2, linear1, 2, linear2, intersections, &
         num_intersections, status)
    case_success = ( &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0 .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Intersection of two quadratic curves.
    ! NOTE: The first candidate is a specialization of [0, 0], [1/2, 1], [1, 1]
    !       onto the interval [1/4, 1] and the second candidate is
    !       specialization of [0, 1], [1/2, 1], [1, 0] onto the interval
    !       [0, 3/4]. We expect them to intersect at s = 1/3, t = 2/3, which is
    !       the point [1/2, 3/4].
    quadratic1(1, :) = [0.25_dp, 0.4375_dp]
    quadratic1(2, :) = [0.625_dp, 1.0_dp]
    quadratic1(3, :) = 1
    quadratic2(1, :) = [0.0_dp, 1.0_dp]
    quadratic2(2, :) = [0.375_dp, 1.0_dp]
    quadratic2(3, :) = [0.75_dp, 0.4375_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         3 * intersections(1, 1) == 1.0_dp .AND. &
         3 * intersections(2, 1) == 2.0_dp .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Tangent curves, with a ``status`` failure due to parallel lines.
    quadratic1(1, :) = 0
    quadratic1(2, :) = [0.375_dp, 0.75_dp]
    quadratic1(3, :) = [0.75_dp, 0.375_dp]
    quadratic2(1, :) = [0.25_dp, 0.625_dp]
    quadratic2(2, :) = [0.625_dp, 0.25_dp]
    quadratic2(3, :) = 1
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         status == Status_PARALLEL)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Tangent curves, which cause the number of candidate pairs
    !         to become to high.
    quadratic1(1, :) = 0
    quadratic1(2, :) = [-0.5_dp, 1.5_dp]
    quadratic1(3, :) = 1
    quadratic2(1, :) = [-1.0_dp, 1.0_dp]
    quadratic2(2, :) = 0.5_dp
    quadratic2(3, :) = [0.0_dp, 2.0_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, status)
    case_success = ( &
         num_intersections == 1 .AND. &
         status == 88)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Badly scaled curves, which cause the subdivision proceass to
    !         take too many iterations before being "almost linear".
    ! NOTE: This is just the quadratic given by [0, 0], [4.5, 9], [9, 0] and
    !       the line given by [0, 8], [6, 0], but scaled up by a factor of
    !       2**14 == 16384. This makes the linearization error large.
    !       This definitely points out a flaw in using absolute vs. relative
    !       linearization error. However, there is an upside in not using
    !       relative error: it is faster to compute without normalization.
    quadratic1(1, :) = 0
    quadratic1(2, :) = [73728.0_dp, 147456.0_dp]
    quadratic1(3, :) = [147456.0_dp, 0.0_dp]
    linear1(1, :) = [0.0_dp, 131072.0_dp]
    linear1(2, :) = [98304.0_dp, 0.0_dp]
    call all_intersections( &
         3, quadratic1, 2, linear1, intersections, &
         num_intersections, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         status == Status_NO_CONVERGE)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Curves where there are duplicate intersections caused by
    !         bounding boxes touching at corners.
    quadratic1(1, :) = 0
    quadratic1(2, :) = [0.5_dp, 1.0_dp]
    quadratic1(3, :) = [1.0_dp, 0.0_dp]
    quadratic2(1, :) = [0.0_dp, 0.75_dp]
    quadratic2(2, :) = [0.5_dp, -0.25_dp]
    quadratic2(3, :) = [1.0_dp, 0.75_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, status)
    case_success = ( &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 2 .AND. &
         intersections(1, 1) == 0.25_dp .AND. &
         intersections(2, 1) == 0.25_dp .AND. &
         intersections(1, 2) == 0.75_dp .AND. &
         intersections(2, 2) == 0.75_dp .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Curves that **almost** touch at their endpoints, but don't
    !         actually cross. This causes a "wiggle fail" in
    !         ``from_linearized()``.
    cubic1(1, :) = [-0.7838204403623438_dp, -0.25519640597397464_dp]
    cubic1(2, :) = [-0.7894577677825452_dp, -0.24259531488131633_dp]
    cubic1(3, :) = [-0.7946421067207265_dp, -0.22976394420044136_dp]
    cubic1(4, :) = [-0.799367666650849_dp, -0.21671303774854855_dp]
    cubic2(1, :) = [-0.7993236103108717_dp, -0.21683567278362156_dp]
    cubic2(2, :) = [-0.8072986524226636_dp, -0.21898490744674426_dp]
    cubic2(3, :) = [-0.8152736945344552_dp, -0.2211341421098668_dp]
    cubic2(4, :) = [-0.8232487366462472_dp, -0.2232833767729893_dp]

    call all_intersections( &
         4, cubic1, 4, cubic2, intersections, &
         num_intersections, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         status == Status_WIGGLE_FAIL)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Same as CASE 7, but swap the inputs.
    call all_intersections( &
         4, cubic2, 4, cubic1, intersections, &
         num_intersections, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         status == Status_WIGGLE_FAIL)
    call print_status(name, case_id, case_success, success)

  end subroutine test_all_intersections

  subroutine test_all_intersections_abi(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: linear1(2, 2), linear2(2, 2)
    real(c_double) :: cubic1(4, 2)
    integer(c_int) :: num_intersections
    real(c_double) :: intersections1(2, 2), intersections2(2, 3)
    integer(c_int) :: status
    integer :: case_id
    character(21) :: name

    case_id = 1
    name = "all_intersections_abi"

    ! CASE 1: **Other** failure (overlapping lines).
    linear1(1, :) = [0.0_dp, 0.0_dp]
    linear1(2, :) = [2.0_dp, 4.0_dp]
    linear2(1, :) = [1.0_dp, 2.0_dp]
    linear2(2, :) = [3.0_dp, 6.0_dp]

    call all_intersections_abi( &
         2, linear1, 2, linear2, 2, intersections1, &
         num_intersections, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         status == Status_PARALLEL)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``intersections`` is not large enough.
    linear1(1, :) = [-3.0_dp, 0.0_dp]
    linear1(2, :) = [5.0_dp, 0.0_dp]
    cubic1(1, :) = [-7.0_dp, -9.0_dp]
    cubic1(2, :) = [9.0_dp, 13.0_dp]
    cubic1(3, :) = [-7.0_dp, -13.0_dp]
    cubic1(4, :) = [9.0_dp, 9.0_dp]

    call all_intersections_abi( &
         2, linear1, 4, cubic1, 2, intersections1, &
         num_intersections, status)
    case_success = ( &
         num_intersections == 3 .AND. &
         status == Status_INSUFFICIENT_SPACE)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Just case 7, but with large enough ``intersections``.
    call all_intersections_abi( &
         2, linear1, 4, cubic1, 3, intersections2, &
         num_intersections, status)
    case_success = ( &
         num_intersections == 3 .AND. &
         intersections2(1, 1) == 0.5_dp .AND. &
         intersections2(2, 1) == 0.5_dp .AND. &
         intersections2(1, 2) == 0.375_dp .AND. &
         intersections2(2, 2) == 0.25_dp .AND. &
         intersections2(1, 3) == 0.625_dp .AND. &
         intersections2(2, 3) == 0.75_dp .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_all_intersections_abi

  subroutine test_set_max_candidates(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer(c_int) :: orig_num_candidates, num_candidates1, num_candidates2
    integer :: case_id
    character(18) :: name

    ! NOTE: This is **also** a test for ``get_max_candidates``.

    case_id = 1
    name = "set_max_candidates"

    ! CASE 1: Check that we can properly set and get the value.
    call get_max_candidates(orig_num_candidates)

    num_candidates1 = 33
    call set_max_candidates(num_candidates1)
    call get_max_candidates(num_candidates2)
    case_success = ( &
         orig_num_candidates == 64 .AND. &
         num_candidates1 == num_candidates2)
    call print_status(name, case_id, case_success, success)
    ! Restore the original value.
    call set_max_candidates(orig_num_candidates)

  end subroutine test_set_max_candidates

  subroutine test_set_similar_ulps(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer(c_int) :: orig_num_bits, num_bits1, num_bits2
    integer :: case_id
    character(16) :: name

    ! NOTE: This is **also** a test for ``get_similar_ulps``.

    case_id = 1
    name = "set_similar_ulps"

    ! CASE 1: Check that we can properly set and get the value.
    call get_similar_ulps(orig_num_bits)

    num_bits1 = 7
    call set_similar_ulps(num_bits1)
    call get_similar_ulps(num_bits2)
    case_success = ( &
         orig_num_bits == 1 .AND. &
         num_bits1 == num_bits2)
    call print_status(name, case_id, case_success, success)
    ! Restore the original value.
    call set_similar_ulps(orig_num_bits)

  end subroutine test_set_similar_ulps

end module test_curve_intersection
