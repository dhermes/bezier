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
       Status_SUCCESS, Status_BAD_MULTIPLICITY, Status_NO_CONVERGE, &
       Status_INSUFFICIENT_SPACE, Status_SINGULAR
  use curve, only: &
       CurveData, evaluate_multi, specialize_curve, subdivide_nodes, &
       curves_equal, subdivide_curve
  use curve_intersection, only: &
       BoxIntersectionType_INTERSECTION, BoxIntersectionType_TANGENT, &
       BoxIntersectionType_DISJOINT, &
       Subdivide_FIRST, Subdivide_SECOND, &
       Subdivide_BOTH, Subdivide_NEITHER, &
       linearization_error, segment_intersection, &
       newton_refine_intersect, bbox_intersect, parallel_lines_parameters, &
       line_line_collide, convex_hull_collide, newton_simple_root, &
       newton_double_root, newton_iterate, full_newton_nonzero, full_newton, &
       from_linearized, &
       bbox_line_intersect, check_lines, add_intersection, &
       add_from_linearized, endpoint_check, tangent_bbox_intersection, &
       add_candidates, intersect_one_round, make_same_degree, &
       add_coincident_parameters, all_intersections, all_intersections_abi, &
       set_max_candidates, get_max_candidates
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private &
       test_linearization_error, test_segment_intersection, &
       test_newton_refine_intersect, test_bbox_intersect, &
       test_parallel_lines_parameters, test_line_line_collide, &
       test_convex_hull_collide, test_newton_simple_root, &
       test_newton_double_root, check_closer, test_newton_iterate, &
       test_full_newton_nonzero, test_full_newton, test_from_linearized, &
       test_bbox_line_intersect, test_check_lines, test_add_intersection, &
       test_add_from_linearized, test_endpoint_check, &
       test_tangent_bbox_intersection, test_add_candidates, &
       test_intersect_one_round, test_make_same_degree, &
       test_add_coincident_parameters, test_all_intersections, &
       test_all_intersections_abi, test_set_max_candidates
  public curve_intersection_all_tests

contains

  subroutine curve_intersection_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_linearization_error(success)
    call test_segment_intersection(success)
    call test_newton_refine_intersect(success)
    call test_bbox_intersect(success)
    call test_parallel_lines_parameters(success)
    call test_line_line_collide(success)
    call test_convex_hull_collide(success)
    call test_newton_simple_root(success)
    call test_newton_double_root(success)
    call test_newton_iterate(success)
    call test_full_newton_nonzero(success)
    call test_full_newton(success)
    call test_from_linearized(success)
    call test_bbox_line_intersect(success)
    call test_check_lines(success)
    call test_add_intersection(success)
    call test_add_from_linearized(success)
    call test_endpoint_check(success)
    call test_tangent_bbox_intersection(success)
    call test_add_candidates(success)
    call test_intersect_one_round(success)
    call test_make_same_degree(success)
    call test_add_coincident_parameters(success)
    call test_all_intersections(success)
    call test_all_intersections_abi(success)
    call test_set_max_candidates(success)

  end subroutine curve_intersection_all_tests

  subroutine test_linearization_error(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 2), nodes2(2, 3), nodes3(2, 5)
    real(c_double) :: nodes4(3, 3), nodes5(2, 4), nodes6(2, 6)
    real(c_double) :: left_nodes2(2, 3), right_nodes2(2, 3)
    real(c_double) :: error1, error2, expected
    integer :: case_id
    character(19) :: name

    case_id = 1
    name = "linearization_error"

    ! CASE 1: Linear curve (i.e. no error).
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 2.0_dp]
    call linearization_error( &
         2, 2, nodes1, error1)
    case_success = (error1 == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Degree-elevated line (i.e. no error) as quadratic.
    nodes2(:, 1) = 0
    nodes2(:, 2) = [0.5_dp, 1.0_dp]
    nodes2(:, 3) = [1.0_dp, 2.0_dp]
    call linearization_error( &
         3, 2, nodes2, error1)
    case_success = (error1 == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Degree-elevated line (i.e. no error) as quartic.
    nodes3(:, 1) = 0
    nodes3(:, 2) = [0.25_dp, 0.5_dp]
    nodes3(:, 3) = [0.5_dp, 1.0_dp]
    nodes3(:, 4) = [0.75_dp, 1.5_dp]
    nodes3(:, 5) = [1.0_dp, 2.0_dp]
    call linearization_error( &
         5, 2, nodes3, error1)
    case_success = (error1 == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Line with bad parameterization.
    ! NOTE: This is the line 3 y = 4 x, but with the parameterization
    !       x(s) = 3 s (4 - 3 s).
    nodes2(:, 1) = 0
    nodes2(:, 2) = [6.0_dp, 8.0_dp]
    nodes2(:, 3) = [3.0_dp, 4.0_dp]
    call linearization_error( &
         3, 2, nodes2, error1)
    ! D^2 v = [-9, -12]
    expected = 0.125_dp * 2 * 1 * 15.0_dp
    case_success = (error1 == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Quadratic curve.
    nodes2(:, 1) = 0
    nodes2(:, 2) = 1
    nodes2(:, 3) = [5.0_dp, 6.0_dp]
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
    nodes4(:, 1) = [1.5_dp, 0.0_dp, 6.25_dp]
    nodes4(:, 2) = [3.5_dp, -5.0_dp, 10.25_dp]
    nodes4(:, 3) = [8.5_dp, 2.0_dp, 10.25_dp]
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
    nodes3(:, 1) = [1.0_dp, 1.25_dp]
    nodes3(:, 2) = [-0.5_dp, 0.5_dp]
    nodes3(:, 3) = [-0.5_dp, 2.0_dp]
    nodes3(:, 4) = [1.0_dp, -1.0_dp]
    nodes3(:, 5) = [4.0_dp, 5.0_dp]
    call linearization_error( &
         5, 2, nodes3, error1)
    ! D^2 v = [1.5, 2.25], [1.5, -4.5], [1.5, 9]
    expected = 0.125_dp * 4 * 3 * sqrt(1.5_dp**2 + 9.0_dp**2)
    case_success = (abs(error1 - expected) <= spacing(expected))
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Cubic curve.
    nodes5(:, 1) = 0
    nodes5(:, 2) = 1
    nodes5(:, 3) = [5.0_dp, 6.0_dp]
    nodes5(:, 4) = [6.0_dp, 7.0_dp]
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
    nodes3(:, 1) = 0
    nodes3(:, 2) = 1
    nodes3(:, 3) = [5.0_dp, 6.0_dp]
    nodes3(:, 4) = [6.0_dp, 7.0_dp]
    nodes3(:, 5) = [4.0_dp, 7.0_dp]
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
    nodes6(:, 1) = 0
    nodes6(:, 2) = 1
    nodes6(:, 3) = [7.0_dp, 3.0_dp]
    nodes6(:, 4) = [11.0_dp, 8.0_dp]
    nodes6(:, 5) = [15.0_dp, 1.0_dp]
    nodes6(:, 6) = [16.0_dp, -3.0_dp]
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
    real(c_double) :: start0(2), end0(2)
    real(c_double) :: start1(2), end1(2)
    real(c_double) :: s, t
    logical(c_bool) :: si_success
    integer :: case_id
    character(20) :: name

    case_id = 1
    name = "segment_intersection"

    ! CASE 1: Segments that intersect.
    start0 = [1.75_dp, 2.125_dp]
    end0 = [-1.25_dp, 1.625_dp]
    start1 = [-0.25_dp, 2.625_dp]
    end1 = [1.75_dp, 1.625_dp]
    ! D0 x D1 == 4.0, so there will be no round-off in answer.
    call segment_intersection( &
         start0, end0, start1, end1, s, t, si_success)
    case_success = (si_success .AND. s == 0.25_dp .AND. t == 0.625_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Parallel segments.
    start0 = [0.0_dp, 0.5_dp]
    end0 = [0.0_dp, -0.5_dp]
    start1 = [0.0_dp, 1.0_dp]
    end1 = [0.0_dp, -1.0_dp]
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
    real(c_double) :: nodes3(2, 3), nodes4(2, 3), nodes5(2, 5)
    real(c_double) :: known_s, known_t
    real(c_double) :: wrong_s, wrong_t
    real(c_double) :: new_s, new_t
    integer(c_int) :: status
    integer :: i
    integer :: case_id
    character(23) :: name

    case_id = 1
    name = "newton_refine_intersect"

    ! CASE 1: Intersection of two lines.
    ! Newton's method is exact on linear problems so will
    ! always converge after one step.
    nodes1(:, 1) = 0
    nodes1(:, 2) = 1
    known_s = 0.75_dp
    nodes2(:, 1) = [1.0_dp, 0.0_dp]
    nodes2(:, 2) = [0.0_dp, 3.0_dp]
    known_t = 0.25_dp
    ! NOTE: By construction, the Jacobian matrix will be
    !           [1, 1], [1, -3]
    !       which has determinant -4.0, hence there will
    !       be no round-off when solving.
    wrong_s = known_s - 0.125_dp
    wrong_t = known_t + 0.125_dp
    call newton_refine_intersect( &
         wrong_s, 2, nodes1, wrong_t, 2, nodes2, new_s, new_t, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         new_s == known_s .AND. &
         new_t == known_t .AND. &
         curves_intersect(nodes1, known_s, nodes2, known_t))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Mixed degree (line and quadratic).
    nodes3(:, 1) = 0
    nodes3(:, 2) = [0.5_dp, 1.0_dp]
    nodes3(:, 3) = [1.0_dp, 0.0_dp]
    known_s = 0.5_dp
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = [0.0_dp, 1.0_dp]
    known_t = 0.5_dp
    ! NOTE: By construction, the Jacobian matrix will be
    !           [1, 1], [1, -1]
    !       which has determinant -2.0, hence there will
    !       be no round-off when solving.
    wrong_s = 0.25_dp
    wrong_t = 0.25_dp
    call newton_refine_intersect( &
         wrong_s, 3, nodes3, wrong_t, 2, nodes1, new_s, new_t, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         new_s == 0.4375_dp .AND. &
         abs(known_s - new_s) < abs(known_s - wrong_s) .AND. &
         new_t == 0.5625_dp .AND. &
         abs(known_t - new_t) < abs(known_t - wrong_t) .AND. &
         curves_intersect(nodes3, known_s, nodes1, known_t))
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Early exit (i.e. already an intersection).
    nodes3(:, 1) = 0
    nodes3(:, 2) = [0.5_dp, 1.0_dp]
    nodes3(:, 3) = [1.0_dp, 0.0_dp]
    known_s = 0.25_dp
    nodes4(:, 1) = [1.0_dp, 0.75_dp]
    nodes4(:, 2) = [0.5_dp, -0.25_dp]
    nodes4(:, 3) = [0.0_dp, 0.75_dp]
    known_t = 0.75_dp
    call newton_refine_intersect( &
         known_s, 3, nodes3, known_t, 3, nodes4, new_s, new_t, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         new_s == known_s .AND. &
         new_t == known_t .AND. &
         curves_intersect(nodes3, known_s, nodes4, known_t))
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Intersect quadratics.
    nodes3(:, 1) = 0
    nodes3(:, 2) = [0.5_dp, 1.0_dp]
    nodes3(:, 3) = [1.0_dp, 0.0_dp]
    known_s = 0.25_dp
    nodes4(:, 1) = [1.0_dp, 0.75_dp]
    nodes4(:, 2) = [0.5_dp, -0.25_dp]
    nodes4(:, 3) = [0.0_dp, 0.75_dp]
    known_t = 0.75_dp
    ! NOTE: By construction, the Jacobian matrix will be
    !           [1, 3/4], [1, -5/4]
    !       which has determinant -2.0, hence there will
    !       be no round-off when solving.
    wrong_s = known_s + 0.0625_dp
    wrong_t = known_t + 0.0625_dp
    call newton_refine_intersect( &
         wrong_s, 3, nodes3, wrong_t, 3, nodes4, new_s, new_t, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         new_s == 0.2421875_dp .AND. &
         abs(known_s - new_s) < abs(known_s - wrong_s) .AND. &
         new_t == 0.7578125_dp .AND. &
         abs(known_t - new_t) < abs(known_t - wrong_t) .AND. &
         curves_intersect(nodes3, known_s, nodes4, known_t))
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Convergence test.
    nodes5(:, 1) = 0
    nodes5(:, 2) = [0.25_dp, 1.0_dp]
    nodes5(:, 3) = [0.5_dp, -0.75_dp]
    nodes5(:, 4) = [0.75_dp, 1.0_dp]
    nodes5(:, 5) = [1.0_dp, 0.0_dp]
    known_s = 0.5_dp
    nodes1(:, 1) = [0.5_dp, 0.0_dp]
    nodes1(:, 2) = [0.5_dp, 1.0_dp]
    known_t = 0.21875_dp
    ! We start with a very wrong guess and update it three times.
    wrong_s = 0.0_dp
    wrong_t = 0.0_dp
    case_success = .TRUE.
    do i = 1, 3
       call newton_refine_intersect( &
            wrong_s, 5, nodes5, wrong_t, 2, nodes1, new_s, new_t, status)
       wrong_s = new_s
       wrong_t = new_t
       if (i == 1) then
          case_success = ( &
               case_success .AND. &
               status == Status_SUCCESS .AND. &
               new_s == known_s .AND. &
               new_t == 2.0_dp)
       else
          case_success = ( &
               case_success .AND. &
               status == Status_SUCCESS .AND. &
               new_s == known_s .AND. &
               new_t == known_t)
       end if
    end do

    case_success = ( &
         case_success .AND. &
         curves_intersect(nodes5, known_s, nodes1, known_t))
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Singular Jacobian.
    nodes3(:, 1) = [0.5_dp, 0.0_dp]
    nodes3(:, 2) = 1
    nodes3(:, 3) = [1.5_dp, 0.0_dp]
    known_s = 0.5_dp
    nodes1(:, 1) = [0.0_dp, 0.5_dp]
    nodes1(:, 2) = [1.0_dp, 0.5_dp]
    known_t = 0.5_dp
    ! NOTE: By construction, the Jacobian matrix
    !           [1, -1], [2(1 - 2s), 0]
    !       will evaluate at ``1/2, 1/2`` to
    !           [1, -1], [0, 0]
    !       which has determinant 0.0 exactly.
    call newton_refine_intersect( &
         known_s, 3, nodes3, known_t, 2, nodes1, new_s, new_t, status)
    case_success = (status == Status_SINGULAR)
    call print_status(name, case_id, case_success, success)

  end subroutine test_newton_refine_intersect

  function curves_intersect(nodes1, s, nodes2, t) result(predicate)
    real(c_double), intent(in) :: nodes1(:, :), nodes2(:, :)
    real(c_double), intent(in) :: s, t
    logical(c_bool) :: predicate
    ! Variables outside of signature.
    real(c_double) :: point1(1, 2), point2(1, 2)

    call evaluate_multi( &
         size(nodes1, 2), size(nodes1, 1), nodes1, 1, [s], point1)
    call evaluate_multi( &
         size(nodes2, 2), size(nodes2, 1), nodes2, 1, [t], point2)
    predicate = all(point1 == point2)

  end function curves_intersect

  subroutine test_bbox_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: unit_square(2, 4), other(2, 4), delta
    integer(c_int) :: enum_, i
    integer :: case_id
    character(14) :: name

    case_id = 1
    name = "bbox_intersect"
    unit_square(:, 1) = 0
    unit_square(:, 2) = [1.0_dp, 0.0_dp]
    unit_square(:, 3) = 1
    unit_square(:, 4) = [0.0_dp, 1.0_dp]

    ! CASE 1: Intersecting bbox-es.
    forall (i = 1:4)
       other(:, i) = unit_square(:, i) + 0.5_dp
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Far apart bbox-es.
    forall (i = 1:4)
       other(:, i) = unit_square(:, i) + 100.0_dp
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_DISJOINT)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Disjoint bbox-es that have an "aligned" edge.
    forall (i = 1:4)
       other(:, i) = unit_square(:, i) + [1.0_dp, 2.0_dp]
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_DISJOINT)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Tangent bbox-es.
    forall (i = 1:4)
       other(:, i) = unit_square(:, i) + [1.0_dp, 0.0_dp]
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_TANGENT)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Almost tangent bbox-es.
    delta = 1.0_dp + spacing(1.0_dp)
    forall (i = 1:4)
       other(:, i) = unit_square(:, i) + [delta, 0.0_dp]
    end forall
    call bbox_intersect( &
         4, unit_square, 4, other, enum_)
    case_success = (enum_ == BoxIntersectionType_DISJOINT)
    call print_status(name, case_id, case_success, success)

  end subroutine test_bbox_intersect

  subroutine test_parallel_lines_parameters(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    logical(c_bool) :: disjoint
    real(c_double) :: start0(2), end0(2)
    real(c_double) :: start1(2), end1(2)
    real(c_double) :: parameters(2, 2)
    integer :: case_id
    character(25) :: name

    case_id = 1
    name = "parallel_lines_parameters"

    ! CASE 1: Same line, but segments don't overlap.
    start0 = 0
    end0 = [3.0_dp, 4.0_dp]
    start1 = [6.0_dp, 8.0_dp]
    end1 = [9.0_dp, 12.0_dp]
    call parallel_lines_parameters( &
         start0, end0, start1, end1, disjoint, parameters)
    case_success = (disjoint)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Same as CASE 1, but reverse the direction of line 1.
    call parallel_lines_parameters( &
         start0, end0, end1, start1, disjoint, parameters)
    case_success = (disjoint)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Same as CASE 1, but swap the lines.
    call parallel_lines_parameters( &
         start1, end1, start0, end0, disjoint, parameters)
    case_success = (disjoint)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Same as CASE 1, but reverse the direction of line 0.
    call parallel_lines_parameters( &
         end0, start0, start1, end1, disjoint, parameters)
    case_success = (disjoint)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Same line, ``start1`` contained in segment "0".
    start0 = [6.0_dp, -3.0_dp]
    end0 = [-7.0_dp, 1.0_dp]
    start1 = [1.125_dp, -1.5_dp]
    end1 = [-5.375_dp, 0.5_dp]
    call parallel_lines_parameters( &
         start0, end0, start1, end1, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [0.375_dp, 0.0_dp]) .AND. &
         all(parameters(:, 2) == [0.875_dp, 1.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Same as CASE 5, but reverse the direction of line 1.
    call parallel_lines_parameters( &
         start0, end0, end1, start1, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [0.875_dp, 0.0_dp]) .AND. &
         all(parameters(:, 2) == [0.375_dp, 1.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Same as CASE 5, but swap the lines.
    call parallel_lines_parameters( &
         start1, end1, start0, end0, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [0.0_dp, 0.375_dp]) .AND. &
         all(parameters(:, 2) == [1.0_dp, 0.875_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Same line, ``end1`` contained in segment "0".
    start0 = [1.0_dp, 2.0_dp]
    end0 = [3.0_dp, 5.0_dp]
    start1 = [-2.0_dp, -2.5_dp]
    end1 = [2.0_dp, 3.5_dp]
    call parallel_lines_parameters( &
         start0, end0, start1, end1, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [0.0_dp, 0.75_dp]) .AND. &
         all(parameters(:, 2) == [0.5_dp, 1.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Same as CASE 8, but reverse the direction of line 1.
    call parallel_lines_parameters( &
         start0, end0, end1, start1, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [0.5_dp, 0.0_dp]) .AND. &
         all(parameters(:, 2) == [0.0_dp, 0.25_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 10: Same as CASE 8, but swap the lines.
    call parallel_lines_parameters( &
         start1, end1, start0, end0, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [0.75_dp, 0.0_dp]) .AND. &
         all(parameters(:, 2) == [1.0_dp, 0.5_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 11: Same line, segment "0" fully contained in segment "1".
    start0 = [-9.0_dp, 0.0_dp]
    end0 = [4.0_dp, 5.0_dp]
    start1 = [23.5_dp, 12.5_dp]
    end1 = [-28.5_dp, -7.5_dp]
    call parallel_lines_parameters( &
         start0, end0, start1, end1, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [1.0_dp, 0.375_dp]) .AND. &
         all(parameters(:, 2) == [0.0_dp, 0.625_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 12: Same as CASE 11, but reverse the direction of line 1.
    call parallel_lines_parameters( &
         start0, end0, end1, start1, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [0.0_dp, 0.375_dp]) .AND. &
         all(parameters(:, 2) == [1.0_dp, 0.625_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 13: Same as CASE 11, but swap the lines.
    call parallel_lines_parameters( &
         start1, end1, start0, end0, disjoint, parameters)
    case_success = ( &
         .NOT. disjoint .AND. &
         all(parameters(:, 1) == [0.625_dp, 0.0_dp]) .AND. &
         all(parameters(:, 2) == [0.375_dp, 1.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 14: Parallel, but different lines.
    start0 = [3.0_dp, 2.0_dp]
    end0 = [3.0_dp, 0.75_dp]
    start1 = 0
    end1 = [0.0_dp, 2.0_dp]
    call parallel_lines_parameters( &
         start0, end0, start1, end1, disjoint, parameters)
    case_success = (disjoint)
    call print_status(name, case_id, case_success, success)

  end subroutine test_parallel_lines_parameters

  subroutine test_line_line_collide(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: line1(2, 2), line2(2, 2)
    logical(c_bool) :: collision
    integer :: case_id
    character(17) :: name

    case_id = 1
    name = "line_line_collide"

    ! Will hold throughout.
    line1(:, 1) = 0
    line1(:, 2) = 1

    ! CASE 1: Lines in general position (i.e. not parallel) that collide.
    line2(:, 1) = [0.0_dp, 1.0_dp]
    line2(:, 2) = [1.0_dp, 0.0_dp]
    call line_line_collide(line1, line2, collision)
    case_success = (collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Lines in general position (i.e. not parallel) that don't collide.
    line2(:, 1) = [3.0_dp, 4.0_dp]
    line2(:, 2) = [3.0_dp, 5.0_dp]
    call line_line_collide(line1, line2, collision)
    case_success = (.NOT. collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Segments that lie on parallel (but distinct) lines.
    line2(:, 1) = [0.0_dp, 1.0_dp]
    line2(:, 2) = [1.0_dp, 2.0_dp]
    call line_line_collide(line1, line2, collision)
    case_success = (.NOT. collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Segments that lie on the same line, but are disjoint.
    line2(:, 1) = 2
    line2(:, 2) = 3
    call line_line_collide(line1, line2, collision)
    case_success = (.NOT. collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Segments that lie on the same line and overlap.
    line2(:, 1) = 0.5_dp
    line2(:, 2) = 1.5_dp
    call line_line_collide(line1, line2, collision)
    case_success = (collision)
    call print_status(name, case_id, case_success, success)

  end subroutine test_line_line_collide

  subroutine test_convex_hull_collide(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 2), nodes2(2, 3)
    real(c_double), allocatable :: polygon1(:, :), polygon2(:, :)
    logical(c_bool) :: collision
    integer :: case_id
    character(19) :: name

    case_id = 1
    name = "convex_hull_collide"

    ! CASE 1: Convex hulls are both lines (Also, both ``polygon1`` and
    !         ``polygon2`` need to be allocated.)
    nodes2(:, 1) = 0
    nodes2(:, 2) = 1
    nodes2(:, 3) = 2
    nodes1(:, 1) = [0.0_dp, 1.0_dp]
    nodes1(:, 2) = [0.0_dp, 2.0_dp]
    call convex_hull_collide( &
         3, nodes2, polygon1, 2, nodes1, polygon2, collision)
    case_success = (.NOT. collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: At least one of the convex hulls is not a line (Also,
    !         ``polygon1`` does not need to be resized but ``polygon2``
    !         does.)
    nodes1(:, 1) = 0
    nodes1(:, 2) = 1
    nodes2(:, 1) = [0.0_dp, 0.5_dp]
    nodes2(:, 2) = 1
    nodes2(:, 3) = [2.0_dp, 0.5_dp]
    call convex_hull_collide( &
         2, nodes1, polygon1, 3, nodes2, polygon2, collision)
    case_success = (collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Same as CASE 2 with inputs swapped (Also, ``polygon2`` does not
    !         need to be resized but ``polygon1`` does.)
    deallocate(polygon1)
    allocate(polygon1(2, 2))
    call convex_hull_collide( &
         3, nodes2, polygon1, 2, nodes1, polygon2, collision)
    case_success = (collision)
    call print_status(name, case_id, case_success, success)

  end subroutine test_convex_hull_collide

  subroutine test_newton_simple_root(success)

    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(18) :: name
    real(c_double) :: nodes1(2, 3), nodes2(2, 2)
    real(c_double) :: first_deriv1(2, 2), first_deriv2(2, 1)
    real(c_double) :: jacobian(2, 2)
    real(c_double) :: func_val(2, 1)

    case_id = 1
    name = "newton_simple_root"

    ! CASE 1: F(s, t) != [0, 0].
    !           B1(s) = [s(s + 2) ]
    !                   [4s(1 - s)]
    !           B2(t) = [1 + 3t]
    !                   [4 - 4t]
    !           DF    = [2 + 2s, -3]
    !                   [4 - 8s,  4]
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 2.0_dp]
    nodes1(:, 3) = [3.0_dp, 0.0_dp]
    first_deriv1(:, 1) = [2.0_dp, 4.0_dp]
    first_deriv1(:, 2) = [4.0_dp, -4.0_dp]
    nodes2(:, 1) = [1.0_dp, 4.0_dp]
    nodes2(:, 2) = [4.0_dp, 0.0_dp]
    first_deriv2(:, 1) = [3.0_dp, -4.0_dp]
    call newton_simple_root( &
         0.5_dp, 3, nodes1, first_deriv1, &
         0.25_dp, 2, nodes2, first_deriv2, jacobian, func_val)
    case_success = ( &
         all(jacobian(:, 1) == [3.0_dp, 0.0_dp]) .AND. &
         all(jacobian(:, 2) == [-3.0_dp, 4.0_dp]) .AND. &
         all(func_val(:, 1) == [-0.5_dp, -2.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: F(s, t) == [0, 0].
    !           B1(s) = [2s(1 + s)]
    !                   [6s(1 - s)]
    !           B2(t) = [21t]
    !                   [ 9t]
    !           DF    = [2 +  4s, -21]
    !                   [6 - 12s,  -9]
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 3.0_dp]
    nodes1(:, 3) = [4.0_dp, 0.0_dp]
    first_deriv1(:, 1) = [2.0_dp, 6.0_dp]
    first_deriv1(:, 2) = [6.0_dp, -6.0_dp]
    nodes2(:, 1) = 0
    nodes2(:, 2) = [21.0_dp, 9.0_dp]
    first_deriv2(:, 1) = [21.0_dp, 9.0_dp]
    call newton_simple_root( &
         0.75_dp, 3, nodes1, first_deriv1, &
         0.125_dp, 2, nodes2, first_deriv2, jacobian, func_val)
    case_success = ( &
         all(func_val == 0.0_dp))
    call print_status(name, case_id, case_success, success)

  end subroutine test_newton_simple_root

  subroutine test_newton_double_root(success)

    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(18) :: name
    real(c_double) :: nodes1(2, 3), nodes2(2, 3)
    real(c_double) :: first_deriv1(2, 2), first_deriv2(2, 2)
    real(c_double) :: second_deriv1(2, 1), second_deriv2(2, 1)
    real(c_double) :: modified_lhs(2, 2)
    real(c_double) :: modified_rhs(2, 1)

    case_id = 1
    name = "newton_double_root"

    ! Same for both.
    !           B1(s) = [s(s + 2) ]
    !                   [4s(1 - s)]
    !           B2(t) = [1 + 3t]
    !                   [4 - 4t]
    !           DF    = [2 + 2s, -3]
    !                   [4 - 8s,  4]
    !           B1(s) = [4s(1 - s)]
    !                   [       2s]
    !           B2(t) = [2(2t^2 - 2t + 1)]
    !                   [              2t]
    !           B1'(s) x B2'(s) = -16(s + t - 1)
    !           DG    = [4 - 8s, 4 - 8t]
    !                   [     2,     -2]
    !                   [   -16,    -16]
    !           DG^T DG = [   4(16s^2 - 16s + 69), 4(16st - 8s - 8t + 67)]
    !                     [4(16st - 8s - 8t + 67),    4(16t^2 - 16t + 69)]
    !           DG^T G  =
    !     [4(8s^3 - 12s^2 + 8st^2 - 8st + 73s - 4t^2 + 67t - 66)]
    !     [4(8s^2t - 4s^2 - 8st + 67s + 8t^3 - 12t^2 + 73t - 66)]
    nodes1(:, 1) = 0
    nodes1(:, 2) = [2.0_dp, 1.0_dp]
    nodes1(:, 3) = [0.0_dp, 2.0_dp]
    first_deriv1(:, 1) = [4.0_dp, 2.0_dp]
    first_deriv1(:, 2) = [-4.0_dp, 2.0_dp]
    second_deriv1(:, 1) = [-8.0_dp, 0.0_dp]
    nodes2(:, 1) = [2.0_dp, 0.0_dp]
    nodes2(:, 2) = [0.0_dp, 1.0_dp]
    nodes2(:, 3) = [2.0_dp, 2.0_dp]
    first_deriv2(:, 1) = [-4.0_dp, 2.0_dp]
    first_deriv2(:, 2) = [4.0_dp, 2.0_dp]
    second_deriv2(:, 1) = [8.0_dp, 0.0_dp]

    ! CASE 1: G(s, t) != [0, 0, 0].
    call newton_double_root( &
         0.75_dp, 3, nodes1, first_deriv1, second_deriv1, &
         0.25_dp, 3, nodes2, first_deriv2, second_deriv2, &
         modified_lhs, modified_rhs)
    case_success = ( &
         all(modified_lhs(:, 1) == [264.0_dp, 248.0_dp]) .AND. &
         all(modified_lhs(:, 2) == [248.0_dp, 264.0_dp]) .AND. &
         all(modified_rhs(:, 1) == [3.0_dp, -3.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: G(s, t) == [0, 0, 0].
    call newton_double_root( &
         0.5_dp, 3, nodes1, first_deriv1, second_deriv1, &
         0.5_dp, 3, nodes2, first_deriv2, second_deriv2, &
         modified_lhs, modified_rhs)
    case_success = ( &
         all(modified_rhs(:, 1) == 0.0_dp))
    call print_status(name, case_id, case_success, success)

  end subroutine test_newton_double_root

  subroutine check_closer( &
       s, current_s, expected_s, t, current_t, expected_t, success)

    ! Checks if ``current_s`` is closer to ``expected_s`` than ``s`` was
    ! but also makes sure it is not as close as ``sqrt(machine precision)``
    ! in relative error.

    real(c_double), intent(in) :: s
    real(c_double), intent(in) :: current_s
    real(c_double), intent(in) :: expected_s
    real(c_double), intent(in) :: t
    real(c_double), intent(in) :: current_t
    real(c_double), intent(in) :: expected_t
    logical, intent(out) :: success
    ! Variables outside of signature.
    real(c_double) :: err_s, err_t

    err_s = abs(expected_s - current_s)
    err_t = abs(expected_t - current_t)

    success = ( &
         err_s < abs(expected_s - s) .AND. &
         err_t < abs(expected_t - t) .AND. &
         err_s > 0.5_dp**26 * expected_s .AND. &
         err_t > 0.5_dp**26 * expected_t)

  end subroutine check_closer

  subroutine test_newton_iterate(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(14) :: name
    real(c_double) :: s, t, new_s, new_t
    logical(c_bool) :: converged
    real(c_double) :: quadratic1(2, 3), quadratic2(2, 3)
    real(c_double) :: first_deriv1(2, 2), first_deriv2(2, 2)
    real(c_double) :: line(2, 2), line_first_deriv(2, 1)
    real(c_double) :: second_deriv1(2, 1), second_deriv2(2, 1)
    real(c_double) :: expected_s, expected_t

    case_id = 1
    name = "newton_iterate"

    ! CASE 1: RHS is exactly zero during "Newton simple" iteration.
    !         B1([10922/32768, 10923/32768]) and B2([16383/16384, 1]) are
    !         linearized and when the segments intersect they produce
    !         t = 109217/109216 > 1.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [4.5_dp, 9.0_dp]
    quadratic1(:, 3) = [9.0_dp, 0.0_dp]
    first_deriv1 = 2 * (quadratic1(:, 2:) - quadratic1(:, :2))
    s = 671023103.0_dp / 2013069312.0_dp

    quadratic2(:, 1) = [11.0_dp, 8.0_dp]
    quadratic2(:, 2) = [7.0_dp, 10.0_dp]
    quadratic2(:, 3) = [3.0_dp, 4.0_dp]
    first_deriv2 = 2 * (quadratic2(:, 2:) - quadratic2(:, :2))
    t = 1789394945.0_dp / 1789394944.0_dp

    call newton_iterate(f_simple, s, t, new_s, new_t, converged)
    case_success = ( &
         3 * new_s == 1.0_dp .AND. &
         new_t == 1.0_dp .AND. &
         converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Jacobian becomes singular during "Newton simple" iteration.
    !         B1([5461/8192, 5462/8192]) and B2([2730/8192, 2731/8192]) are
    !         linearized and the segments are parallel. The curves intersect
    !         at the point B1(2/3) = [1/2, 1/2] = B2(1/3) and they have
    !         parallel tangent vectors B1'(2/3) = [3/4, 0] = B2'(1/3).
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [0.375_dp, 0.75_dp]
    quadratic1(:, 3) = [0.75_dp, 0.375_dp]
    first_deriv1 = 2 * (quadratic1(:, 2:) - quadratic1(:, :2))
    s = 10923.0_dp / 16384.0_dp

    quadratic2(:, 1) = [0.25_dp, 0.625_dp]
    quadratic2(:, 2) = [0.625_dp, 0.25_dp]
    quadratic2(:, 3) = 1
    first_deriv2 = 2 * (quadratic2(:, 2:) - quadratic2(:, :2))
    t = 5461.0_dp / 16384.0_dp

    call newton_iterate(f_simple, s, t, new_s, new_t, converged)
    case_success = ( &
         3 * new_s == 2.0_dp + 0.5_dp**14 .AND. &
         3 * new_t == 1.0_dp - 0.5_dp**14 .AND. &
         .NOT. converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: The method converges linearly (rather than quadratically)
    !         during "Newton simple" iteration.
    !         B1([2730/8192, 2731/8192]) and B2([1/2, 2049/4096]) are
    !         linearized and when the segments intersect they produce
    !         t = -1/6 < 0.
    quadratic1(:, 1) = [0.5_dp, 0.125_dp]
    quadratic1(:, 2) = [1.25_dp, -0.25_dp]
    quadratic1(:, 3) = [2.0_dp, 0.5_dp]
    first_deriv1 = 2 * (quadratic1(:, 2:) - quadratic1(:, :2))
    s = 12287.0_dp / 36864.0_dp

    quadratic2(:, 1) = [0.5_dp, -0.125_dp]
    quadratic2(:, 2) = [1.0_dp, 0.125_dp]
    quadratic2(:, 3) = [1.5_dp, -0.125_dp]
    first_deriv2 = 2 * (quadratic2(:, 2:) - quadratic2(:, :2))
    t = 12287.0_dp / 24576.0_dp

    call newton_iterate(f_simple, s, t, new_s, new_t, converged)
    call check_closer( &
         s, new_s, 1.0_dp / 3.0_dp, t, new_t, 0.5_dp, case_success)
    case_success = ( &
         case_success .AND. &
         .NOT. converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: The method converges because relative error is below threshold
    !         during "Newton double" iteration.
    !         B1([2730/8192, 2731/8192]) and B2([2047/4096, 1/2]) are
    !         linearized and when the segments intersect they produce
    !         t = 11/10 > 1.
    !         This re-uses the nodes from CASE 3.

    second_deriv1 = (first_deriv1(:, 2:) - first_deriv1(:, :1))
    second_deriv2 = (first_deriv2(:, 2:) - first_deriv2(:, :1))

    ! NOTE: These ``s-t`` values come after the simple root case exits
    !       due to linear convergence, having started from
    !       s = 6827 / 20480 and t = 20481 / 40960 and updating 4 times.
    s = 109227.0_dp / 327680.0_dp
    t = 327681.0_dp / 655360.0_dp

    call newton_iterate(f_double, s, t, new_s, new_t, converged)
    expected_s = 1.0_dp / 3.0_dp
    case_success = ( &
         expected_s - new_s == spacing(expected_s) .AND. &
         new_t == 0.5_dp .AND. &
         converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: The method converges because relative error is below threshold
    !         during "Newton simple" iteration.
    !         B1([12287/16384, 3/4]) and B2([2457/8192, 2458/8192]) are
    !         linearized and when the segments intersect they produce
    !         s = 33555797/33551701 > 1.
    quadratic1(:, 1) = [1.0_dp, 0.0_dp]
    quadratic1(:, 2) = [-1.0_dp, 0.25_dp]
    quadratic1(:, 3) = [1.0_dp, 0.5_dp]
    first_deriv1 = 2 * (quadratic1(:, 2:) - quadratic1(:, :2))
    s = 25163776.0_dp / 33551701.0_dp

    quadratic2(:, 1) = [-0.125_dp, -0.28125_dp]
    quadratic2(:, 2) = [0.5_dp, 1.28125_dp]
    quadratic2(:, 3) = [1.125_dp, -0.28125_dp]
    first_deriv2 = 2 * (quadratic2(:, 2:) - quadratic2(:, :2))
    t = 41228331827.0_dp / 137427767296.0_dp

    call newton_iterate(f_simple, s, t, new_s, new_t, converged)
    expected_t = 3.0_dp / 10.0_dp
    case_success = ( &
         new_s == 0.75_dp .AND. &
         new_t - expected_t == spacing(expected_t) .AND. &
         converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Jacobian becomes singular during "Newton double" iteration.
    !         The curves are tangent and have the same curvature (i.e.
    !         triple root).
    quadratic1(:, 1) = [12.0_dp, 4.0_dp]
    quadratic1(:, 2) = -4
    quadratic1(:, 3) = [-4.0_dp, 4.0_dp]
    first_deriv1 = 2 * (quadratic1(:, 2:) - quadratic1(:, :2))
    second_deriv1 = (first_deriv1(:, 2:) - first_deriv1(:, :1))
    s = 1125899532873825.0_dp * 0.5_dp**51

    quadratic2(:, 1) = [6.0_dp, 1.0_dp]
    quadratic2(:, 2) = [-2.0_dp, -1.0_dp]
    quadratic2(:, 3) = [-2.0_dp, 1.0_dp]
    first_deriv2 = 2 * (quadratic2(:, 2:) - quadratic2(:, :2))
    second_deriv2 = (first_deriv2(:, 2:) - first_deriv2(:, :1))
    t = 4503596635588511.0_dp * 0.5_dp**53

    call newton_iterate(f_double, s, t, new_s, new_t, converged)
    call check_closer( &
         s, new_s, 0.5_dp, t, new_t, 0.5_dp, case_success)
    case_success = ( &
         case_success .AND. &
         .NOT. converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: The method converges linearly (rather than quadratically)
    !         during "Newton simple" iteration.
    !         B1([16387/32768, 16388/32768]) and B2([8195/16384, 8196/16384])
    !         are linearized and when the segments intersect they produce
    !         s = t = -9/7 < 0.
    !         This re-uses the nodes from CASE 6.

    ! NOTE: These ``s-t`` values come after the simple root case exits
    !       due to linear convergence, having started from
    !       s = 28675 / 57344 and t = 14339 / 28672.
    s = 4503629435419537.0_dp * 0.5_dp**53
    t = 2251829621738309.0_dp * 0.5_dp**52

    call newton_iterate(f_double, s, t, new_s, new_t, converged)
    call check_closer( &
         s, new_s, 0.5_dp, t, new_t, 0.5_dp, case_success)
    case_success = ( &
         case_success .AND. &
         .NOT. converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: A sequence of ``s-t`` parameters converges linearly via
    !         ``sn = tn = 0.5 - k`` ==> ``s{n+1} = t{n+1} = 0.5 - 0.5 k``.
    !         This will be used to show an "early" exit when the root gets
    !         within ``sqrt(machine precision)`` of the actual root
    !         ``s = t = 0.5``.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [0.5_dp, 1.0_dp]
    quadratic1(:, 3) = [1.0_dp, 0.0_dp]
    first_deriv1 = 2 * (quadratic1(:, 2:) - quadratic1(:, :2))
    s = -0.5_dp

    line(:, 1) = [0.0_dp, 0.5_dp]
    line(:, 2) = [1.0_dp, 0.5_dp]
    line_first_deriv = line(:, 2:) - line(:, :1)
    t = -0.5_dp

    call newton_iterate(f_curve_line, s, t, new_s, new_t, converged)
    case_success = ( &
         new_s == 0.5_dp - 0.5_dp**4 .AND. &
         new_t == 0.5_dp - 0.5_dp**4 .AND. &
         .NOT. converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Same as CASE 8, but closer to the solution (but still not close
    !         enough to "converge" prematurely).
    s = 0.5_dp - 0.5_dp**20
    t = 0.5_dp - 0.5_dp**20
    call newton_iterate(f_curve_line, s, t, new_s, new_t, converged)
    case_success = ( &
         new_s == 0.5_dp - 0.5_dp**24 .AND. &
         new_t == 0.5_dp - 0.5_dp**24 .AND. &
         .NOT. converged)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: Same as CASE 8 and 9 and "converges" prematurely, i.e.
    !          ``F(s, t)`` evaluates to zero earlier than it should.
    s = 0.5_dp - 0.5_dp**26
    t = 0.5_dp - 0.5_dp**26
    call newton_iterate(f_curve_line, s, t, new_s, new_t, converged)
    case_success = ( &
         new_s == 0.5_dp - 0.5_dp**28 .AND. &
         new_t == 0.5_dp - 0.5_dp**28 .AND. &
         converged)
    call print_status(name, case_id, case_success, success)

  contains

    ! NOTE: This is a closure around several variables in the
    !       scope above which may change:
    !       * ``quadratic1``
    !       * ``first_deriv1``
    !       * ``quadratic2``
    !       * ``first_deriv2``
    subroutine f_simple(s, t, jacobian, func_val)
      real(c_double), intent(in) :: s
      real(c_double), intent(in) :: t
      real(c_double), intent(out) :: jacobian(2, 2)
      real(c_double), intent(out) :: func_val(2, 1)

      call newton_simple_root( &
           s, 3, quadratic1, first_deriv1, &
           t, 3, quadratic2, first_deriv2, jacobian, func_val)

    end subroutine f_simple

    ! NOTE: This is a closure around several variables in the
    !       scope above which may change:
    !       * ``quadratic1``
    !       * ``first_deriv1``
    !       * ``second_deriv1``
    !       * ``quadratic2``
    !       * ``first_deriv2``
    !       * ``second_deriv2``
    subroutine f_double(s, t, jacobian, func_val)
      real(c_double), intent(in) :: s
      real(c_double), intent(in) :: t
      real(c_double), intent(out) :: jacobian(2, 2)
      real(c_double), intent(out) :: func_val(2, 1)

      call newton_double_root( &
           s, 3, quadratic1, first_deriv1, second_deriv1, &
           t, 3, quadratic2, first_deriv2, second_deriv2, &
           jacobian, func_val)

    end subroutine f_double

    ! NOTE: This is a closure around several variables in the
    !       scope above which may change:
    !       * ``quadratic1``
    !       * ``first_deriv1``
    !       * ``line``
    !       * ``line_first_deriv``
    subroutine f_curve_line(s, t, jacobian, func_val)
      real(c_double), intent(in) :: s
      real(c_double), intent(in) :: t
      real(c_double), intent(out) :: jacobian(2, 2)
      real(c_double), intent(out) :: func_val(2, 1)

      call newton_simple_root( &
           s, 3, quadratic1, first_deriv1, &
           t, 2, line, line_first_deriv, jacobian, func_val)

    end subroutine f_curve_line

  end subroutine test_newton_iterate

  subroutine test_full_newton_nonzero(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(19) :: name
    real(c_double) :: s, t, new_s, new_t
    real(c_double) :: quadratic1(2, 3), quadratic2(2, 3), line(2, 2)
    integer(c_int) :: status
    real(c_double) :: expected_s, expected_t

    case_id = 1
    name = "full_newton_nonzero"

    ! CASE 1: Simple root, converges.
    !         B1([4095/8192, 1/2]) and B2([1365/8192, 1366/8192]) are
    !         linearized and when the segments intersect they produce
    !         s = 24580/24579 > 1.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [0.375_dp, 0.75_dp]
    quadratic1(:, 3) = [0.75_dp, 0.375_dp]
    s = 100675585.0_dp / 201351168.0_dp

    quadratic2(:, 1) = [0.25_dp, 0.5625_dp]
    quadratic2(:, 2) = [0.625_dp, 0.1875_dp]
    quadratic2(:, 3) = [1.0_dp, 0.9375_dp]
    t = 33558529.0_dp / 201351168.0_dp

    call full_newton_nonzero( &
         s, 3, quadratic1, t, 3, quadratic2, new_s, new_t, status)
    expected_s = 0.5_dp
    expected_t = 1.0_dp / 6.0_dp
    case_success = ( &
         new_s - expected_s == spacing(expected_s) .AND. &
         new_t - expected_t == 3 * spacing(expected_t) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Double root, converges.
    !         B1([5461/8192, 5462/8192]) and B2([2730/8192, 2731/8192]) are
    !         linearized and the segments are parallel. The curves intersect
    !         at the point B1(2/3) = [1/2, 1/2] = B2(1/3) and they have
    !         parallel tangent vectors B1'(2/3) = [3/4, 0] = B2'(1/3).
    !         This re-uses the ``quadratic1`` from CASE 1.

    s = 10923.0_dp / 16384.0_dp

    quadratic2(:, 1) = [0.25_dp, 0.625_dp]
    quadratic2(:, 2) = [0.625_dp, 0.25_dp]
    quadratic2(:, 3) = 1
    t = 5461.0_dp / 16384.0_dp

    call full_newton_nonzero( &
         s, 3, quadratic1, t, 3, quadratic2, new_s, new_t, status)
    expected_s = 2.0_dp / 3.0_dp
    expected_t = 1.0_dp / 3.0_dp
    case_success = ( &
         new_s - expected_s == spacing(expected_s) .AND. &
         new_t == expected_t .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Triple root, does not converge.
    !         B1([16382/32768, 16383/32768]) and B2([8190/16384, 8191/16384])
    !         are linearized and when the segments intersect they produce
    !         s = t = 4/3 > 1.
    quadratic1(:, 1) = [12.0_dp, 4.0_dp]
    quadratic1(:, 2) = -4
    quadratic1(:, 3) = [-4.0_dp, 4.0_dp]
    s = 24575.0_dp / 49152.0_dp

    quadratic2(:, 1) = [6.0_dp, 1.0_dp]
    quadratic2(:, 2) = [-2.0_dp, -1.0_dp]
    quadratic2(:, 3) = [-2.0_dp, 1.0_dp]
    t = 12287.0_dp / 24576.0_dp

    call full_newton_nonzero( &
         s, 3, quadratic1, t, 3, quadratic2, new_s, new_t, status)
    call check_closer( &
         s, new_s, 0.5_dp, t, new_t, 0.5_dp, case_success)
    case_success = ( &
         case_success .AND. &
         status == Status_BAD_MULTIPLICITY)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Intersect a line and a curve (and line has special handling
    !         for the second derivative).
    !         B1([5461/16384, 5462/16384]) and B2([0, 1]) are linearized
    !         and when the segments intersect they produce s = -1/3 < 0.
    quadratic1(:, 1) = [0.0_dp, 2.25_dp]
    quadratic1(:, 2) = [1.5_dp, -2.25_dp]
    quadratic1(:, 3) = [3.0_dp, 2.25_dp]
    s = 8191.0_dp / 24576.0_dp

    line(:, 1) = [-0.5_dp, 1.75_dp]
    line(:, 2) = [4.0_dp, -2.75_dp]
    t = 12287.0_dp / 36864.0_dp

    call full_newton_nonzero( &
         s, 3, quadratic1, t, 2, line, new_s, new_t, status)
    expected_s = 1.0_dp / 3.0_dp
    expected_t = 1.0_dp / 3.0_dp
    case_success = ( &
         new_s - expected_s == spacing(expected_s) .AND. &
         new_t == expected_t .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_full_newton_nonzero

  subroutine test_full_newton(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(11) :: name
    real(c_double) :: s, t, new_s, new_t
    real(c_double) :: quadratic1(2, 3), quadratic2(2, 3)
    integer(c_int) :: status
    real(c_double) :: expected_s, expected_t

    case_id = 1
    name = "full_newton"

    ! CASE 1: Both parameters near zero.
    !         B1([0, 1/8192]) and B2([0, 1/8192]) are linearized and the
    !         segments are parallel, and the root is a double root.
    quadratic1(:, 1) = [1.0_dp, 0.0_dp]
    quadratic1(:, 2) = 1
    quadratic1(:, 3) = [0.0_dp, 1.0_dp]
    s = 1.0_dp / 16384.0_dp

    quadratic2(:, 1) = [1.0_dp, 0.0_dp]
    quadratic2(:, 2) = [1.0_dp, 1.5_dp]
    quadratic2(:, 3) = [-0.5_dp, 1.5_dp]
    t = 1.0_dp / 16384.0_dp

    call full_newton( &
         s, 3, quadratic1, t, 3, quadratic2, new_s, new_t, status)
    case_success = ( &
         new_s == 0.0_dp .AND. &
         new_t == 0.0_dp .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``s`` near zero, ``t`` not.
    !         B1([0, 1/8192]) and B2([1/4, 2049/8192]) are linearized and the
    !         segments are parallel, and the root is a double root.
    quadratic1(:, 1) = [1.0_dp, 0.0_dp]
    quadratic1(:, 2) = [1.0_dp, 1.5_dp]
    quadratic1(:, 3) = [-0.5_dp, 1.5_dp]
    s = 1.0_dp / 16384.0_dp

    quadratic2(:, 1) = [0.9375_dp, -0.5625_dp]
    quadratic2(:, 2) = [1.1875_dp, 0.6875_dp]
    quadratic2(:, 3) = [0.4375_dp, 0.9375_dp]
    t = 4097.0_dp / 16384.0_dp

    call full_newton( &
         s, 3, quadratic1, t, 3, quadratic2, new_s, new_t, status)
    case_success = ( &
         new_s == 0.0_dp .AND. &
         new_t == 0.25_dp .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: ``t`` near zero, ``s`` not. We re-use the nodes and parameters
    !         from CASE 2 and just swap the order.
    call full_newton( &
         t, 3, quadratic2, s, 3, quadratic1, new_s, new_t, status)
    case_success = ( &
         new_s == 0.25_dp .AND. &
         new_t == 0.0_dp .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Neither parameter near zero.
    !         B1([6826/8192, 6827/8192]) and B2([1/2, 4097/8192]) are
    !         linearized and when the segments intersect they produce
    !         t = -1/24579 < 0. The root is a simple root.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [0.375_dp, 0.75_dp]
    quadratic1(:, 3) = [0.75_dp, 0.375_dp]
    s = 167792639.0_dp / 201351168.0_dp

    quadratic2(:, 1) = [0.25_dp, 0.5625_dp]
    quadratic2(:, 2) = [0.625_dp, 0.1875_dp]
    quadratic2(:, 3) = [1.0_dp, 0.9375_dp]
    t = 100675583.0_dp / 201351168.0_dp

    expected_s = 5.0_dp / 6.0_dp
    expected_t = 0.5_dp
    call full_newton( &
         s, 3, quadratic1, t, 3, quadratic2, new_s, new_t, status)
    case_success = ( &
         expected_s - new_s == 2 * spacing(expected_s) .AND. &
         expected_t - new_t == 0.5_dp * spacing(expected_t) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_full_newton

  subroutine test_from_linearized(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    type(CurveData) :: curve1, curve2, curve3, curve4
    real(c_double) :: root_nodes1(2, 3), root_nodes2(2, 3)
    real(c_double) :: refined_s, refined_t
    logical(c_bool) :: does_intersect
    integer(c_int) :: status
    real(c_double) :: expected_s, expected_t
    integer :: case_id
    character(15) :: name

    case_id = 1
    name = "from_linearized"

    ! Will hold throughout.
    curve1%start = 0.0_dp
    curve1%end_ = 1.0_dp
    curve2%start = 0.0_dp
    curve2%end_ = 1.0_dp

    ! CASE 1: Basic test of not-very-linearized quadratics that do intersect.
    allocate(curve1%nodes(2, 3))
    curve1%nodes(:, 1) = 0
    curve1%nodes(:, 2) = [0.5_dp, 1.0_dp]
    curve1%nodes(:, 3) = 1

    allocate(curve2%nodes(2, 3))
    curve2%nodes(:, 1) = [0.0_dp, 1.0_dp]
    curve2%nodes(:, 2) = [0.5_dp, 1.0_dp]
    curve2%nodes(:, 3) = [1.0_dp, 0.0_dp]

    call from_linearized( &
         curve1, 3, curve1%nodes, &
         curve2, 3, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         does_intersect .AND. &
         refined_s == 0.5_dp .AND. &
         refined_t == 0.5_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Curves have overlapping bounding boxes, but do not intersect. The
    !         linearized intersection has ``s = 0.5, t == 1.25``, i.e. it is
    !         outside the unit interval. In this case, we fall back to the
    !         convex hulls (which are disjoint).
    curve1%nodes(:, 1) = 0
    curve1%nodes(:, 2) = [0.5_dp, 0.0_dp]
    curve1%nodes(:, 3) = 1

    deallocate(curve2%nodes)
    allocate(curve2%nodes(2, 3))
    curve2%nodes(:, 1) = [1.75_dp, -0.75_dp]
    curve2%nodes(:, 2) = [1.25_dp, -0.75_dp]
    curve2%nodes(:, 3) = [0.75_dp, 0.25_dp]

    call from_linearized( &
         curve1, 3, curve1%nodes, &
         curve2, 3, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         .NOT. does_intersect)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Linearized segments are parallel, so we fall back to the convex
    !         hulls (which are disjoint).
    deallocate(curve1%nodes)
    allocate(curve1%nodes(2, 2))
    curve1%nodes(:, 1) = 0
    curve1%nodes(:, 2) = 1

    ! NOTE: This is a "bad" parameterization of ``y == x``.
    curve2%nodes(:, 1) = 2
    curve2%nodes(:, 2) = 2.5009765625_dp
    curve2%nodes(:, 3) = 3

    call from_linearized( &
         curve1, 2, curve1%nodes, &
         curve2, 3, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         .NOT. does_intersect)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Similar to CASE 2, this is just the linearized segments as lines.
    curve1%nodes(:, 1) = 0
    curve1%nodes(:, 2) = 1

    deallocate(curve2%nodes)
    allocate(curve2%nodes(2, 2))
    curve2%nodes(:, 1) = [1.75_dp, -0.75_dp]
    curve2%nodes(:, 2) = [0.75_dp, 0.25_dp]

    call from_linearized( &
         curve1, 2, curve1%nodes, &
         curve2, 2, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         .NOT. does_intersect)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Parallel lines (should have been handled already by
    !         ``check_lines``).
    curve1%nodes(:, 1) = 0
    curve1%nodes(:, 2) = 1

    curve2%nodes(:, 1) = [0.0_dp, 1.0_dp]
    curve2%nodes(:, 2) = [1.0_dp, 2.0_dp]

    call from_linearized( &
         curve1, 2, curve1%nodes, &
         curve2, 2, curve2%nodes, &
         refined_s, refined_t, does_intersect, status)
    case_success = ( &
         status == Status_SUCCESS .AND. &
         .NOT. does_intersect)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Curves have parallel linearized segments caused by a tangent
    !         intersection.
    curve3%start = 5461.0_dp / 8192
    curve3%end_ = 5462.0_dp / 8192
    allocate(curve3%nodes(2, 3))
    root_nodes1(:, 1) = 0
    root_nodes1(:, 2) = [0.375_dp, 0.75_dp]
    root_nodes1(:, 3) = [0.75_dp, 0.375_dp]
    call specialize_curve( &
         3, 2, root_nodes1, curve3%start, curve3%end_, curve3%nodes)

    curve4%start = 2730.0_dp / 8192
    curve4%end_ = 2731.0_dp / 8192
    allocate(curve4%nodes(2, 3))
    root_nodes2(:, 1) = [0.25_dp, 0.625_dp]
    root_nodes2(:, 2) = [0.625_dp, 0.25_dp]
    root_nodes2(:, 3) = 1
    call specialize_curve( &
         3, 2, root_nodes2, curve4%start, curve4%end_, curve4%nodes)

    call from_linearized( &
         curve3, 3, root_nodes1, &
         curve4, 3, root_nodes2, &
         refined_s, refined_t, does_intersect, status)
    expected_s = 2.0_dp / 3.0_dp
    expected_t = 1.0_dp / 3.0_dp
    case_success = ( &
         status == Status_SUCCESS .AND. &
         does_intersect .AND. &
         refined_s - expected_s == spacing(expected_s) .AND. &
         refined_t == expected_t)
    call print_status(name, case_id, case_success, success)

  end subroutine test_from_linearized

  subroutine test_bbox_line_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: unit_square(2, 4)
    real(c_double) :: line_start(2), line_end(2)
    integer(c_int) :: enum_
    integer :: case_id
    character(19) :: name

    case_id = 1
    name = "bbox_line_intersect"
    unit_square(:, 1) = 0
    unit_square(:, 2) = [1.0_dp, 0.0_dp]
    unit_square(:, 3) = 1
    unit_square(:, 4) = [0.0_dp, 1.0_dp]

    ! CASE 1: Line starts inside bounding box.
    line_start = 0.5_dp
    line_end = [0.5_dp, 1.5_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Line ends (but does not start) in bounding box.
    line_start = [-1.0_dp, 0.5_dp]
    line_end = 0.5_dp
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Line segment "pierces" bbox from the bottom.
    line_start = [0.5_dp, -0.5_dp]
    line_end = [0.5_dp, 1.5_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Line segment "pierces" bbox from the right.
    line_start = [-0.5_dp, 0.5_dp]
    line_end = [1.5_dp, 0.5_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Line segment "pierces" bbox from the top.
    line_start = [-0.25_dp, 0.5_dp]
    line_end = [0.5_dp, 1.25_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_INTERSECTION)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Line segment is disjoint from bbox.
    line_start = 2
    line_end = [2.0_dp, 5.0_dp]
    call bbox_line_intersect( &
         4, unit_square, line_start, line_end, enum_)
    case_success = (enum_ == BoxIntersectionType_DISJOINT)
    call print_status(name, case_id, case_success, success)

  end subroutine test_bbox_line_intersect

  subroutine test_check_lines(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 3), nodes2(2, 3), nodes3(2, 2), nodes4(2, 2)
    logical(c_bool) :: both_linear
    logical(c_bool) :: coincident
    real(c_double), allocatable :: intersections(:, :)
    integer(c_int) :: num_intersections
    integer :: case_id
    character(19) :: name

    case_id = 1
    name = "check_lines"

    ! CASE 1: Test first curve not linearized.
    nodes1(:, 1) = 0
    nodes1(:, 2) = 1
    nodes1(:, 3) = [2.0_dp, 0.0_dp]
    nodes2(:, 1) = [0.0_dp, 1.0_dp]
    nodes2(:, 2) = [1.0_dp, 0.0_dp]
    nodes2(:, 3) = [2.0_dp, 1.0_dp]

    call check_lines( &
         3, nodes1, 3, nodes2, &
         both_linear, coincident, intersections, num_intersections)
    case_success = ( &
         .NOT. both_linear .AND. &
         .NOT. coincident .AND. &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Test second curve not linearized.
    nodes3(:, 1) = 0
    nodes3(:, 2) = 1

    call check_lines( &
         2, nodes3, 3, nodes2, &
         both_linear, coincident, intersections, num_intersections)
    case_success = ( &
         .NOT. both_linear .AND. &
         .NOT. coincident .AND. &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Lines aren't parallel, but do not intersect in their parameter
    !         ranges.
    nodes4(:, 1) = [2.0_dp, 3.0_dp]
    nodes4(:, 2) = [3.0_dp, 2.0_dp]

    call check_lines( &
         2, nodes3, 2, nodes4, &
         both_linear, coincident, intersections, num_intersections)
    case_success = ( &
         both_linear .AND. &
         .NOT. coincident .AND. &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Lines do intersect (aren't parallel). Also, ``intersections``
    !         gets allocated.
    nodes4(:, 1) = [0.0_dp, 1.0_dp]
    nodes4(:, 2) = [1.0_dp, 0.0_dp]

    call check_lines( &
         2, nodes3, 2, nodes4, &
         both_linear, coincident, intersections, num_intersections)
    case_success = ( &
         both_linear .AND. &
         .NOT. coincident .AND. &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         all(intersections(:, 1) == 0.5_dp))
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Lines do intersect (aren't parallel), but ``intersections``
    !         is already allocated.
    nodes4(:, 1) = [0.0_dp, 1.0_dp]
    nodes4(:, 2) = [2.0_dp, -1.0_dp]

    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]))
    call check_lines( &
         2, nodes3, 2, nodes4, &
         both_linear, coincident, intersections, num_intersections)
    case_success = ( &
         case_success .AND. &
         both_linear .AND. &
         .NOT. coincident .AND. &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         all(intersections(:, 1) == [0.5_dp, 0.25_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Lines are parallel, but do not shared a coincident segment.
    !         This also checks that even in the "failure" case (i.e. not
    !         coincident), ``intersections`` will get allocated. It **also**
    !         shows that it increases the size when too small.
    nodes4(:, 1) = [0.0_dp, 1.0_dp]
    nodes4(:, 2) = [1.0_dp, 2.0_dp]

    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]))
    call check_lines( &
         2, nodes3, 2, nodes4, &
         both_linear, coincident, intersections, num_intersections)
    case_success = ( &
         case_success .AND. &
         both_linear .AND. &
         .NOT. coincident .AND. &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 0)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Lines are parallel and do share a coincident segment (and
    !         ``intersections`` is "big enough").
    nodes4(:, 1) = [0.5_dp, 0.5_dp]
    nodes4(:, 2) = [4.5_dp, 4.5_dp]

    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]))
    call check_lines( &
         2, nodes3, 2, nodes4, &
         both_linear, coincident, intersections, num_intersections)
    case_success = ( &
         case_success .AND. &
         both_linear .AND. &
         coincident .AND. &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 2 .AND. &
         all(intersections(:, 1) == [0.5_dp, 0.0_dp]) .AND. &
         all(intersections(:, 2) == [1.0_dp, 0.125_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Same as CASE 8, but ``intersections`` does not start
    !         out allocated.
    deallocate(intersections)
    call check_lines( &
         2, nodes3, 2, nodes4, &
         both_linear, coincident, intersections, num_intersections)
    case_success = ( &
         case_success .AND. &
         both_linear .AND. &
         coincident .AND. &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 2 .AND. &
         all(intersections(:, 1) == [0.5_dp, 0.0_dp]) .AND. &
         all(intersections(:, 2) == [1.0_dp, 0.125_dp]))
    call print_status(name, case_id, case_success, success)

  end subroutine test_check_lines

  subroutine test_add_intersection(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    real(c_double) :: s_val, t_val
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

    ! CASE 6: ``s``-value is below the ``ZERO_THRESHOLD``.
    intersections(1, 1) = 0.0_dp
    intersections(2, 1) = 0.75_dp
    num_intersections = 1

    s_val = 0.5_dp**46
    call add_intersection( &
         s_val, 0.75_dp, num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 1)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: ``s``-value is below the ``ZERO_THRESHOLD``.
    intersections(1, 1) = 0.125_dp
    intersections(2, 1) = 0.0_dp
    num_intersections = 1

    t_val = 0.5_dp**45
    call add_intersection( &
         0.125_dp, t_val, num_intersections, intersections)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 1)
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_intersection

  subroutine test_add_from_linearized(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer(c_int) :: status
    real(c_double) :: nodes1(2, 2), nodes2(2, 3)
    real(c_double) :: root_nodes1(2, 3), root_nodes2(2, 3)
    type(CurveData) :: first, second
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    integer :: case_id
    character(19) :: name

    case_id = 1
    name = "add_from_linearized"
    num_intersections = 0

    ! CASE 1: Lines that intersect. Since lines, there are 0 subdivisions.
    first%start = 0.0_dp
    first%end_ = 1.0_dp
    nodes1(:, 1) = 0
    nodes1(:, 2) = 1
    first%nodes = nodes1

    second%start = 0.0_dp
    second%end_ = 1.0_dp
    nodes1(:, 1) = [0.0_dp, 1.0_dp]
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    second%nodes = nodes1

    call add_from_linearized( &
         first, first%nodes, second, second%nodes, &
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
    nodes1(:, 1) = 0
    nodes1(:, 2) = 1
    first%nodes = nodes1

    second%start = 0.0_dp
    second%end_ = 1.0_dp
    nodes1(:, 1) = [3.0_dp, 0.0_dp]
    nodes1(:, 2) = [2.0_dp, 1.0_dp]
    second%nodes = nodes1

    call add_from_linearized( &
         first, first%nodes, second, second%nodes, &
         num_intersections, intersections, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Quadratic curves that intersect after many (12) subdivisions.
    root_nodes1(:, 1) = [0.25_dp, 0.4375_dp]
    root_nodes1(:, 2) = [0.625_dp, 1.0_dp]
    root_nodes1(:, 3) = 1
    first%start = 1365.0_dp / 4096.0_dp
    first%end_ = 1366.0_dp / 4096.0_dp
    ! NOTE: This is the result of
    !       call specialize_curve( &
    !            3, 2, root_nodes1, first%start, first%end_, ...)
    nodes2(:, 1) = [134201344.0_dp, 201310207.0_dp]
    nodes2(:, 2) = [134225920.0_dp, 201334786.0_dp]
    nodes2(:, 3) = [134250496.0_dp, 201359356.0_dp]
    first%nodes = 0.5_dp**28 * nodes2

    root_nodes2(:, 1) = [0.0_dp, 1.0_dp]
    root_nodes2(:, 2) = [0.375_dp, 1.0_dp]
    root_nodes2(:, 3) = [0.75_dp, 0.4375_dp]
    second%start = 2730.0_dp / 4096.0_dp
    second%end_ = 2731.0_dp / 4096.0_dp
    ! NOTE: This is the result of
    !       call specialize_curve( &
    !            3, 2, root_nodes2, second%start, second%end_, ...)
    nodes2(:, 1) = [134184960.0_dp, 201359356.0_dp]
    nodes2(:, 2) = [134209536.0_dp, 201334786.0_dp]
    nodes2(:, 3) = [134234112.0_dp, 201310207.0_dp]
    second%nodes = 0.5_dp**28 * nodes2

    call add_from_linearized( &
         first, root_nodes1, second, root_nodes2, &
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

    ! CASE 4: A line that is tangent to a curve, which requires a "full"
    !         Newton iteration.
    first%start = 5461.0_dp / 16384.0_dp
    first%end_ = 5462.0_dp / 16384.0_dp
    root_nodes1(:, 1) = [0.0_dp, 2.25_dp]
    root_nodes1(:, 2) = [1.5_dp, -2.25_dp]
    root_nodes1(:, 3) = [3.0_dp, 2.25_dp]
    call specialize_curve( &
         3, 2, root_nodes1, first%start, first%end_, first%nodes)

    second%start = 0.0_dp
    second%end_ = 1.0_dp
    nodes1(:, 1) = [-0.5_dp, 1.75_dp]
    nodes1(:, 2) = [4.0_dp, -2.75_dp]
    second%nodes = nodes1

    call add_from_linearized( &
         first, root_nodes1, second, second%nodes, &
         num_intersections, intersections, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         status == Status_SUCCESS .AND. &
         3.0_dp * intersections(1, 1) == 1.0_dp .AND. &
         3.0_dp * intersections(2, 1) == 1.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Same as CASE 4, with arguments swapped.
    call add_from_linearized( &
         second, second%nodes, first, root_nodes1, &
         num_intersections, intersections, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         status == Status_SUCCESS .AND. &
         3.0_dp * intersections(1, 1) == 1.0_dp .AND. &
         3.0_dp * intersections(2, 1) == 1.0_dp)
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_from_linearized

  subroutine test_endpoint_check(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    type(CurveData) :: first, second
    real(c_double) :: node_first(2), node_second(2)
    real(c_double) :: s, t
    real(c_double) :: nodes1(2, 2), nodes2(2, 3)
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
    nodes1(:, 1) = 0
    nodes1(:, 2) = 1
    first%nodes = nodes1
    s = 1.0_dp
    node_first = first%nodes(:, 2)

    nodes1(:, 1) = 1
    nodes1(:, 2) = [2.0_dp, 1.0_dp]
    second%nodes = nodes1
    t = 0.0_dp
    node_second = second%nodes(:, 1)

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
    nodes2(:, 1) = 0
    nodes2(:, 2) = [0.25_dp, 0.5_dp]
    nodes2(:, 3) = 0.5_dp
    first%nodes = nodes2
    node_first = first%nodes(:, 3)

    second%start = 0.5_dp
    second%end_ = 1.0_dp
    nodes2(:, 1) = 0.5_dp
    nodes2(:, 2) = [0.5_dp, 0.0_dp]
    nodes2(:, 3) = [1.0_dp, -0.5_dp]
    second%nodes = nodes2
    node_second = second%nodes(:, 1)

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
    real(c_double) :: nodes1(2, 3), nodes2(2, 2)
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    integer :: case_id
    character(25) :: name

    case_id = 1
    name = "tangent_bbox_intersection"
    num_intersections = 0

    ! CASE 1: Standard case of two quadratics that touch at a single endpoint.
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 2.0_dp]
    nodes1(:, 3) = [2.0_dp, 0.0_dp]
    first%nodes = nodes1

    nodes1(:, 1) = [2.0_dp, 0.0_dp]
    nodes1(:, 2) = [3.0_dp, 2.0_dp]
    nodes1(:, 3) = [4.0_dp, 0.0_dp]
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
    nodes1(:, 1) = 0
    nodes1(:, 2) = [-1.0_dp, 0.5_dp]
    nodes1(:, 3) = [0.0_dp, 1.0_dp]
    first%nodes = nodes1

    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.5_dp]
    nodes1(:, 3) = [0.0_dp, 1.0_dp]
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
    nodes2(:, 1) = 0
    nodes2(:, 2) = [2.0_dp, 1.0_dp]
    first%nodes = nodes2

    nodes2(:, 1) = [0.5_dp, 1.0_dp]
    nodes2(:, 2) = [2.5_dp, 2.0_dp]
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
    allocate(curve1%nodes(2, 3))
    curve1%nodes(:, 1) = [0.0_dp, 1.0_dp]
    curve1%nodes(:, 2) = [1.0_dp, 2.0_dp]
    curve1%nodes(:, 3) = [2.0_dp, 1.0_dp]
    allocate(curve2%nodes(2, 2))
    curve2%nodes(:, 1) = [0.0_dp, 1.25_dp]
    curve2%nodes(:, 2) = [2.0_dp, 1.25_dp]
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
    real(c_double) :: fixed_quadratic1(2, 3), fixed_quadratic2(2, 3)
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
    fixed_quadratic1(:, 1) = [0.25_dp, 0.4375_dp]
    fixed_quadratic1(:, 2) = [0.625_dp, 1.0_dp]
    fixed_quadratic1(:, 3) = 1
    ! NOTE: ``fixed_quadratic2`` is a specialization of
    !       [0, 1], [1/2, 1], [1, 0] onto the interval [0, 3/4].
    fixed_quadratic2(:, 1) = [0.0_dp, 1.0_dp]
    fixed_quadratic2(:, 2) = [0.375_dp, 1.0_dp]
    fixed_quadratic2(:, 3) = [0.75_dp, 0.4375_dp]

    fixed_line1(:, 1) = 0
    fixed_line1(:, 2) = 1

    fixed_line2(:, 1) = [0.0_dp, 1.0_dp]
    fixed_line2(:, 2) = [1.0_dp, 0.0_dp]

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

    ! CASE 5: Disjoint bounding boxes.
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve with a line.
    candidates(1, 1)%nodes = fixed_quadratic1
    ! Populate the "second" curve with a line.
    allocate(candidates(2, 1)%nodes(2, 2))
    candidates(2, 1)%nodes(:, 1) = [1.0_dp, 1.25_dp]
    candidates(2, 1)%nodes(:, 2) = [0.0_dp, 2.0_dp]

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

    ! CASE 6: Tangent bounding boxes (**with** an intersection), noting
    !         that tangency is only allowed for a pair with both curves
    !         can't be linearized (since tangency can be resolved in the
    !         linearized case by checking endpoints).
    num_candidates = 1
    allocate(candidates(2, num_candidates))
    ! Populate the "first" curve.
    allocate(candidates(1, 1)%nodes(2, 3))
    candidates(1, 1)%nodes(:, 1) = 0
    candidates(1, 1)%nodes(:, 2) = [0.5_dp, 1.0_dp]
    candidates(1, 1)%nodes(:, 3) = [1.0_dp, 0.0_dp]
    ! [0 <= x <= 1.0], [0.0 <= y <= 1.0]
    ! Populate the "second" curve.
    allocate(candidates(2, 1)%nodes(2, 3))
    candidates(2, 1)%nodes(:, 1) = [1.0_dp, 0.0_dp]
    candidates(2, 1)%nodes(:, 2) = [1.5_dp, 0.5_dp]
    candidates(2, 1)%nodes(:, 3) = [2.0_dp, -0.25_dp]
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

  subroutine test_make_same_degree(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    real(c_double) :: nodes1(2, 3), nodes2(2, 3), nodes3(2, 2), nodes4(2, 4)
    integer(c_int) :: num_nodes
    real(c_double), allocatable :: elevated1(:, :)
    real(c_double), allocatable :: elevated2(:, :)
    character(16) :: name

    case_id = 1
    name = "make_same_degree"

    ! CASE 1: Same degree.
    nodes1(:, 1) = [0.0_dp, 1.0_dp]
    nodes1(:, 2) = [1.0_dp, 2.0_dp]
    nodes1(:, 3) = [2.0_dp, 1.0_dp]
    nodes2(:, 1) = [1.0_dp, 2.0_dp]
    nodes2(:, 2) = [2.0_dp, 1.0_dp]
    nodes2(:, 3) = [4.0_dp, 2.0_dp]
    call make_same_degree( &
         3, nodes1, 3, nodes2, num_nodes, elevated1, elevated2)
    case_success = ( &
         num_nodes == 3 .AND. &
         allocated(elevated1) .AND. &
         all(shape(elevated1) == [2, 3]) .AND. &
         all(elevated1 == nodes1) .AND. &
         allocated(elevated2) .AND. &
         all(shape(elevated2) == [2, 3]) .AND. &
         all(elevated2 == nodes2))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Elevate once (first argument).
    nodes3(:, 1) = 0
    nodes3(:, 2) = 1
    nodes1(:, 1) = [1.0_dp, 2.0_dp]
    nodes1(:, 2) = 2
    nodes1(:, 3) = 0
    call make_same_degree( &
         2, nodes3, 3, nodes1, num_nodes, elevated1, elevated2)
    case_success = ( &
         num_nodes == 3 .AND. &
         allocated(elevated1) .AND. &
         all(shape(elevated1) == [2, 3]) .AND. &
         all(elevated1(:, 1) == 0) .AND. &
         all(elevated1(:, 2) == 0.5_dp) .AND. &
         all(elevated1(:, 3) == 1) .AND. &
         allocated(elevated2) .AND. &
         all(shape(elevated2) == [2, 3]) .AND. &
         all(elevated2 == nodes1))
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Same as CASE 2, but swap arguments.
    call make_same_degree( &
         3, nodes1, 2, nodes3, num_nodes, elevated1, elevated2)
    case_success = ( &
         num_nodes == 3 .AND. &
         allocated(elevated1) .AND. &
         all(shape(elevated1) == [2, 3]) .AND. &
         all(elevated1 == nodes1) .AND. &
         allocated(elevated2) .AND. &
         all(shape(elevated2) == [2, 3]) .AND. &
         all(elevated2(:, 1) == 0) .AND. &
         all(elevated2(:, 2) == 0.5_dp) .AND. &
         all(elevated2(:, 3) == 1))
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Elevate twice (first argument).
    nodes3(:, 1) = 0
    nodes3(:, 2) = 3
    nodes4(:, 1) = [0.0_dp, 1.0_dp]
    nodes4(:, 2) = [1.0_dp, 2.0_dp]
    nodes4(:, 3) = [3.0_dp, 2.0_dp]
    nodes4(:, 4) = [4.0_dp, 2.0_dp]
    call make_same_degree( &
         2, nodes3, 4, nodes4, num_nodes, elevated1, elevated2)
    case_success = ( &
         num_nodes == 4 .AND. &
         allocated(elevated1) .AND. &
         all(shape(elevated1) == [2, 4]) .AND. &
         all(elevated1(:, 1) == 0) .AND. &
         all(elevated1(:, 2) == 1) .AND. &
         all(elevated1(:, 3) == 2) .AND. &
         all(elevated1(:, 4) == 3) .AND. &
         allocated(elevated2) .AND. &
         all(shape(elevated2) == [2, 4]) .AND. &
         all(elevated2 == nodes4))
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Same as CASE 4, but swap arguments.
    call make_same_degree( &
         4, nodes4, 2, nodes3, num_nodes, elevated1, elevated2)
    case_success = ( &
         num_nodes == 4 .AND. &
         allocated(elevated1) .AND. &
         all(shape(elevated1) == [2, 4]) .AND. &
         all(elevated1 == nodes4) .AND. &
         allocated(elevated2) .AND. &
         all(shape(elevated2) == [2, 4]) .AND. &
         all(elevated2(:, 1) == 0) .AND. &
         all(elevated2(:, 2) == 1) .AND. &
         all(elevated2(:, 3) == 2) .AND. &
         all(elevated2(:, 4) == 3))
    call print_status(name, case_id, case_success, success)

  end subroutine test_make_same_degree

  subroutine test_add_coincident_parameters(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: linear1(2, 2)
    real(c_double) :: quadratic1(2, 3), quadratic2(2, 3)
    real(c_double) :: cubic1(2, 4), cubic2(2, 4)
    real(c_double) :: params(2, 8), start, end_
    real(c_double) :: expected(2, 2)
    integer(c_int) :: num_intersections
    real(c_double), allocatable :: intersections(:, :)
    logical(c_bool) :: coincident
    integer(c_int) :: i
    integer :: case_id
    character(25) :: name

    case_id = 1
    name = "add_coincident_parameters"
    num_intersections = 0
    allocate(intersections(2, 2))

    ! CASE 1: Actual curve degrees differ (i.e. one is a line, the other is
    !         a quadratic).
    linear1(:, 1) = 0
    linear1(:, 2) = 1
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = 1
    quadratic1(:, 3) = [2.0_dp, 0.0_dp]

    call add_coincident_parameters( &
         2, linear1, 3, quadratic1, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Actual curve degrees are the same, but one is elevated (and is
    !         the same curve segment).
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = 3
    quadratic1(:, 3) = [6.0_dp, 0.0_dp]
    cubic1(:, 1) = 0
    cubic1(:, 2) = 2
    cubic1(:, 3) = [4.0_dp, 2.0_dp]
    cubic1(:, 4) = [6.0_dp, 0.0_dp]

    call add_coincident_parameters( &
         3, quadratic1, 4, cubic1, &
         num_intersections, intersections, coincident)
    case_success = ( &
         all(shape(intersections) == [2, 2]) .AND. &
         all(intersections(:, 1) == 0) .AND. &
         all(intersections(:, 2) == 1) .AND. &
         num_intersections == 2 .AND. &
         coincident)
    call print_status(name, case_id, case_success, success)
    num_intersections = 0

    ! CASE 3: Fail due to invalid point (i.e curve has a self-crossing).
    cubic1(:, 1) = [0.0_dp, 16.0_dp]
    cubic1(:, 2) = [-8.0_dp, 0.0_dp]
    cubic1(:, 3) = [8.0_dp, 8.0_dp]
    cubic1(:, 4) = [-6.0_dp, 13.0_dp]
    linear1(:, 1) = [-2.0_dp, 11.0_dp]
    linear1(:, 2) = [-2.0_dp, 16.0_dp]
    call add_coincident_parameters( &
         4, cubic1, 2, linear1, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Same as CASE 3, with arguments swapped.
    call add_coincident_parameters( &
         2, linear1, 4, cubic1, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Segments are part of the same algebraic curve, but the segments
    !         only touch at endpoints (i.e. they don't overlap).
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [1.0_dp, 2.0_dp]
    quadratic1(:, 3) = [3.0_dp, 2.0_dp]
    params(1, :4) = [1.0_dp, 2.0_dp, -1.0_dp, 0.0_dp]
    params(2, :4) = [2.0_dp, 1.0_dp, 0.0_dp, -1.0_dp]

    case_success = .TRUE.
    do i = 1, 4
       call specialize_curve( &
            3, 2, quadratic1, params(1, i), params(2, i), quadratic2)
       call add_coincident_parameters( &
            3, quadratic1, 3, quadratic2, &
            num_intersections, intersections, coincident)
       case_success = ( &
            case_success .AND. &
            num_intersections == 0 .AND. &
            .NOT. coincident)
    end do
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Segments are part of the same algebraic curve, but the segments
    !         are totally disjoint.
    cubic1(:, 1) = [1.0_dp, 0.0_dp]
    cubic1(:, 2) = [1.0_dp, 2.0_dp]
    cubic1(:, 3) = [3.0_dp, 2.0_dp]
    cubic1(:, 4) = [4.0_dp, 0.0_dp]
    params(1, :4) = [1.5_dp, 2.0_dp, -1.0_dp, -0.5_dp]
    params(2, :4) = [2.0_dp, 1.5_dp, -0.5_dp, -1.0_dp]

    case_success = .TRUE.
    do i = 1, 4
       call specialize_curve( &
            4, 2, cubic1, params(1, i), params(2, i), cubic2)
       call add_coincident_parameters( &
            4, cubic1, 4, cubic2, &
            num_intersections, intersections, coincident)
       case_success = ( &
            case_success .AND. &
            num_intersections == 0 .AND. &
            .NOT. coincident)
    end do
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Curve endpoints touch, but curves are not coincident.
    quadratic1(:, 1) = [1.0_dp, 8.0_dp]
    quadratic1(:, 2) = [3.0_dp, 6.0_dp]
    quadratic1(:, 3) = [3.0_dp, 5.0_dp]
    quadratic2(:, 1) = [3.0_dp, 5.0_dp]
    quadratic2(:, 2) = [3.0_dp, -1.0_dp]
    quadratic2(:, 3) = [4.0_dp, -4.0_dp]
    call add_coincident_parameters( &
         3, quadratic1, 3, quadratic2, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Identical curves.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [2.0_dp, 3.0_dp]
    quadratic1(:, 3) = 5
    call add_coincident_parameters( &
         3, quadratic1, 3, quadratic1, &
         num_intersections, intersections, coincident)
    case_success = ( &
         all(shape(intersections) == [2, 2]) .AND. &
         all(intersections(:, 1) == 0) .AND. &
         all(intersections(:, 2) == 1) .AND. &
         num_intersections == 2 .AND. &
         coincident)
    call print_status(name, case_id, case_success, success)
    num_intersections = 0

    ! CASE 9: Segments are part of the same algebraic curve, one segment
    !         is interior to the other **AND** they touch at endpoints.
    quadratic1(:, 1) = [4.0_dp, 1.0_dp]
    quadratic1(:, 2) = [6.0_dp, 3.0_dp]
    quadratic1(:, 3) = [2.0_dp, 1.0_dp]
    params(1, :) = [ &
         0.0_dp, 0.5_dp, 0.5_dp, 1.0_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp]
    params(2, :) = [ &
         0.5_dp, 0.0_dp, 1.0_dp, 0.5_dp, 2.0_dp, 0.0_dp, 1.0_dp, -1.0_dp]
    case_success = .TRUE.
    do i = 1, 8
       start = params(1, i)
       end_ = params(2, i)
       call specialize_curve( &
            3, 2, quadratic1, start, end_, quadratic2)
       call add_coincident_parameters( &
            3, quadratic1, 3, quadratic2, &
            num_intersections, intersections, coincident)
       if (start == 2.0_dp .OR. end_ == 2.0_dp) then
          expected(:, 1) = [0.0_dp, -start / (end_ - start)]
          expected(:, 2) = [1.0_dp, (end_ - 1.0_dp) / (end_ - start)]
       else if (start == -1.0_dp .OR. end_ == -1.0_dp) then
          expected(:, 1) = [0.0_dp, -start / (end_ - start)]
          expected(:, 2) = [1.0_dp, (1.0_dp - start) / (end_ - start)]
       else
          expected(:, 1) = [start, 0.0_dp]
          expected(:, 2) = [end_, 1.0_dp]
       end if
       case_success = ( &
            case_success .AND. &
            all(shape(intersections) == [2, 2]) .AND. &
            all(intersections == expected) .AND. &
            num_intersections == 2 .AND. &
            coincident)
       num_intersections = 0
    end do
    call print_status(name, case_id, case_success, success)

    ! CASE 10: Segments are part of the same algebraic curve, one segment
    !          is fully interior to the other (i.e. same as CASE 9 but with
    !          no endpoint contact).
    quadratic1(:, 1) = -1
    quadratic1(:, 2) = [0.0_dp, 2.0_dp]
    quadratic1(:, 3) = [2.0_dp, 0.0_dp]
    params(1, :4) = [0.25_dp, 0.75_dp, -0.5_dp, 1.5_dp]
    params(2, :4) = [0.75_dp, 0.25_dp, 1.5_dp, -0.5_dp]
    case_success = .TRUE.
    do i = 1, 4
       start = params(1, i)
       end_ = params(2, i)
       call specialize_curve( &
            3, 2, quadratic1, start, end_, quadratic2)
       call add_coincident_parameters( &
            3, quadratic1, 3, quadratic2, &
            num_intersections, intersections, coincident)
       if (start == -0.5_dp .OR. end_ == -0.5_dp) then
          expected(:, 1) = [0.0_dp, (end_ - 1.0_dp) / (end_ - start)]
          expected(:, 2) = [1.0_dp, end_ / (end_ - start)]
       else
          expected(:, 1) = [start, 0.0_dp]
          expected(:, 2) = [end_, 1.0_dp]
       end if
       case_success = ( &
            case_success .AND. &
            all(shape(intersections) == [2, 2]) .AND. &
            all(intersections == expected) .AND. &
            num_intersections == 2 .AND. &
            coincident)
       num_intersections = 0
    end do
    call print_status(name, case_id, case_success, success)

    ! CASE 11: Segments are part of the same algebraic curve and they overlap
    !          in a staggered fashion. I.e curve A starts somewhere in the
    !          middle of curve B and curve B ends somewhere in the middle of
    !          curve A.
    cubic1(:, 1) = [0.0_dp, -1.0_dp]
    cubic1(:, 2) = [1.0_dp, 2.0_dp]
    cubic1(:, 3) = [1.0_dp, 0.0_dp]
    cubic1(:, 4) = [3.0_dp, 2.0_dp]
    params(1, :4) = [0.5_dp, 1.5_dp, -0.5_dp, 0.5_dp]
    params(2, :4) = [1.5_dp, 0.5_dp, 0.5_dp, -0.5_dp]
    case_success = .TRUE.
    do i = 1, 4
       start = params(1, i)
       end_ = params(2, i)
       call specialize_curve( &
            4, 2, cubic1, start, end_, cubic2)
       call add_coincident_parameters( &
            4, cubic1, 4, cubic2, &
            num_intersections, intersections, coincident)
       if (start == 1.5_dp .OR. end_ == 1.5_dp) then
          expected(:, 1) = [0.5_dp, start - 0.5_dp]
          expected(:, 2) = [1.0_dp, 0.5_dp]
       else
          expected(:, 1) = [0.0_dp, 0.5_dp]
          expected(:, 2) = [0.5_dp, end_ + 0.5_dp]
       end if
       case_success = ( &
            case_success .AND. &
            all(shape(intersections) == [2, 2]) .AND. &
            all(intersections == expected) .AND. &
            num_intersections == 2 .AND. &
            coincident)
       num_intersections = 0
    end do
    call print_status(name, case_id, case_success, success)

    ! CASE 12: Exactly one endpoint from each curve lies on the other curve,
    !          but the curves are not coincident.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [16.0_dp, 0.0_dp]
    quadratic1(:, 3) = 16
    quadratic2(:, 1) = [7.0_dp, 1.0_dp]
    quadratic2(:, 2) = [7.0_dp, 17.0_dp]
    quadratic2(:, 3) = [23.0_dp, 17.0_dp]
    call add_coincident_parameters( &
         3, quadratic1, 3, quadratic2, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

    ! CASE 13: Both endpoints from one curve lie on the other curve,
    !          but the curves are not coincident.
    cubic1(:, 1) = 0
    cubic1(:, 2) = 32
    cubic1(:, 3) = [96.0_dp, 32.0_dp]
    cubic1(:, 4) = [128.0_dp, 0.0_dp]
    cubic2(:, 1) = [29.0_dp, 18.0_dp]
    cubic2(:, 2) = [49.0_dp, 38.0_dp]
    cubic2(:, 3) = [79.0_dp, 38.0_dp]
    cubic2(:, 4) = [99.0_dp, 18.0_dp]
    call add_coincident_parameters( &
         4, cubic1, 4, cubic2, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

    ! CASE 14: Same as CASE 13, with arguments swapped.
    call add_coincident_parameters( &
         4, cubic2, 4, cubic1, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

    ! CASE 15: Exactly one endpoint from curve A lies on curve B, but not at
    !          one of the endpoints of curve B.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = 2
    quadratic1(:, 3) = [4.0_dp, 0.0_dp]
    quadratic2(:, 1) = [2.0_dp, 1.0_dp]
    quadratic2(:, 2) = [4.0_dp, 3.0_dp]
    quadratic2(:, 3) = [6.0_dp, 2.0_dp]
    call add_coincident_parameters( &
         3, quadratic1, 3, quadratic2, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

    ! CASE 16: Same as CASE 15, with arguments swapped.
    call add_coincident_parameters( &
         3, quadratic2, 3, quadratic1, &
         num_intersections, intersections, coincident)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident)
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_coincident_parameters

  subroutine test_all_intersections(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double), allocatable :: intersections(:, :)
    real(c_double) :: linear1(2, 2), linear2(2, 2)
    real(c_double) :: quadratic1(2, 3), quadratic2(2, 3)
    real(c_double) :: cubic1(2, 4), cubic2(2, 4)
    integer(c_int) :: num_intersections
    logical(c_bool) :: coincident
    real(c_double) :: expected_s, expected_t
    integer(c_int) :: status
    integer :: case_id
    character(17) :: name

    case_id = 1
    name = "all_intersections"

    ! CASE 1: No intersections.
    linear1(:, 1) = 0
    linear1(:, 2) = 1
    linear2(:, 1) = 3
    linear2(:, 2) = [4.0_dp, 3.0_dp]

    call all_intersections( &
         2, linear1, 2, linear2, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         .NOT. allocated(intersections) .AND. &
         num_intersections == 0 .AND. &
         .NOT. coincident .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Intersection of two quadratic curves.
    ! NOTE: The first candidate is a specialization of [0, 0], [1/2, 1], [1, 1]
    !       onto the interval [1/4, 1] and the second candidate is
    !       specialization of [0, 1], [1/2, 1], [1, 0] onto the interval
    !       [0, 3/4]. We expect them to intersect at s = 1/3, t = 2/3, which is
    !       the point [1/2, 3/4].
    quadratic1(:, 1) = [0.25_dp, 0.4375_dp]
    quadratic1(:, 2) = [0.625_dp, 1.0_dp]
    quadratic1(:, 3) = 1
    quadratic2(:, 1) = [0.0_dp, 1.0_dp]
    quadratic2(:, 2) = [0.375_dp, 1.0_dp]
    quadratic2(:, 3) = [0.75_dp, 0.4375_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         3 * intersections(1, 1) == 1.0_dp .AND. &
         3 * intersections(2, 1) == 2.0_dp .AND. &
         .NOT. coincident .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Tangent curves which have to do a "full" Newton iteration.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [0.375_dp, 0.75_dp]
    quadratic1(:, 3) = [0.75_dp, 0.375_dp]
    quadratic2(:, 1) = [0.25_dp, 0.625_dp]
    quadratic2(:, 2) = [0.625_dp, 0.25_dp]
    quadratic2(:, 3) = 1
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, coincident, status)
    expected_s = 2.0_dp / 3.0_dp
    expected_t = 1.0_dp / 3.0_dp
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         intersections(1, 1) - expected_s == spacing(expected_s) .AND. &
         intersections(2, 1) == expected_t .AND. &
         .NOT. coincident .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Tangent curves, which cause the number of candidate pairs
    !         to become too high and will trigger a "prune" step.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [-0.5_dp, 1.5_dp]
    quadratic1(:, 3) = 1
    quadratic2(:, 1) = [-1.0_dp, 1.0_dp]
    quadratic2(:, 2) = 0.5_dp
    quadratic2(:, 3) = [0.0_dp, 2.0_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 1]) .AND. &
         num_intersections == 1 .AND. &
         all(intersections(:, 1) == 0.5_dp) .AND. &
         .NOT. coincident .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Badly scaled curves, which cause the subdivision process to
    !         take too many iterations before being "almost linear".
    ! NOTE: This is just the quadratic given by [0, 0], [4.5, 9], [9, 0] and
    !       the line given by [0, 8], [6, 0], but scaled up by a factor of
    !       2**14 == 16384. This makes the linearization error large.
    !       This definitely points out a flaw in using absolute vs. relative
    !       linearization error. However, there is an upside in not using
    !       relative error: it is faster to compute without normalization.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [73728.0_dp, 147456.0_dp]
    quadratic1(:, 3) = [147456.0_dp, 0.0_dp]
    linear1(:, 1) = [0.0_dp, 131072.0_dp]
    linear1(:, 2) = [98304.0_dp, 0.0_dp]
    call all_intersections( &
         3, quadratic1, 2, linear1, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident .AND. &
         status == Status_NO_CONVERGE)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Curves where there are duplicate intersections caused by
    !         bounding boxes touching at corners.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [0.5_dp, 1.0_dp]
    quadratic1(:, 3) = [1.0_dp, 0.0_dp]
    quadratic2(:, 1) = [0.0_dp, 0.75_dp]
    quadratic2(:, 2) = [0.5_dp, -0.25_dp]
    quadratic2(:, 3) = [1.0_dp, 0.75_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 2 .AND. &
         intersections(1, 1) == 0.25_dp .AND. &
         intersections(2, 1) == 0.25_dp .AND. &
         intersections(1, 2) == 0.75_dp .AND. &
         intersections(2, 2) == 0.75_dp .AND. &
         .NOT. coincident .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Curves that **almost** touch at their endpoints, but don't
    !         actually cross. This causes a "wiggle fail" in
    !         ``from_linearized()`` but the convex hulls are disjoint.
    cubic1(:, 1) = [-0.7838204403623438_dp, -0.25519640597397464_dp]
    cubic1(:, 2) = [-0.7894577677825452_dp, -0.24259531488131633_dp]
    cubic1(:, 3) = [-0.7946421067207265_dp, -0.22976394420044136_dp]
    cubic1(:, 4) = [-0.799367666650849_dp, -0.21671303774854855_dp]
    cubic2(:, 1) = [-0.7993236103108717_dp, -0.21683567278362156_dp]
    cubic2(:, 2) = [-0.8072986524226636_dp, -0.21898490744674426_dp]
    cubic2(:, 3) = [-0.8152736945344552_dp, -0.2211341421098668_dp]
    cubic2(:, 4) = [-0.8232487366462472_dp, -0.2232833767729893_dp]

    call all_intersections( &
         4, cubic1, 4, cubic2, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         .NOT. coincident .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Same as CASE 4, except the nodes have been re-specialized to
    !         [0, 2] and [-1, 1], respectively. This way, the first round
    !         reduces to CASE 4, so the "prune" stage happens at an odd
    !         iteration rather than en even iteration.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [-1.0_dp, 3.0_dp]
    quadratic1(:, 3) = [6.0_dp, -2.0_dp]
    quadratic2(:, 1) = [-6.0_dp, 4.0_dp]
    quadratic2(:, 2) = [1.0_dp, -1.0_dp]
    quadratic2(:, 3) = [0.0_dp, 2.0_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 2]) .AND. &
         num_intersections == 1 .AND. &
         all(intersections(:, 1) == [0.25_dp, 0.75_dp]) .AND. &
         .NOT. coincident .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Coincident curves (i.e. segments on a common curve). In this
    !         case, the second curve is just the first curve specialized to
    !         [1/2, 1]. In this case, "convex hull pruning" is of no use
    !         because at a certain point the subdivided curve segments will
    !         be identical.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [0.5_dp, 0.25_dp]
    quadratic1(:, 3) = [1.0_dp, 0.0_dp]
    quadratic2(:, 1) = [0.5_dp, 0.125_dp]
    quadratic2(:, 2) = [0.75_dp, 0.125_dp]
    quadratic2(:, 3) = [1.0_dp, 0.0_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 16]) .AND. &
         num_intersections == 2 .AND. &
         all(intersections(:, 1) == [0.5_dp, 0.0_dp]) .AND. &
         all(intersections(:, 2) == 1) .AND. &
         coincident .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: Curve are tangent and have the same curvature, but
    !          there are more than ``MAX_CANDIDATES`` candidates
    !          that can't be pruned.
    quadratic1(:, 1) = [12.0_dp, 4.0_dp]
    quadratic1(:, 2) = [-12.0_dp, -8.0_dp]
    quadratic1(:, 3) = [0.0_dp, 16.0_dp]
    quadratic2(:, 1) = [6.0_dp, 1.0_dp]
    quadratic2(:, 2) = [-6.0_dp, -2.0_dp]
    quadratic2(:, 3) = [0.0_dp, 4.0_dp]
    call all_intersections( &
         3, quadratic1, 3, quadratic2, intersections, &
         num_intersections, coincident, status)
    case_success = ( &
         allocated(intersections) .AND. &
         all(shape(intersections) == [2, 16]) .AND. &
         num_intersections == 0 .AND. &
         .NOT. coincident .AND. &
         status == 74)
    call print_status(name, case_id, case_success, success)

  end subroutine test_all_intersections

  subroutine test_all_intersections_abi(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: linear1(2, 2)
    real(c_double) :: quadratic1(2, 3), quadratic2(2, 3)
    real(c_double) :: cubic1(2, 4)
    integer(c_int) :: num_intersections
    logical(c_bool) :: coincident
    real(c_double) :: intersections1(2, 2), intersections2(2, 3)
    integer(c_int) :: status
    integer :: case_id
    character(21) :: name

    case_id = 1
    name = "all_intersections_abi"

    ! CASE 1: **Other** failure (tangent curves with the same curvature).
    quadratic1(:, 1) = [12.0_dp, 4.0_dp]
    quadratic1(:, 2) = -4
    quadratic1(:, 3) = [-4.0_dp, 4.0_dp]
    quadratic2(:, 1) = [6.0_dp, 1.0_dp]
    quadratic2(:, 2) = [-2.0_dp, -1.0_dp]
    quadratic2(:, 3) = [-2.0_dp, 1.0_dp]

    call all_intersections_abi( &
         3, quadratic1, 3, quadratic2, 2, intersections1, &
         num_intersections, coincident, status)
    case_success = ( &
         num_intersections == 1 .AND. &
         .NOT. coincident .AND. &
         status == Status_BAD_MULTIPLICITY)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``intersections`` is not large enough.
    linear1(:, 1) = [-3.0_dp, 0.0_dp]
    linear1(:, 2) = [5.0_dp, 0.0_dp]
    cubic1(:, 1) = [-7.0_dp, -9.0_dp]
    cubic1(:, 2) = [9.0_dp, 13.0_dp]
    cubic1(:, 3) = [-7.0_dp, -13.0_dp]
    cubic1(:, 4) = [9.0_dp, 9.0_dp]

    call all_intersections_abi( &
         2, linear1, 4, cubic1, 2, intersections1, &
         num_intersections, coincident, status)
    case_success = ( &
         num_intersections == 3 .AND. &
         .NOT. coincident .AND. &
         status == Status_INSUFFICIENT_SPACE)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Just case 7, but with large enough ``intersections``.
    call all_intersections_abi( &
         2, linear1, 4, cubic1, 3, intersections2, &
         num_intersections, coincident, status)
    case_success = ( &
         num_intersections == 3 .AND. &
         intersections2(1, 1) == 0.5_dp .AND. &
         intersections2(2, 1) == 0.5_dp .AND. &
         intersections2(1, 2) == 0.375_dp .AND. &
         intersections2(2, 2) == 0.25_dp .AND. &
         intersections2(1, 3) == 0.625_dp .AND. &
         intersections2(2, 3) == 0.75_dp .AND. &
         .NOT. coincident .AND. &
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

end module test_curve_intersection
