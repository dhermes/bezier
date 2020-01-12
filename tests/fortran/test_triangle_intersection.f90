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
  use status, only: &
       Status_SUCCESS, Status_BAD_MULTIPLICITY, Status_NO_CONVERGE, &
       Status_INSUFFICIENT_SPACE, Status_SAME_CURVATURE, Status_BAD_INTERIOR, &
       Status_EDGE_END, Status_UNKNOWN
  use curve, only: CurveData, LOCATE_MISS
  use surface_intersection, only: &
       Intersection, CurvedPolygonSegment, IntersectionClassification_UNSET, &
       IntersectionClassification_FIRST, IntersectionClassification_SECOND, &
       IntersectionClassification_OPPOSED, &
       IntersectionClassification_TANGENT_FIRST, &
       IntersectionClassification_TANGENT_SECOND, &
       IntersectionClassification_IGNORED_CORNER, &
       IntersectionClassification_TANGENT_BOTH, &
       IntersectionClassification_COINCIDENT, &
       IntersectionClassification_COINCIDENT_UNUSED, &
       SurfaceContained_NEITHER, &
       SurfaceContained_FIRST, SurfaceContained_SECOND, newton_refine, &
       locate_point, classify_intersection, update_edge_end_unused, &
       find_corner_unused, add_st_vals, should_keep, &
       surfaces_intersection_points, is_first, is_second, &
       get_next, to_front, add_segment, &
       interior_combine, surfaces_intersect, surfaces_intersect_abi
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private &
       test_newton_refine, test_locate_point, test_classify_intersection, &
       test_update_edge_end_unused, test_find_corner_unused, &
       test_add_st_vals, test_should_keep, test_surfaces_intersection_points, &
       test_is_first, test_is_second, &
       intersection_check, intersection_equal, test_get_next, test_to_front, &
       test_add_segment, test_interior_combine, segment_check, &
       test_surfaces_intersect, test_surfaces_intersect_abi
  public surface_intersection_all_tests

contains

  subroutine surface_intersection_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_newton_refine(success)
    call test_locate_point(success)
    call test_classify_intersection(success)
    call test_update_edge_end_unused(success)
    call test_find_corner_unused(success)
    call test_add_st_vals(success)
    call test_should_keep(success)
    call test_surfaces_intersection_points(success)
    call test_is_first(success)
    call test_is_second(success)
    call test_get_next(success)
    call test_to_front(success)
    call test_add_segment(success)
    call test_interior_combine(success)
    call test_surfaces_intersect(success)
    call test_surfaces_intersect_abi(success)

  end subroutine surface_intersection_all_tests

  subroutine test_newton_refine(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(23) :: name
    real(c_double) :: nodes(2, 6)
    real(c_double) :: updated_s, updated_t

    case_id = 1
    name = "newton_refine (Surface)"

    ! CASE 1: Quadratic surface.
    ! NOTE: This surface is given by
    !       [(4 s - t^2) / 4, (4 s^2 + 4 s t - t^2 - 4 s + 8 t) / 8]
    nodes(:, 1) = 0
    nodes(:, 2) = [0.5_dp, -0.25_dp]
    nodes(:, 3) = [1.0_dp, 0.0_dp]
    nodes(:, 4) = [0.0_dp, 0.5_dp]
    nodes(:, 5) = 0.5_dp
    nodes(:, 6) = [-0.25_dp, 0.875_dp]
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
    nodes(:, 1) = 0
    nodes(:, 2) = [0.5_dp, 0.0_dp]
    nodes(:, 3) = [1.0_dp, 0.0_dp]
    nodes(:, 4) = [0.0_dp, 0.5_dp]
    nodes(:, 5) = 0.5_dp
    nodes(:, 6) = [0.0_dp, 1.0_dp]
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
    real(c_double) :: x_val, y_val, s_val, t_val
    real(c_double) :: expected_s, expected_t
    character(23) :: name

    case_id = 1
    name = "newton_refine (Surface)"

    ! CASE 1: B(s, t) = (s, t), match.
    allocate(nodes(2, 3))
    nodes(:, 1) = 0
    nodes(:, 2) = [1.0_dp, 0.0_dp]
    nodes(:, 3) = [0.0_dp, 1.0_dp]
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
    allocate(nodes(2, 6))
    nodes(:, 1) = 0
    nodes(:, 2) = [1.0_dp, 0.0_dp]
    nodes(:, 3) = [2.0_dp, 0.0_dp]
    nodes(:, 4) = [2.0_dp, 1.0_dp]
    nodes(:, 5) = 2
    nodes(:, 6) = [0.0_dp, 2.0_dp]
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

  subroutine test_classify_intersection(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(21) :: name
    type(Intersection) :: intersection_
    type(CurveData) :: edges_first(3), edges_second(3)
    real(c_double) :: nodes1(2, 2), nodes2(2, 3)
    integer(c_int) :: enum_, status

    case_id = 1
    name = "classify_intersection"

    ! CASE 1: Intersection is on "the end" of the first edge.
    intersection_%s = 1.0_dp
    intersection_%t = 0.5_dp
    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = (status == Status_EDGE_END)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Intersection is on "the end" of the second edge.
    intersection_%s = 0.5_dp
    intersection_%t = 1.0_dp
    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = (status == Status_EDGE_END)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Intersection is on interior, "SECOND".
    nodes1(:, 1) = 0
    nodes1(:, 2) = 1
    edges_first(1)%nodes = nodes1

    nodes1(:, 1) = [0.25_dp, 0.0_dp]
    nodes1(:, 2) = [0.75_dp, 1.0_dp]
    edges_second(1)%nodes = nodes1

    intersection_ = Intersection(0.5_dp, 0.5_dp, 1, 1)

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Intersection is on interior, "FIRST".
    call classify_intersection( &
         edges_second, edges_first, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Intersection is tangent, "FIRST".
    nodes2(:, 1) = 0
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [3.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes2

    nodes2(:, 1) = [1.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [2.0_dp, 0.0_dp]
    edges_second(1)%nodes = nodes2

    intersection_ = Intersection(0.5_dp, 0.5_dp, 1, 1)

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Intersection is tangent, "SECOND".
    call classify_intersection( &
         edges_second, edges_first, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Intersection is tangent, but only to machine precision. The
    !         cross product of the tangent vectors is ``-2^{-53}``, rather
    !         than ``0.0``.
    nodes2(:, 1) = 0
    nodes2(:, 2) = [0.5_dp, -1.0_dp / 5.0_dp]
    nodes2(:, 3) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes2

    nodes2(:, 1) = [0.5_dp, -0.125_dp]
    nodes2(:, 2) = [2.0625_dp, 0.1875_dp]
    nodes2(:, 3) = [3.625_dp, 0.5_dp]
    edges_second(1)%nodes = nodes2

    intersection_ = Intersection(0.75_dp, 2.0_dp / 25.0_dp - 0.5_dp**56, 1, 1)
    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Intersection is an ignored corner.
    ! NOTE: The curves in ``edges_first`` are the edges of the
    !       "unit simplex" ...
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes1
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = [0.0_dp, 1.0_dp]
    edges_first(2)%nodes = nodes1
    nodes1(:, 1) = [0.0_dp, 1.0_dp]
    nodes1(:, 2) = 0
    edges_first(3)%nodes = nodes1
    ! ... and those in ``edges_second`` are the edges in the surface given by
    !           [ 0,  0]
    !           [-1,  0]
    !           [ 0, -1]
    nodes1(:, 1) = [-1.0_dp, 0.0_dp]
    nodes1(:, 2) = [0.0_dp, -1.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(:, 1) = [0.0_dp, -1.0_dp]
    nodes1(:, 2) = 0
    edges_second(2)%nodes = nodes1
    nodes1(:, 1) = 0
    nodes1(:, 2) = [-1.0_dp, 0.0_dp]
    edges_second(3)%nodes = nodes1

    intersection_ = Intersection(0.0_dp, 0.0_dp, 1, 3)

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_IGNORED_CORNER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Intersection is a corner.
    nodes1(:, 1) = [0.0_dp, 0.5_dp]
    nodes1(:, 2) = 1
    edges_first(1)%nodes = nodes1
    nodes1(:, 1) = 1
    nodes1(:, 2) = 0
    edges_first(2)%nodes = nodes1
    deallocate(edges_first(3)%nodes)  ! Unset.

    deallocate(edges_second(1)%nodes)  ! Unset.
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = [1.0_dp, 2.0_dp]
    edges_second(2)%nodes = nodes1
    deallocate(edges_second(3)%nodes)  ! Unset.

    intersection_ = Intersection(0.0_dp, 0.5_dp, 2, 2)

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: Intersection is tangent, use curvature of "FIRST".
    deallocate(edges_first(1)%nodes)  ! Unset.
    nodes2(:, 1) = [2.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [1.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(:, 1) = [3.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = 0
    edges_second(2)%nodes = nodes2

    intersection_ = Intersection(0.5_dp, 0.5_dp, 2, 2)

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 11: Intersection is tangent, use curvature of "SECOND".
    nodes2(:, 1) = [1.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [2.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(:, 1) = 0
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [3.0_dp, 0.0_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 12: Intersection is tangent, same curvature and direction.
    nodes2(:, 1) = [1.0_dp, 0.25_dp]
    nodes2(:, 2) = [-0.5_dp, -0.25_dp]
    nodes2(:, 3) = [0.0_dp, 0.25_dp]
    edges_first(2)%nodes = nodes2

    nodes2(:, 1) = [0.75_dp, 0.25_dp]
    nodes2(:, 2) = -0.25_dp
    nodes2(:, 3) = [-0.25_dp, 0.25_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = (status == Status_SAME_CURVATURE)
    call print_status(name, case_id, case_success, success)

    ! CASE 13: Intersection is tangent, same curvature, opposite direction.
    nodes2(:, 1) = [0.0_dp, 0.25_dp]
    nodes2(:, 2) = [-0.5_dp, -0.25_dp]
    nodes2(:, 3) = [1.0_dp, 0.25_dp]
    edges_first(2)%nodes = nodes2

    nodes2(:, 1) = [0.75_dp, 0.25_dp]
    nodes2(:, 2) = -0.25_dp
    nodes2(:, 3) = [-0.25_dp, 0.25_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = (status == Status_SAME_CURVATURE)
    call print_status(name, case_id, case_success, success)

    ! CASE 14: Intersection is tangent, opposite direction, curvatures have
    !          same sign and there is no overlap.
    nodes2(:, 1) = [2.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [1.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(:, 1) = 1
    nodes2(:, 2) = [1.5_dp, 0.0_dp]
    nodes2(:, 3) = [2.0_dp, 1.0_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_OPPOSED .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 15: Intersection is tangent, opposite direction, curvatures have
    !          same sign and there **is** overlap.
    nodes2(:, 1) = [1.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [2.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(:, 1) = [2.0_dp, 1.0_dp]
    nodes2(:, 2) = [1.5_dp, 0.0_dp]
    nodes2(:, 3) = 1
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_BOTH .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 16: Intersection is tangent, opposite direction, curvatures have
    !          opposite sign and there is no overlap.
    nodes2(:, 1) = [2.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [1.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(:, 1) = 0
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [3.0_dp, 0.0_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_OPPOSED .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 17: Intersection is tangent, opposite direction, curvatures have
    !          opposite sign and there **is** overlap.
    nodes2(:, 1) = [1.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = [2.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(:, 1) = [3.0_dp, 0.0_dp]
    nodes2(:, 2) = [1.5_dp, 1.0_dp]
    nodes2(:, 3) = 0
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_BOTH .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 18: Intersection at corner, but the corner "staddles" an edge.
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes1

    nodes1(:, 1) = 1
    nodes1(:, 2) = [0.5_dp, 0.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(:, 1) = [0.5_dp, 0.0_dp]
    nodes1(:, 2) = [1.5_dp, -1.0_dp]
    edges_second(2)%nodes = nodes1

    intersection_ = Intersection(0.5_dp, 0.0_dp, 2, 2)

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 19: Intersection at corner of **both** edges, the corners just
    !          "kiss", hence are ignored.
    nodes1(:, 1) = [0.5_dp, 1.0_dp]
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes1
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = [1.5_dp, 0.25_dp]
    edges_first(2)%nodes = nodes1

    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = [0.0_dp, 1.0_dp]
    edges_second(2)%nodes = nodes1

    intersection_ = Intersection(0.0_dp, 0.0_dp, 2, 2)

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_IGNORED_CORNER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 20: Intersection at corner of **both** edges, "FIRST" is interior.
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes1
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = [0.0_dp, 1.0_dp]
    edges_first(2)%nodes = nodes1

    nodes1(:, 1) = [0.5_dp, 0.25_dp]
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = 1
    edges_second(2)%nodes = nodes1

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 21: Intersection at corner of **both** edges, "SECOND" is interior.
    nodes1(:, 1) = [0.5_dp, 0.25_dp]
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes1
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = 1
    edges_first(2)%nodes = nodes1

    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(:, 1) = [1.0_dp, 0.0_dp]
    nodes1(:, 2) = [0.0_dp, 1.0_dp]
    edges_second(2)%nodes = nodes1

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 22: Intersection at corner of **both** edges, one corner is
    !          completely contained in the area of the other corner.
    nodes1(:, 1) = [0.5_dp, 1.0_dp]
    nodes1(:, 2) = 0
    edges_first(1)%nodes = nodes1
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.5_dp]
    edges_first(2)%nodes = nodes1

    nodes1(:, 1) = [0.0_dp, 1.0_dp]
    nodes1(:, 2) = 0
    edges_second(1)%nodes = nodes1
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    edges_second(2)%nodes = nodes1

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_classify_intersection

  subroutine test_update_edge_end_unused(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(22) :: name
    type(Intersection) :: intersections(1)
    integer(c_int) :: num_intersections, unused

    case_id = 1
    name = "update_edge_end_unused"

    unused = IntersectionClassification_COINCIDENT_UNUSED

    ! CASE 1: ``s`` is the corner, no matching intersection.
    num_intersections = 0
    call update_edge_end_unused( &
         1.0_dp, 1, 0.5_dp, 1, intersections, num_intersections)
    case_success = ( &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.0_dp, 0.5_dp, 2, 1, unused))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``s`` is the corner, matching intersection.
    num_intersections = 1
    intersections(1)%s = 0.0_dp
    intersections(1)%t = 0.25_dp
    intersections(1)%index_first = 3
    intersections(1)%index_second = 2
    call update_edge_end_unused( &
         1.0_dp, 2, 0.25_dp, 2, intersections, num_intersections)
    case_success = ( &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.0_dp, 0.25_dp, 3, 2, unused))
    call print_status(name, case_id, case_success, success)

    ! CASE 3: ``t`` is the corner, no matching intersection.
    num_intersections = 0
    call update_edge_end_unused( &
         0.5_dp, 1, 1.0_dp, 1, intersections, num_intersections)
    case_success = ( &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.5_dp, 0.0_dp, 1, 2, unused))
    call print_status(name, case_id, case_success, success)

    ! CASE 4: ``t`` is the corner, matching intersection.
    num_intersections = 1
    intersections(1)%s = 0.25_dp
    intersections(1)%t = 0.0_dp
    intersections(1)%index_first = 2
    intersections(1)%index_second = 3
    call update_edge_end_unused( &
         0.25_dp, 2, 1.0_dp, 2, intersections, num_intersections)
    case_success = ( &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.25_dp, 0.0_dp, 2, 3, unused))
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Both ``s`` and ``t`` are corners, no matching intersection.
    num_intersections = 0
    call update_edge_end_unused( &
         1.0_dp, 3, 1.0_dp, 3, intersections, num_intersections)
    case_success = ( &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.0_dp, 0.0_dp, 1, 1, unused))
    call print_status(name, case_id, case_success, success)

  end subroutine test_update_edge_end_unused

  subroutine test_find_corner_unused(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(18) :: name
    type(Intersection) :: intersections(1)
    logical(c_bool) :: found

    case_id = 1
    name = "find_corner_unused"

    ! CASE 1: ``s == 0``, not found.
    intersections(1)%index_first = 2
    call find_corner_unused( &
         0.0_dp, 1, 1, intersections, 1, found)
    case_success = (.NOT. found)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``s == 0``, found.
    intersections(1)%s = 0.0_dp
    intersections(1)%t = 0.5_dp
    intersections(1)%index_first = 2
    intersections(1)%index_second = 3
    intersections(1)%interior_curve = ( &
         IntersectionClassification_COINCIDENT_UNUSED)
    call find_corner_unused( &
         0.0_dp, 2, 3, intersections, 1, found)
    case_success = found
    call print_status(name, case_id, case_success, success)

    ! CASE 3: ``t == 0`` (i.e. ``s /= 0``), not found.
    intersections(1)%index_second = 2
    call find_corner_unused( &
         0.5_dp, 1, 1, intersections, 1, found)
    case_success = (.NOT. found)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: ``t == 0`` (i.e. ``s /= 0``), found.
    intersections(1)%s = 0.5_dp
    intersections(1)%t = 0.0_dp
    intersections(1)%index_first = 2
    intersections(1)%index_second = 3
    intersections(1)%interior_curve = ( &
         IntersectionClassification_COINCIDENT_UNUSED)
    call find_corner_unused( &
         0.5_dp, 2, 3, intersections, 1, found)
    case_success = found
    call print_status(name, case_id, case_success, success)

  end subroutine test_find_corner_unused

  subroutine test_add_st_vals(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(11) :: name
    real(c_double), allocatable :: st_vals(:, :)
    type(Intersection), allocatable :: intersections(:)
    integer(c_int) :: num_intersections
    type(CurveData) :: edges_first(3), edges_second(3)
    integer(c_int) :: first, second, coincident, known_enum
    integer(c_int) :: status

    case_id = 1
    name = "add_st_vals"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND
    coincident = IntersectionClassification_COINCIDENT
    known_enum = IntersectionClassification_UNSET

    ! CASE 1: ``intersections`` must be allocated.
    num_intersections = 0
    allocate(st_vals(2, 1))
    st_vals(:, 1) = [0.25_dp, 0.75_dp]
    case_success = .NOT. allocated(intersections)
    ! Set the relevant edges as lines.
    allocate(edges_first(2)%nodes(2, 2))
    allocate(edges_second(3)%nodes(2, 2))
    edges_first(2)%nodes(:, 1) = 0
    edges_first(2)%nodes(:, 2) = 1
    edges_second(3)%nodes(:, 1) = [1.0_dp, -0.5_dp]
    edges_second(3)%nodes(:, 2) = [0.0_dp, 0.5_dp]
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, known_enum, &
         2, 3, intersections, num_intersections, status)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         size(intersections) == 1 .AND. &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.25_dp, 0.75_dp, 2, 3, second) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``intersections`` must be re-allocated.
    st_vals(:, 1) = [0.0_dp, 0.5_dp]
    ! Set the relevant edges as lines.
    allocate(edges_first(3)%nodes(2, 2))
    allocate(edges_second(1)%nodes(2, 2))
    edges_first(3)%nodes(:, 1) = 0
    edges_first(3)%nodes(:, 2) = 1
    edges_second(1)%nodes(:, 1) = [1.0_dp, -1.0_dp]
    edges_second(1)%nodes(:, 2) = [-1.0_dp, 1.0_dp]
    case_success = ( &
         allocated(intersections) .AND. &
         num_intersections == 1)
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, known_enum, &
         3, 1, intersections, num_intersections, status)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         size(intersections) == 2 .AND. &
         num_intersections == 2 .AND. &
         intersection_check( &
         intersections(2), 0.0_dp, 0.5_dp, 3, 1, second) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: ``intersections`` is large enough.
    num_intersections = 0
    st_vals(:, 1) = [0.875_dp, 0.125_dp]
    ! Set the relevant edges as lines.
    allocate(edges_second(2)%nodes(2, 2))
    edges_first(2)%nodes(:, 1) = 0
    edges_first(2)%nodes(:, 2) = 1
    edges_second(2)%nodes(:, 1) = [0.75_dp, 0.625_dp]
    edges_second(2)%nodes(:, 2) = [1.75_dp, 2.625_dp]
    case_success = allocated(intersections)
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, known_enum, &
         2, 2, intersections, num_intersections, status)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         size(intersections) == 2 .AND. &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.875_dp, 0.125_dp, 2, 2, second) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Intersection causes error.
    st_vals(:, 1) = 0.5_dp
    deallocate(edges_first(2)%nodes)
    allocate(edges_first(2)%nodes(2, 3))
    edges_first(2)%nodes(:, 1) = [1.0_dp, 0.25_dp]
    edges_first(2)%nodes(:, 2) = [-0.5_dp, -0.25_dp]
    edges_first(2)%nodes(:, 3) = [0.0_dp, 0.25_dp]
    deallocate(edges_second(2)%nodes)
    allocate(edges_second(2)%nodes(2, 3))
    edges_second(2)%nodes(:, 1) = [0.75_dp, 0.25_dp]
    edges_second(2)%nodes(:, 2) = -0.25_dp
    edges_second(2)%nodes(:, 3) = [-0.25_dp, 0.25_dp]
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, known_enum, &
         2, 2, intersections, num_intersections, status)
    case_success = (status == Status_SAME_CURVATURE)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: ``intersections`` has edge ends, gets "over-allocated".
    deallocate(intersections)
    deallocate(st_vals)
    allocate(st_vals(2, 2))
    num_intersections = 0
    st_vals(:, 1) = [1.0_dp, 0.125_dp]
    st_vals(:, 2) = 0.5_dp
    ! Set the relevant edges as lines.
    allocate(edges_first(1)%nodes(2, 2))
    edges_first(1)%nodes(:, 1) = 0
    edges_first(1)%nodes(:, 2) = 1
    edges_second(1)%nodes(:, 1) = [0.0_dp, 1.0_dp]
    edges_second(1)%nodes(:, 2) = [1.0_dp, 0.0_dp]
    call add_st_vals( &
         edges_first, edges_second, 2, st_vals, known_enum, &
         1, 1, intersections, num_intersections, status)
    case_success = ( &
         allocated(intersections) .AND. &
         size(intersections) == 2 .AND. &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.5_dp, 0.5_dp, 1, 1, first) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Intersection is coincident.
    known_enum = IntersectionClassification_COINCIDENT
    num_intersections = 0
    st_vals(:, 1) = [0.5_dp, 0.5_dp]
    ! Set the relevant edges as lines.
    edges_first(1)%nodes(:, 1) = 0
    edges_first(1)%nodes(:, 2) = 1
    edges_second(1)%nodes(:, 1) = [1.0_dp, 0.0_dp]
    edges_second(1)%nodes(:, 2) = [0.0_dp, 1.0_dp]
    case_success = ( &
         allocated(intersections) .AND. &
         size(intersections) == 2)
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, known_enum, &
         1, 1, intersections, num_intersections, status)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         size(intersections) == 2 .AND. &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.5_dp, 0.5_dp, 1, 1, coincident) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_st_vals

  subroutine test_should_keep(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(11) :: name
    type(Intersection) :: intersection_
    logical(c_bool) :: predicate

    case_id = 1
    name = "should_keep"

    ! CASE 1: Classified as an "acceptable" value.
    intersection_%interior_curve = IntersectionClassification_FIRST
    predicate = should_keep(intersection_)
    case_success = predicate
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Classified as tangent, in a corner.
    intersection_%s = 0.0_dp
    intersection_%t = 0.5_dp
    intersection_%interior_curve = IntersectionClassification_TANGENT_SECOND
    predicate = should_keep(intersection_)
    case_success = predicate
    call print_status(name, case_id, case_success, success)

    ! CASE 3: In a corner, not classified as tangent.
    intersection_%s = 0.0_dp
    intersection_%t = 0.5_dp
    intersection_%interior_curve = IntersectionClassification_IGNORED_CORNER
    predicate = should_keep(intersection_)
    case_success = .NOT. predicate
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Classified as tangent, not in a corner.
    intersection_%s = 0.25_dp
    intersection_%t = 0.25_dp
    intersection_%interior_curve = IntersectionClassification_TANGENT_FIRST
    predicate = should_keep(intersection_)
    case_success = .NOT. predicate
    call print_status(name, case_id, case_success, success)

    ! CASE 5: "Always unused" classification type.
    intersection_%interior_curve = IntersectionClassification_OPPOSED
    predicate = should_keep(intersection_)
    case_success = .NOT. predicate
    call print_status(name, case_id, case_success, success)

  end subroutine test_should_keep

  subroutine test_surfaces_intersection_points(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(28) :: name
    real(c_double) :: linear1(2, 3), linear2(2, 3)
    real(c_double) :: quadratic1(2, 6), quadratic2(2, 6)
    type(Intersection), allocatable :: intersections(:)
    integer(c_int) :: num_intersections
    integer(c_int) :: all_types, status
    integer(c_int) :: first, second, coincident
    real(c_double) :: expected_s, expected_t

    case_id = 1
    name = "surfaces_intersection_points"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND
    coincident = IntersectionClassification_COINCIDENT

    ! CASE 1: Overlapping triangles (i.e. degree 1 surfaces).
    linear1(:, 1) = 0
    linear1(:, 2) = [2.0_dp, 0.0_dp]
    linear1(:, 3) = [1.0_dp, 2.0_dp]
    linear2(:, 1) = [0.0_dp, 1.0_dp]
    linear2(:, 2) = [2.0_dp, 1.0_dp]
    linear2(:, 3) = [1.0_dp, -1.0_dp]
    call surfaces_intersection_points( &
         3, linear1, 1, 3, linear2, 1, &
         intersections, num_intersections, all_types, status)
    case_success = ( &
         allocated(intersections) .AND. &
         size(intersections) == 6 .AND. &
         num_intersections == 6 .AND. &
         intersection_check( &
         intersections(1), 0.75_dp, 0.5_dp, 1, 2, first) .AND. &
         intersection_check( &
         intersections(2), 0.25_dp, 0.5_dp, 1, 3, second) .AND. &
         intersection_check( &
         intersections(3), 0.5_dp, 0.75_dp, 2, 1, first) .AND. &
         intersection_check( &
         intersections(4), 0.25_dp, 0.25_dp, 2, 2, second) .AND. &
         intersection_check( &
         intersections(5), 0.5_dp, 0.25_dp, 3, 1, second) .AND. &
         intersection_check( &
         intersections(6), 0.75_dp, 0.75_dp, 3, 3, first) .AND. &
         all_types == 2**first + 2**second .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Tangent intersection, which requires a "full" Newton iteration.
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [6.0_dp, 12.0_dp]
    quadratic1(:, 3) = [12.0_dp, 6.0_dp]
    quadratic1(:, 4) = [4.0_dp, -8.0_dp]
    quadratic1(:, 5) = [10.0_dp, -5.0_dp]
    quadratic1(:, 6) = [8.0_dp, -16.0_dp]
    quadratic2(:, 1) = [4.0_dp, 10.0_dp]
    quadratic2(:, 2) = [10.0_dp, 4.0_dp]
    quadratic2(:, 3) = [16.0_dp, 16.0_dp]
    quadratic2(:, 4) = [6.0_dp, 21.0_dp]
    quadratic2(:, 5) = [12.0_dp, 24.0_dp]
    quadratic2(:, 6) = [8.0_dp, 32.0_dp]
    call surfaces_intersection_points( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         intersections, num_intersections, all_types, status)
    expected_s = 2.0_dp / 3.0_dp
    expected_t = 1.0_dp / 3.0_dp
    ! Account for an error of 1 ULP in computed ``s``.
    expected_s = expected_s + spacing(expected_s)
    case_success = ( &
         allocated(intersections) .AND. &
         size(intersections) == 6 .AND. &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), expected_s, expected_t, 1, 1, second) .AND. &
         all_types == 2**second .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Tangent intersection with the same curvature.
    num_intersections = 0
    quadratic1(:, 1) = [0.0_dp, 2.0_dp]
    quadratic1(:, 2) = [-4.0_dp, -2.0_dp]
    quadratic1(:, 3) = [8.0_dp, 2.0_dp]
    quadratic1(:, 4) = [4.0_dp, 4.0_dp]
    quadratic1(:, 5) = [8.0_dp, 4.0_dp]
    quadratic1(:, 6) = [8.0_dp, 6.0_dp]
    quadratic2(:, 1) = [-2.0_dp, 2.0_dp]
    quadratic2(:, 2) = [-2.0_dp, -2.0_dp]
    quadratic2(:, 3) = [6.0_dp, 2.0_dp]
    quadratic2(:, 4) = [1.0_dp, 5.0_dp]
    quadratic2(:, 5) = [5.0_dp, 5.0_dp]
    quadratic2(:, 6) = [4.0_dp, 8.0_dp]
    call surfaces_intersection_points( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         intersections, num_intersections, all_types, status)
    case_success = ( &
         allocated(intersections) .AND. &
         size(intersections) == 6 .AND. &
         num_intersections == 0 .AND. &
         all_types == 0 .AND. &
         status == Status_BAD_MULTIPLICITY)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Coincident intersection with shared interior, which sets
    !         ``known_enum`` to ``COINCIDENT``.
    num_intersections = 0
    quadratic1(:, 1) = [1.0_dp, 0.0_dp]
    quadratic1(:, 2) = [0.5_dp, 0.25_dp]
    quadratic1(:, 3) = 0
    quadratic1(:, 4) = [0.75_dp, -0.375_dp]
    quadratic1(:, 5) = [0.25_dp, -0.375_dp]
    quadratic1(:, 6) = [0.5_dp, -0.75_dp]
    quadratic2(:, 1) = [0.75_dp, 0.09375_dp]
    quadratic2(:, 2) = [0.25_dp, 0.21875_dp]
    quadratic2(:, 3) = [-0.25_dp, -0.15625_dp]
    quadratic2(:, 4) = [0.625_dp, -0.328125_dp]
    quadratic2(:, 5) = [0.125_dp, -0.453125_dp]
    quadratic2(:, 6) = [0.5_dp, -0.75_dp]
    call surfaces_intersection_points( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         intersections, num_intersections, all_types, status)
    case_success = ( &
         allocated(intersections) .AND. &
         size(intersections) == 6 .AND. &
         num_intersections == 3 .AND. &
         intersection_check( &
         intersections(1), 0.25_dp, 0.0_dp, 1, 1, coincident) .AND. &
         intersection_check( &
         intersections(2), 0.0_dp, 0.75_dp, 2, 1, first) .AND. &
         intersection_check( &
         intersections(3), 0.0_dp, 0.0_dp, 3, 3, second) .AND. &
         all_types == (2**first + 2**second + 2**coincident) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_surfaces_intersection_points

  subroutine test_is_first(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(8) :: name
    logical(c_bool) :: predicate

    case_id = 1
    name = "is_first"

    ! CASE 1: Success.
    predicate = is_first(IntersectionClassification_FIRST)
    case_success = predicate
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Failure.
    predicate = is_first(IntersectionClassification_COINCIDENT)
    case_success = .NOT. predicate
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Success, and tangent.
    predicate = is_first(IntersectionClassification_TANGENT_FIRST)
    case_success = predicate
    call print_status(name, case_id, case_success, success)

  end subroutine test_is_first

  subroutine test_is_second(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(9) :: name
    logical(c_bool) :: predicate

    case_id = 1
    name = "is_second"

    ! CASE 1: Success.
    predicate = is_second(IntersectionClassification_SECOND)
    case_success = predicate
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Failure.
    predicate = is_second(IntersectionClassification_TANGENT_BOTH)
    case_success = .NOT. predicate
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Success, and tangent.
    predicate = is_second(IntersectionClassification_TANGENT_SECOND)
    case_success = predicate
    call print_status(name, case_id, case_success, success)

  end subroutine test_is_second

  function intersection_check( &
       intersection_, s, t, index_first, index_second, &
       interior_curve) result(same)

    type(Intersection), intent(in) :: intersection_
    real(c_double), intent(in) :: s, t
    integer(c_int), intent(in) :: index_first, index_second
    integer(c_int), intent(in) :: interior_curve
    logical(c_bool) :: same

    same = ( &
         intersection_%s == s .AND. &
         intersection_%t == t .AND. &
         intersection_%index_first == index_first .AND. &
         intersection_%index_second == index_second .AND. &
         intersection_%interior_curve == interior_curve)

  end function intersection_check

  function intersection_equal(intersection1, intersection2) result(same)

    type(Intersection), intent(in) :: intersection1, intersection2
    logical(c_bool) :: same

    same = intersection_check( &
         intersection2, &
         intersection1%s, &
         intersection1%t, &
         intersection1%index_first, &
         intersection1%index_second, &
         intersection1%interior_curve)

  end function intersection_equal

  subroutine test_get_next(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(8) :: name
    integer(c_int) :: first, second, coincident
    type(Intersection) :: intersections(5)
    type(Intersection) :: curr_node, next_node
    logical(c_bool) :: at_start
    integer(c_int) :: unused(4)
    integer(c_int) :: remaining

    case_id = 1
    name = "get_next"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND
    coincident = IntersectionClassification_COINCIDENT

    ! CASE 1: On an edge of first surface, no other intersections on that edge.
    curr_node = Intersection(0.125_dp, -1.0_dp, 2, -1, first)
    intersections(1) = Intersection(0.5_dp, -1.0_dp, 3, -1, first)
    remaining = 4
    call get_next( &
         1, intersections(:1), unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_check(next_node, 1.0_dp, -1.0_dp, 2, -1, first) .AND. &
         remaining == 4 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: On an edge of first surface, **multiple** other intersections
    !         on same edge.
    curr_node = Intersection(0.25_dp, -1.0_dp, 3, -1, first)
    ! An "acceptable" intersection that will be overtaken by the
    ! next since 0.25 < 0.5 < 0.875.
    intersections(1) = Intersection(0.875_dp, -1.0_dp, 3, -1, second)
    ! Closer to ``curr_node`` then previous.
    intersections(2) = Intersection(0.5_dp, -1.0_dp, 3, -1, first)
    ! On a different edge.
    intersections(3) = Intersection(-1.0_dp, -1.0_dp, 2, -1)
    ! Same edge, but parameter comes **before** ``curr_node``.
    intersections(4) = Intersection(0.125_dp, -1.0_dp, 3, -1, first)
    ! Past the already accepted intersection: 0.25 < 0.5 < 0.625
    intersections(5) = Intersection(0.625_dp, -1.0_dp, 3, -1, first)
    ! Populate the list of indices to be removed.
    unused = [1, 2, 3, 5]
    remaining = 4
    call get_next( &
         5, intersections, unused, remaining, &
         2, curr_node, next_node, at_start)
    case_success = ( &
         intersection_equal(next_node, intersections(2)) .AND. &
         remaining == 3 .AND. &
         all(unused == [1, 3, 5, 5]) .AND. &
         at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: On an edge of second surface, no other intersections on
    !         that edge.
    curr_node = Intersection(-1.0_dp, 0.625_dp, -1, 3, second)
    intersections(1) = Intersection(-1.0_dp, 0.5_dp, -1, 1, first)
    remaining = 2
    call get_next( &
         1, intersections(:1), unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_check(next_node, -1.0_dp, 1.0_dp, -1, 3, second) .AND. &
         remaining == 2 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: On an edge of second surface, **multiple** other intersections
    !         on same edge.
    curr_node = Intersection(-1.0_dp, 0.125_dp, -1, 2, second)
    ! An "acceptable" intersection that will be overtaken by the
    ! next since 0.125 < 0.625 < 0.75.
    intersections(1) = Intersection(-1.0_dp, 0.75_dp, -1, 2, first)
    ! On a different edge.
    intersections(2) = Intersection(-1.0_dp, -1.0_dp, -1, 3)
    ! Closer to ``curr_node`` then previous accepted.
    intersections(3) = Intersection(-1.0_dp, 0.625_dp, -1, 2, second)
    ! Same edge, but parameter comes **before** ``curr_node``.
    intersections(4) = Intersection(-1.0_dp, 0.0625_dp, -1, 2, first)
    ! Past the already accepted intersection: 0.125 < 0.625 < 0.6875
    intersections(5) = Intersection(-1.0_dp, 0.6875_dp, -1, 2, second)
    ! Populate the list of indices to be removed, missing 3.
    unused(:3) = [1, 2, 4]
    remaining = 3
    call get_next( &
         5, intersections, unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_equal(next_node, intersections(3)) .AND. &
         remaining == 3 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: On an edge of both surfaces (i.e. along a ``COINCIDENT``
    !         segment). Next node found along edge of first surface.
    curr_node = Intersection(0.25_dp, 0.0_dp, 1, 3, coincident)
    ! NOTE: ``intersections(1)`` is "nonsense". A coincident edge should
    !       only have one other intesection along it, the "nonsense"
    !       intesection is just there to hit a "just in case" code path
    !       that checks the smallest ``s`` value after ``curr_node``.
    intersections(1) = Intersection(0.875_dp, -1.0_dp, 1, 2, first)
    intersections(2) = Intersection(0.75_dp, 0.0_dp, 1, 1, second)
    unused(:2) = [2, 3]
    remaining = 2
    call get_next( &
         2, intersections(:2), unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_check(next_node, 0.75_dp, 0.0_dp, 1, 1, second) .AND. &
         remaining == 1 .AND. &
         all(unused(:2) == [3, 3]) .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: On an edge of both surfaces (i.e. along a ``COINCIDENT``
    !         segment). Next node found along edge of second surface.
    curr_node = Intersection(0.625_dp, 0.0_dp, 2, 2, coincident)
    ! NOTE: ``intersections(1)`` is "nonsense". A coincident edge should
    !       only have one other intesection along it, the "nonsense"
    !       intesection is just there to hit a "just in case" code path
    !       that checks the smallest ``t`` value after ``curr_node``.
    intersections(1) = Intersection(-1.0_dp, 0.75_dp, 1, 2, first)
    intersections(2) = Intersection(0.0_dp, 0.5_dp, 3, 2, second)
    unused(:2) = [2, 5]
    remaining = 2
    call get_next( &
         2, intersections(:2), unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_check(next_node, 0.0_dp, 0.5_dp, 3, 2, second) .AND. &
         remaining == 1 .AND. &
         all(unused(:2) == [5, 5]) .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: On an edge of both surfaces (i.e. along a ``COINCIDENT``
    !         segment). Next node found at corner of edges.
    curr_node = Intersection(0.625_dp, 0.5_dp, 1, 3, coincident)
    intersections(1) = Intersection(0.0_dp, 0.0_dp, 2, 1, first)
    unused(:2) = [1, 5]
    remaining = 2
    call get_next( &
         1, intersections(:1), unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_check( &
         next_node, 1.0_dp, 1.0_dp, 1, 3, coincident) .AND. &
         remaining == 2 .AND. &
         .NOT. at_start)

  end subroutine test_get_next

  subroutine test_to_front(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(8) :: name
    integer(c_int) :: first, second
    type(Intersection) :: intersections(1)
    type(Intersection) :: curr_node, next_node
    logical(c_bool) :: at_start
    integer(c_int) :: unused(4)
    integer(c_int) :: remaining

    case_id = 1
    name = "to_front"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND

    ! CASE 1: Unchanged node (i.e. not at the end).
    curr_node = Intersection(0.5_dp, -1.0_dp, 3, -1, first)
    remaining = 4
    call to_front( &
         0, intersections(:0), unused, remaining, &
         -1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_equal(curr_node, next_node) .AND. &
         remaining == 4 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Change the ``s`` parameter, no corresponding node.
    curr_node = Intersection(1.0_dp, -1.0_dp, 3, -1, first)
    remaining = 3
    call to_front( &
         0, intersections(:0), unused, remaining, &
         -1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_check(next_node, 0.0_dp, -1.0_dp, 1, -1, first) .AND. &
         remaining == 3 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Change the ``s`` parameter, new corner node exists.
    curr_node = Intersection(1.0_dp, -1.0_dp, 2, -1, first)
    intersections(1) = Intersection( &
         0.0_dp, 0.5_dp, 3, 1, IntersectionClassification_TANGENT_FIRST)
    unused = [1, 3, 4, 6]
    remaining = 4
    call to_front( &
         1, intersections, unused, remaining, &
         -1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_equal(next_node, intersections(1)) .AND. &
         remaining == 3 .AND. &
         all(unused == [3, 4, 6, 6]) .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Change the ``t`` parameter, no corresponding node.
    curr_node = Intersection(-1.0_dp, 1.0_dp, -1, 2, second)
    remaining = 4
    call to_front( &
         0, intersections(:0), unused, remaining, &
         -1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_check(next_node, -1.0_dp, 0.0_dp, -1, 3, second) .AND. &
         remaining == 4 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Change the ``t`` parameter, new corner node exists.
    curr_node = Intersection(-1.0_dp, 1.0_dp, -1, 1, second)
    intersections(1) = Intersection( &
         0.5_dp, 0.0_dp, 3, 2, IntersectionClassification_TANGENT_SECOND)
    remaining = 0
    call to_front( &
         1, intersections, unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         intersection_equal(next_node, intersections(1)) .AND. &
         remaining == 0 .AND. &
         at_start)
    call print_status(name, case_id, case_success, success)

  end subroutine test_to_front

  subroutine test_add_segment(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(11) :: name
    type(Intersection) :: curr_node, next_node
    type(CurvedPolygonSegment), allocatable :: segments(:)
    integer(c_int) :: count

    case_id = 1
    name = "add_segment"

    ! CASE 1: ``segments`` is not allocated; intersection is classified as
    !         SECOND.
    count = 0
    curr_node = Intersection( &
         -1.0_dp, 0.25_dp, -1, 1, IntersectionClassification_SECOND)
    next_node = Intersection(t=0.5_dp)
    case_success = .NOT. allocated(segments)
    call add_segment(curr_node, next_node, count, segments)
    case_success = ( &
         case_success .AND. &
         count == 1 .AND. &
         allocated(segments) .AND. &
         size(segments) == 1 .AND. &
         segments(1)%start == curr_node%t .AND. &
         segments(1)%end_ == next_node%t .AND. &
         segments(1)%edge_index == curr_node%index_second + 3)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``segments`` is allocated, but not large enough; intersection
    !         is classified as FIRST.
    count = 1
    curr_node = Intersection( &
         0.5_dp, -1.0_dp, 2, -1, IntersectionClassification_FIRST)
    next_node = Intersection(s=1.0_dp)
    case_success = ( &
         allocated(segments) .AND. &
         size(segments) == 1)
    call add_segment(curr_node, next_node, count, segments)
    case_success = ( &
         case_success .AND. &
         count == 2 .AND. &
         allocated(segments) .AND. &
         size(segments) == 2 .AND. &
         segments(2)%start == curr_node%s .AND. &
         segments(2)%end_ == next_node%s .AND. &
         segments(2)%edge_index == curr_node%index_first)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: ``segments`` is allocated and does not need to be resized;
    !         intersection is classified as FIRST.
    count = 1
    curr_node = Intersection( &
         0.125_dp, -1.0_dp, 3, -1, IntersectionClassification_FIRST)
    next_node = Intersection(s=0.75_dp)
    case_success = ( &
         allocated(segments) .AND. &
         size(segments) == 2)
    call add_segment(curr_node, next_node, count, segments)
    case_success = ( &
         case_success .AND. &
         count == 2 .AND. &
         allocated(segments) .AND. &
         size(segments) == 2 .AND. &
         segments(2)%start == curr_node%s .AND. &
         segments(2)%end_ == next_node%s .AND. &
         segments(2)%edge_index == curr_node%index_first)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Intersection is classified as COINCIDENT, nodes match
    !         on first edge.
    count = 0
    curr_node = Intersection( &
         0.125_dp, -1.0_dp, 3, -1, IntersectionClassification_COINCIDENT)
    next_node = Intersection(s=0.625_dp, index_first=3)
    case_success = ( &
         allocated(segments) .AND. &
         size(segments) == 2)
    call add_segment(curr_node, next_node, count, segments)
    case_success = ( &
         case_success .AND. &
         count == 1 .AND. &
         allocated(segments) .AND. &
         size(segments) == 2 .AND. &
         segments(1)%start == curr_node%s .AND. &
         segments(1)%end_ == next_node%s .AND. &
         segments(1)%edge_index == curr_node%index_first)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Intersection is classified as COINCIDENT, nodes match
    !         on second edge.
    count = 0
    curr_node = Intersection( &
         -1.0_dp, 0.25_dp, -1, 1, IntersectionClassification_COINCIDENT)
    next_node = Intersection(-1.0_dp, 1.0_dp, -2, 1)
    case_success = ( &
         allocated(segments) .AND. &
         size(segments) == 2)
    call add_segment(curr_node, next_node, count, segments)
    case_success = ( &
         case_success .AND. &
         count == 1 .AND. &
         allocated(segments) .AND. &
         size(segments) == 2 .AND. &
         segments(1)%start == curr_node%t .AND. &
         segments(1)%end_ == next_node%t .AND. &
         segments(1)%edge_index == curr_node%index_second + 3)
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_segment

  subroutine test_interior_combine(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(16) :: name
    integer(c_int) :: first, second
    type(Intersection) :: intersections(4)
    integer(c_int), allocatable :: segment_ends(:)
    type(CurvedPolygonSegment), allocatable :: segments(:)
    integer(c_int) :: num_intersected, contained, status

    case_id = 1
    name = "interior_combine"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND

    ! CASE 1: ``segments`` and ``segment_ends`` are not allocated,
    !         from 1L-6L (ID: 23).
    intersections(1) = Intersection(0.25_dp, 0.75_dp, 1, 3, second)
    intersections(2) = Intersection(0.75_dp, 0.25_dp, 3, 1, first)

    case_success = ( &
         .NOT. allocated(segments) .AND. &
         .NOT. allocated(segment_ends))
    call interior_combine( &
         2, intersections(:2), &
         num_intersected, segment_ends, segments, &
         contained, status)
    case_success = ( &
         case_success .AND. &
         num_intersected == 1 .AND. &
         allocated(segment_ends) .AND. &
         all(segment_ends == [4]) .AND. &
         allocated(segments) .AND. &
         size(segments) == 4 .AND. &
         segment_check(segments(1), 0.75_dp, 1.0_dp, 3) .AND. &
         segment_check(segments(2), 0.0_dp, 0.25_dp, 1) .AND. &
         segment_check(segments(3), 0.75_dp, 1.0_dp, 6) .AND. &
         segment_check(segments(4), 0.0_dp, 0.25_dp, 4) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``segments`` and ``segment_ends`` must be re-allocated
    !         (``segment_ends`` is re-allocated in ``finalize_segment()``);
    !         two disjoint intersections, from 13L-38Q (ID: 39).
    intersections(1) = Intersection(0.625_dp, 0.25_dp, 1, 1, first)
    intersections(2) = Intersection(0.375_dp, 0.75_dp, 1, 1, second)
    intersections(3) = Intersection(0.3125_dp, 0.25_dp, 1, 2, first)
    intersections(4) = Intersection(0.6875_dp, 0.75_dp, 1, 3, second)

    case_success = ( &
         allocated(segment_ends) .AND. &
         size(segment_ends) == 1 .AND. &
         allocated(segments) .AND. &
         size(segments) == 4)
    call interior_combine( &
         4, intersections, &
         num_intersected, segment_ends, segments, &
         contained, status)
    case_success = ( &
         case_success .AND. &
         num_intersected == 2 .AND. &
         allocated(segment_ends) .AND. &
         all(segment_ends == [3, 6]) .AND. &
         allocated(segments) .AND. &
         size(segments) == 6 .AND. &
         segment_check(segments(1), 0.75_dp, 1.0_dp, 6) .AND. &
         segment_check(segments(2), 0.0_dp, 0.25_dp, 4) .AND. &
         segment_check(segments(3), 0.625_dp, 0.6875_dp, 1) .AND. &
         segment_check(segments(4), 0.3125_dp, 0.375_dp, 1) .AND. &
         segment_check(segments(5), 0.75_dp, 1.0_dp, 4) .AND. &
         segment_check(segments(6), 0.0_dp, 0.25_dp, 5) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: ``segments`` has been allocated but is not large enough (at the
    !         outset of the subroutine ``segments`` is resized based on the
    !         size of ``intersections``); intersection terminates after node
    !         is "rotated" via ``to_front()``, from 1L-28Q (ID: 25) with
    !         rounded parameters from ``s = 0.40587727318528016`` and
    !         ``t == 0.59412272681471978``.
    intersections(1) = Intersection(0.375_dp, 0.5625_dp, 2, 2, second)
    intersections(2) = Intersection(0.1875_dp, 0.0_dp, 1, 3, first)
    ! Note that the corner is the last intersection, which will be used
    ! as ``start``.
    deallocate(segments)
    allocate(segments(1))  ! Too small.
    case_success = ( &
         allocated(segment_ends) .AND. &
         size(segment_ends) == 2 .AND. &
         allocated(segments) .AND. &
         size(segments) == 1)
    call interior_combine( &
         2, intersections(:2), &
         num_intersected, segment_ends, segments, &
         contained, status)
    case_success = ( &
         case_success .AND. &
         num_intersected == 1 .AND. &
         allocated(segment_ends) .AND. &
         size(segment_ends) == 2 .AND. &  ! Unchanged
         segment_ends(1) == 3 .AND. &
         allocated(segments) .AND. &
         size(segments) == 3 .AND. &
         segment_check(segments(1), 0.1875_dp, 1.0_dp, 1) .AND. &
         segment_check(segments(2), 0.0_dp, 0.375_dp, 2) .AND. &
         segment_check(segments(3), 0.5625_dp, 1.0_dp, 5) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: ``segments`` is larger than needed; inspired by 11Q-12Q (ID: 16)
    !         with rounded parameters from ``0.14644660940672624`` and
    !         ``0.8535533905932737``.
    intersections(1) = Intersection(0.140625_dp, 0.84375_dp, 1, 1, first)
    intersections(2) = Intersection(0.84375_dp, 0.140625_dp, 1, 1, second)
    case_success = ( &
         allocated(segment_ends) .AND. &
         size(segment_ends) == 2 .AND. &
         allocated(segments) .AND. &
         size(segments) == 3)
    call interior_combine( &
         2, intersections(:2), &
         num_intersected, segment_ends, segments, &
         contained, status)
    case_success = ( &
         case_success .AND. &
         num_intersected == 1 .AND. &
         allocated(segment_ends) .AND. &
         size(segment_ends) == 2 .AND. &  ! Unchanged
         segment_ends(1) == 2 .AND. &
         allocated(segments) .AND. &
         size(segments) == 3 .AND. &  ! Unchanged
         segment_check(segments(1), 0.140625_dp, 0.84375_dp, 4) .AND. &
         segment_check(segments(2), 0.140625_dp, 0.84375_dp, 1) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: First surface is completely contained in the other; from
    !         13Q-35Q (ID: 37).
    intersections(1) = Intersection(0.0_dp, 0.5_dp, 3, 3, first)
    call interior_combine( &
         1, intersections(:1), &
         num_intersected, segment_ends, segments, &
         contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         contained == SurfaceContained_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Second surface is completely contained in the other; from
    !         1L-9L (ID: 30).
    intersections(1) = Intersection(0.875_dp, 0.0_dp, 1, 2, second)
    intersections(2) = Intersection(0.75_dp, 0.0_dp, 2, 3, second)
    intersections(3) = Intersection(0.875_dp, 0.0_dp, 3, 1, second)
    call interior_combine( &
         3, intersections(:3), &
         num_intersected, segment_ends, segments, &
         contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         contained == SurfaceContained_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Intersection only contains "whole" edges of each surface, but
    !         is not either one of the surface; from 1L-39Q (ID: 40).
    intersections(1) = Intersection(0.0_dp, 0.0_dp, 2, 1, second)
    intersections(2) = Intersection(0.0_dp, 0.0_dp, 3, 2, first)
    call interior_combine( &
         2, intersections(:2), &
         num_intersected, segment_ends, segments, &
         contained, status)
    case_success = ( &
         num_intersected == 1 .AND. &
         allocated(segment_ends) .AND. &
         size(segment_ends) == 2 .AND. &
         segment_ends(1) == 3 .AND. &
         allocated(segments) .AND. &
         size(segments) == 3 .AND. &
         segment_check(segments(1), 0.0_dp, 1.0_dp, 3) .AND. &
         segment_check(segments(2), 0.0_dp, 1.0_dp, 1) .AND. &
         segment_check(segments(3), 0.0_dp, 1.0_dp, 4) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_interior_combine

  function segment_check(segment, start, end_, edge_index) result(same)

    type(CurvedPolygonSegment), intent(in) :: segment
    real(c_double), intent(in) :: start, end_
    integer(c_int), intent(in) :: edge_index
    logical(c_bool) :: same

    same = ( &
         segment%start == start .AND. &
         segment%end_ == end_ .AND. &
         segment%edge_index == edge_index)

  end function segment_check

  subroutine test_surfaces_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(18) :: name
    integer(c_int) :: num_intersected, contained, status
    real(c_double) :: linear1(2, 3), linear2(2, 3)
    real(c_double) :: quadratic1(2, 6), quadratic2(2, 6)
    real(c_double) :: cubic1(2, 10), cubic2(2, 10)
    integer(c_int), allocatable :: segment_ends(:)
    type(CurvedPolygonSegment), allocatable :: segments(:)
    real(c_double) :: start, end_

    case_id = 1
    name = "surfaces_intersect"

    ! CASE 1: Disjoint bounding boxes.
    linear1(:, 1) = 0
    linear1(:, 2) = [1.0_dp, 0.0_dp]
    linear1(:, 3) = [0.0_dp, 1.0_dp]
    linear2(:, 1) = 10
    linear2(:, 2) = [11.0_dp, 10.0_dp]
    linear2(:, 3) = [10.0_dp, 11.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         .NOT. allocated(segment_ends) .AND. &
         .NOT. allocated(segments) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Intersecting edges results in error.
    quadratic1(:, 1) = [0.0_dp, 2.0_dp]
    quadratic1(:, 2) = [-4.0_dp, -2.0_dp]
    quadratic1(:, 3) = [8.0_dp, 2.0_dp]
    quadratic1(:, 4) = [4.0_dp, 4.0_dp]
    quadratic1(:, 5) = [8.0_dp, 4.0_dp]
    quadratic1(:, 6) = [8.0_dp, 6.0_dp]
    quadratic2(:, 1) = [-2.0_dp, 2.0_dp]
    quadratic2(:, 2) = [-2.0_dp, -2.0_dp]
    quadratic2(:, 3) = [6.0_dp, 2.0_dp]
    quadratic2(:, 4) = [1.0_dp, 5.0_dp]
    quadratic2(:, 5) = [5.0_dp, 5.0_dp]
    quadratic2(:, 6) = [4.0_dp, 8.0_dp]

    call surfaces_intersect( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         .NOT. allocated(segments) .AND. &
         .NOT. allocated(segment_ends) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_BAD_MULTIPLICITY)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Surface 2 contained in surface 1.
    linear1(:, 1) = 0
    linear1(:, 2) = [5.0_dp, 0.0_dp]
    linear1(:, 3) = [0.0_dp, 5.0_dp]
    linear2(:, 1) = [2.0_dp, 1.0_dp]
    linear2(:, 2) = [3.0_dp, 1.0_dp]
    linear2(:, 3) = 2

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         .NOT. allocated(segment_ends) .AND. &
         .NOT. allocated(segments) .AND. &
         contained == SurfaceContained_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Surface 1 contained in surface 2 (re-uses all data from CASE 3).
    call surfaces_intersect( &
         3, linear2, 1, 3, linear1, 1, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         .NOT. allocated(segment_ends) .AND. &
         .NOT. allocated(segments) .AND. &
         contained == SurfaceContained_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Surfaces disjoint with overlapping bounding boxes (re-uses
    !         ``linear1`` from previous cases).
    linear2(:, 1) = [7.0_dp, 2.0_dp]
    linear2(:, 2) = [7.0_dp, 3.0_dp]
    linear2(:, 3) = [6.0_dp, 3.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         .NOT. allocated(segment_ends) .AND. &
         .NOT. allocated(segments) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Surfaces intersect.
    linear1(:, 1) = 0
    linear1(:, 2) = [8.0_dp, 0.0_dp]
    linear1(:, 3) = [0.0_dp, 8.0_dp]
    linear2(:, 1) = [4.0_dp, 5.0_dp]
    linear2(:, 2) = [-4.0_dp, 5.0_dp]
    linear2(:, 3) = [4.0_dp, -3.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 1 .AND. &
         allocated(segment_ends) .AND. &
         all(segment_ends == [6]) .AND. &
         allocated(segments) .AND. &
         size(segments) == 6 .AND. &
         segment_check(segments(1), 0.5_dp, 0.625_dp, 5) .AND. &
         segment_check(segments(2), 0.125_dp, 0.5_dp, 1) .AND. &
         segment_check(segments(3), 0.375_dp, 0.875_dp, 6) .AND. &
         segment_check(segments(4), 0.5_dp, 0.625_dp, 2) .AND. &
         segment_check(segments(5), 0.125_dp, 0.5_dp, 4) .AND. &
         segment_check(segments(6), 0.375_dp, 0.875_dp, 3) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: The only intersection(s) are ``OPPOSED``.
    quadratic1(:, 1) = [4.0_dp, 0.0_dp]
    quadratic1(:, 2) = [2.0_dp, 4.0_dp]
    quadratic1(:, 3) = 0
    quadratic1(:, 4) = [3.0_dp, -2.0_dp]
    quadratic1(:, 5) = [1.0_dp, -2.0_dp]
    quadratic1(:, 6) = [2.0_dp, -4.0_dp]
    quadratic2(:, 1) = [0.0_dp, 4.0_dp]
    quadratic2(:, 2) = [2.0_dp, 0.0_dp]
    quadratic2(:, 3) = 4
    quadratic2(:, 4) = [1.0_dp, 6.0_dp]
    quadratic2(:, 5) = [3.0_dp, 6.0_dp]
    quadratic2(:, 6) = [2.0_dp, 8.0_dp]

    call surfaces_intersect( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         allocated(segment_ends) .AND. &  ! Though unused, not de-allocated.
         allocated(segments) .AND. &  ! Though unused, not de-allocated.
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: The only intersection(s) are ``IGNORED_CORNER``.
    linear1(:, 1) = 0
    linear1(:, 2) = [1.0_dp, 0.0_dp]
    linear1(:, 3) = [0.0_dp, 1.0_dp]
    linear2(:, 1) = 0
    linear2(:, 2) = [-2.0_dp, 3.0_dp]
    linear2(:, 3) = [1.0_dp, -3.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         allocated(segment_ends) .AND. &  ! Though unused, not de-allocated.
         allocated(segments) .AND. &  ! Though unused, not de-allocated.
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: The only intersection(s) are ``TANGENT_FIRST``.
    quadratic1(:, 1) = [2.0_dp, 1.25_dp]
    quadratic1(:, 2) = [4.0_dp, 0.75_dp]
    quadratic1(:, 3) = [6.0_dp, 1.25_dp]
    quadratic1(:, 4) = 3
    quadratic1(:, 5) = [5.0_dp, 3.0_dp]
    quadratic1(:, 6) = [4.0_dp, 5.0_dp]
    quadratic2(:, 1) = 0
    quadratic2(:, 2) = [4.0_dp, 2.0_dp]
    quadratic2(:, 3) = [8.0_dp, 0.0_dp]
    quadratic2(:, 4) = [1.5_dp, 7.5_dp]
    quadratic2(:, 5) = [11.0_dp, 6.0_dp]
    quadratic2(:, 6) = 8

    call surfaces_intersect( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         allocated(segment_ends) .AND. &  ! Though unused, not de-allocated.
         allocated(segments) .AND. &  ! Though unused, not de-allocated.
         contained == SurfaceContained_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: The only intersection(s) are ``TANGENT_SECOND`` (re-uses
    !          **all** data from CASE 9).
    call surfaces_intersect( &
         6, quadratic2, 2, 6, quadratic1, 2, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         allocated(segment_ends) .AND. &  ! Though unused, not de-allocated.
         allocated(segments) .AND. &  ! Though unused, not de-allocated.
         contained == SurfaceContained_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 11: In a previous implementation, round-off issues in
    !          ``curve_intersection.f90::add_intersection()`` caused a
    !          "phantom" extra intersection point caused issues in
    !          ``interior_combine()``. This case is kept around to guard
    !          against regression.
    cubic1(:, 1) = [-0.519247936646441_dp, 0.008196262233806585_dp]
    cubic1(:, 2) = [-0.5206153113582636_dp, 0.01587839530234961_dp]
    cubic1(:, 3) = [-0.5220361984817866_dp, 0.023564120023660037_dp]
    cubic1(:, 4) = [-0.5234919947680754_dp, 0.03126920999203689_dp]
    cubic1(:, 5) = [-0.5306013153328831_dp, 0.010994185722894958_dp]
    cubic1(:, 6) = [-0.5320185860555878_dp, 0.018720567953255777_dp]
    cubic1(:, 7) = [-0.533471521483344_dp, 0.026443532602826243_dp]
    cubic1(:, 8) = [-0.5418564772120571_dp, 0.013842126833196683_dp]
    cubic1(:, 9) = [-0.5433055706504355_dp, 0.02160135650021924_dp]
    cubic1(:, 10) = [-0.5530139922771271_dp, 0.016726767154940626_dp]
    cubic2(:, 1) = [-0.5492475273303934_dp, 0.004806678750684627_dp]
    cubic2(:, 2) = [-0.5507539103531026_dp, 0.011215423321262663_dp]
    cubic2(:, 3) = [-0.552283211157318_dp, 0.017577411042302475_dp]
    cubic2(:, 4) = [-0.553868947686436_dp, 0.023906392050982415_dp]
    cubic2(:, 5) = [-0.5543287770635862_dp, 0.0032189095236835417_dp]
    cubic2(:, 6) = [-0.5558905319225977_dp, 0.00960140358887212_dp]
    cubic2(:, 7) = [-0.5574932750785362_dp, 0.015938552164569394_dp]
    cubic2(:, 8) = [-0.5596130517429451_dp, 0.0014977481523963628_dp]
    cubic2(:, 9) = [-0.5612419767391262_dp, 0.007849282849502192_dp]
    cubic2(:, 10) = [-0.5650788564504334_dp, -0.0003406439314109452_dp]

    call surfaces_intersect( &
         10, cubic1, 3, 10, cubic2, 3, &
         segment_ends, segments, &
         num_intersected, contained, status)
    start = 0.60937510406326056_dp
    end_ = 0.64417397312262270_dp
    case_success = ( &
         num_intersected == 1 .AND. &
         allocated(segment_ends) .AND. &
         all(segment_ends == [3]) .AND. &
         allocated(segments) .AND. &
         segment_check(segments(1), start, end_, 4) .AND. &
         segment_check(segments(2), 0.97192953004411586_dp, 1.0_dp, 2) .AND. &
         segment_check(segments(3), 0.0_dp, 0.029255079571205871_dp, 3) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 12: Two surfaces share an edge but no interior, will result in
    !          ``COINCIDENT_UNUSED``. This is 29Q-42Q (ID: 43).
    quadratic1(:, 1) = -0.5_dp
    quadratic1(:, 2) = [0.3125_dp, -0.25_dp]
    quadratic1(:, 3) = [1.25_dp, 0.0_dp]
    quadratic1(:, 4) = [-0.25_dp, 0.3125_dp]
    quadratic1(:, 5) = 0.125_dp
    quadratic1(:, 6) = [0.0_dp, 1.25_dp]
    quadratic2(:, 1) = [0.0_dp, 1.25_dp]
    quadratic2(:, 2) = 0.125_dp
    quadratic2(:, 3) = [1.25_dp, 0.0_dp]
    quadratic2(:, 4) = [0.625_dp, 1.25_dp]
    quadratic2(:, 5) = [1.25_dp, 0.625_dp]
    quadratic2(:, 6) = 1.25_dp
    call surfaces_intersect( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 13: Same as CASE 12, with arguments swapped.
    call surfaces_intersect( &
         6, quadratic2, 2, 6, quadratic1, 2, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 14: Surfaces are disjoint, but bounding boxes are not.
    linear1(:, 1) = 0
    linear1(:, 2) = [8.0_dp, 0.0_dp]
    linear1(:, 3) = [0.0_dp, 8.0_dp]
    linear2(:, 1) = [6.0_dp, 3.0_dp]
    linear2(:, 2) = 6
    linear2(:, 3) = [3.0_dp, 6.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, &
         segment_ends, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_surfaces_intersect

  subroutine test_surfaces_intersect_abi(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(22) :: name
    real(c_double) :: linear1(2, 3), linear2(2, 3)
    real(c_double) :: quadratic1(2, 6), quadratic2(2, 6)
    integer(c_int) :: num_intersected, contained, status
    integer(c_int) :: segment_ends(2)
    type(CurvedPolygonSegment) :: segments(6)

    case_id = 1
    name = "surfaces_intersect_abi"

    ! CASE 1: Intersection results in error.
    quadratic1(:, 1) = [0.0_dp, 2.0_dp]
    quadratic1(:, 2) = [-4.0_dp, -2.0_dp]
    quadratic1(:, 3) = [8.0_dp, 2.0_dp]
    quadratic1(:, 4) = [4.0_dp, 4.0_dp]
    quadratic1(:, 5) = [8.0_dp, 4.0_dp]
    quadratic1(:, 6) = [8.0_dp, 6.0_dp]
    quadratic2(:, 1) = [-2.0_dp, 2.0_dp]
    quadratic2(:, 2) = [-2.0_dp, -2.0_dp]
    quadratic2(:, 3) = [6.0_dp, 2.0_dp]
    quadratic2(:, 4) = [1.0_dp, 5.0_dp]
    quadratic2(:, 5) = [5.0_dp, 5.0_dp]
    quadratic2(:, 6) = [4.0_dp, 8.0_dp]
    call surfaces_intersect_abi( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         2, segment_ends, 6, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_BAD_MULTIPLICITY)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``segment_ends`` is not large enough; two disjoint intersections
    !         from 13L-38Q (ID: 39).
    linear1(:, 1) = [-8.0_dp, 0.0_dp]
    linear1(:, 2) = [8.0_dp, 0.0_dp]
    linear1(:, 3) = [0.0_dp, 8.0_dp]
    quadratic1(:, 1) = [4.0_dp, 3.0_dp]
    quadratic1(:, 2) = [0.0_dp, -5.0_dp]
    quadratic1(:, 3) = [-4.0_dp, 3.0_dp]
    quadratic1(:, 4) = [2.0_dp, -3.0_dp]
    quadratic1(:, 5) = [-2.0_dp, -3.0_dp]
    quadratic1(:, 6) = [0.0_dp, -9.0_dp]
    segment_ends(:2) = [-1337, -42]
    call surfaces_intersect_abi( &
         3, linear1, 1, 6, quadratic1, 2, &
         1, segment_ends, 3, segments(:3), &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 2 .AND. &
         all(segment_ends(:2) == [-1337, -42]) .AND. &  ! Unchanged.
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_INSUFFICIENT_SPACE)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: ``segments`` is not large enough; two disjoint intersections
    !         from 13L-38Q (ID: 39).
    call surfaces_intersect_abi( &
         3, linear1, 1, 6, quadratic1, 2, &
         2, segment_ends, 3, segments(:3), &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 2 .AND. &
         all(segment_ends(:2) == [3, 6]) .AND. &  ! Changed.
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_INSUFFICIENT_SPACE)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Successful intersection; two disjoint intersections
    !         from 13L-38Q (ID: 39).
    call surfaces_intersect_abi( &
         3, linear1, 1, 6, quadratic1, 2, &
         2, segment_ends, 6, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 2 .AND. &
         all(segment_ends(:2) == [3, 6]) .AND. &
         segment_check(segments(1), 0.75_dp, 1.0_dp, 6) .AND. &
         segment_check(segments(2), 0.0_dp, 0.25_dp, 4) .AND. &
         segment_check(segments(3), 0.625_dp, 0.6875_dp, 1) .AND. &
         segment_check(segments(4), 0.3125_dp, 0.375_dp, 1) .AND. &
         segment_check(segments(5), 0.75_dp, 1.0_dp, 4) .AND. &
         segment_check(segments(6), 0.0_dp, 0.25_dp, 5) .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: No intersections.
    linear1(:, 1) = 0
    linear1(:, 2) = [1.0_dp, 0.0_dp]
    linear1(:, 3) = [0.0_dp, 1.0_dp]
    linear2(:, 1) = 3
    linear2(:, 2) = [4.0_dp, 3.0_dp]
    linear2(:, 3) = [3.0_dp, 4.0_dp]
    call surfaces_intersect_abi( &
         3, linear1, 1, 3, linear2, 1, &
         2, segment_ends, 6, segments, &
         num_intersected, contained, status)
    case_success = ( &
         num_intersected == 0 .AND. &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_surfaces_intersect_abi

end module test_surface_intersection
