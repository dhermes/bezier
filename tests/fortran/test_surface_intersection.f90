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
       Status_SUCCESS, Status_PARALLEL, Status_SAME_CURVATURE, &
       Status_BAD_TANGENT, Status_EDGE_END, Status_UNKNOWN
  use curve, only: CurveData, LOCATE_MISS
  use surface_intersection, only: &
       Intersection, SegmentNode, IntersectionClassification_FIRST, &
       IntersectionClassification_SECOND, IntersectionClassification_OPPOSED, &
       IntersectionClassification_TANGENT_FIRST, &
       IntersectionClassification_TANGENT_SECOND, &
       IntersectionClassification_IGNORED_CORNER, SurfaceContained_NEITHER, &
       SurfaceContained_FIRST, SurfaceContained_SECOND, newton_refine, &
       locate_point, classify_intersection, add_st_vals, &
       surfaces_intersection_points, remove_node, get_next, to_front, &
       surfaces_intersect
  use types, only: dp
  use unit_test_helpers, only: print_status
  implicit none
  private &
       test_newton_refine, test_locate_point, test_classify_intersection, &
       test_add_st_vals, test_surfaces_intersection_points, &
       intersection_check, test_remove_node, test_get_next, test_to_front, &
       test_surfaces_intersect
  public surface_intersection_all_tests

contains

  subroutine surface_intersection_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_newton_refine(success)
    call test_locate_point(success)
    call test_classify_intersection(success)
    call test_add_st_vals(success)
    call test_surfaces_intersection_points(success)
    call test_remove_node(success)
    call test_get_next(success)
    call test_to_front(success)
    call test_surfaces_intersect(success)

  end subroutine surface_intersection_all_tests

  subroutine test_newton_refine(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(23) :: name
    real(c_double) :: nodes(6, 2)
    real(c_double) :: updated_s, updated_t

    case_id = 1
    name = "newton_refine (Surface)"

    ! CASE 1: Quadratic surface.
    ! NOTE: This surface is given by
    !       [(4 s - t^2) / 4, (4 s^2 + 4 s t - t^2 - 4 s + 8 t) / 8]
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
    real(c_double) :: x_val, y_val, s_val, t_val
    real(c_double) :: expected_s, expected_t
    character(23) :: name

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

  subroutine test_classify_intersection(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(21) :: name
    type(Intersection) :: intersection_
    type(CurveData) :: edges_first(3), edges_second(3)
    real(c_double) :: nodes1(2, 2), nodes2(3, 2)
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

    ! CASE 3: Intersection is on interior, "second".
    nodes1(1, :) = 0
    nodes1(2, :) = 1
    edges_first(1)%nodes = nodes1

    nodes1(1, :) = [0.25_dp, 0.0_dp]
    nodes1(2, :) = [0.75_dp, 1.0_dp]
    edges_second(1)%nodes = nodes1

    intersection_%s = 0.5_dp
    intersection_%index_first = 1
    intersection_%t = 0.5_dp
    intersection_%index_second = 1

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Intersection is on interior, "first".
    call classify_intersection( &
         edges_second, edges_first, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Intersection is tangent, "first".
    nodes2(1, :) = 0
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [3.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes2

    nodes2(1, :) = [1.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [2.0_dp, 0.0_dp]
    edges_second(1)%nodes = nodes2

    intersection_%s = 0.5_dp
    intersection_%index_first = 1
    intersection_%t = 0.5_dp
    intersection_%index_second = 1

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Intersection is tangent, "second".
    call classify_intersection( &
         edges_second, edges_first, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Intersection is an ignored corner.
    ! NOTE: The curves in ``edges_first`` are the edges of the
    !       "unit simplex" ...
    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes1
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [0.0_dp, 1.0_dp]
    edges_first(2)%nodes = nodes1
    nodes1(1, :) = [0.0_dp, 1.0_dp]
    nodes1(2, :) = 0
    edges_first(3)%nodes = nodes1
    ! ... and those in ``edges_second`` are the edges in the surface given by
    !           [ 0,  0]
    !           [-1,  0]
    !           [ 0, -1]
    nodes1(1, :) = [-1.0_dp, 0.0_dp]
    nodes1(2, :) = [0.0_dp, -1.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(1, :) = [0.0_dp, -1.0_dp]
    nodes1(2, :) = 0
    edges_second(2)%nodes = nodes1
    nodes1(1, :) = 0
    nodes1(2, :) = [-1.0_dp, 0.0_dp]
    edges_second(3)%nodes = nodes1

    intersection_%s = 0.0_dp
    intersection_%index_first = 1
    intersection_%t = 0.0_dp
    intersection_%index_second = 3

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_IGNORED_CORNER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Intersection is a corner.
    nodes1(1, :) = [0.0_dp, 0.5_dp]
    nodes1(2, :) = 1
    edges_first(1)%nodes = nodes1
    nodes1(1, :) = 1
    nodes1(2, :) = 0
    edges_first(2)%nodes = nodes1
    deallocate(edges_first(3)%nodes)  ! Unset.

    deallocate(edges_second(1)%nodes)  ! Unset.
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 2.0_dp]
    edges_second(2)%nodes = nodes1
    deallocate(edges_second(3)%nodes)  ! Unset.

    intersection_%s = 0.0_dp
    intersection_%index_first = 2
    intersection_%t = 0.5_dp
    intersection_%index_second = 2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Intersection is tangent, use curvature of "first".
    deallocate(edges_first(1)%nodes)  ! Unset.
    nodes2(1, :) = [2.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [1.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(1, :) = [3.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [0.0_dp, 0.0_dp]
    edges_second(2)%nodes = nodes2

    intersection_%s = 0.5_dp
    intersection_%index_first = 2
    intersection_%t = 0.5_dp
    intersection_%index_second = 2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: Intersection is tangent, use curvature of "second".
    nodes2(1, :) = [1.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [2.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(1, :) = [0.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [3.0_dp, 0.0_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_TANGENT_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 11: Intersection is tangent, same curvature and direction.
    nodes2(1, :) = [1.0_dp, 0.25_dp]
    nodes2(2, :) = [-0.5_dp, -0.25_dp]
    nodes2(3, :) = [0.0_dp, 0.25_dp]
    edges_first(2)%nodes = nodes2

    nodes2(1, :) = [0.75_dp, 0.25_dp]
    nodes2(2, :) = [-0.25_dp, -0.25_dp]
    nodes2(3, :) = [-0.25_dp, 0.25_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = (status == Status_SAME_CURVATURE)
    call print_status(name, case_id, case_success, success)

    ! CASE 12: Intersection is tangent, same curvature, opposite direction.
    nodes2(1, :) = [0.0_dp, 0.25_dp]
    nodes2(2, :) = [-0.5_dp, -0.25_dp]
    nodes2(3, :) = [1.0_dp, 0.25_dp]
    edges_first(2)%nodes = nodes2

    nodes2(1, :) = [0.75_dp, 0.25_dp]
    nodes2(2, :) = [-0.25_dp, -0.25_dp]
    nodes2(3, :) = [-0.25_dp, 0.25_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = (status == Status_SAME_CURVATURE)
    call print_status(name, case_id, case_success, success)

    ! CASE 13: Intersection is tangent, opposite direction, curvatures have
    !          same sign and there is no overlap.
    nodes2(1, :) = [2.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [1.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(1, :) = 1
    nodes2(2, :) = [1.5_dp, 0.0_dp]
    nodes2(3, :) = [2.0_dp, 1.0_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_OPPOSED .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 14: Intersection is tangent, opposite direction, curvatures have
    !          same sign and there **is** overlap.
    nodes2(1, :) = [1.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [2.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(1, :) = [2.0_dp, 1.0_dp]
    nodes2(2, :) = [1.5_dp, 0.0_dp]
    nodes2(3, :) = 1
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = (status == Status_BAD_TANGENT)
    call print_status(name, case_id, case_success, success)

    ! CASE 15: Intersection is tangent, opposite direction, curvatures have
    !          opposite sign and there is no overlap.
    nodes2(1, :) = [2.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [1.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(1, :) = [0.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [3.0_dp, 0.0_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_OPPOSED .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 16: Intersection is tangent, opposite direction, curvatures have
    !          opposite sign and there **is** overlap.
    nodes2(1, :) = [1.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [2.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes2

    nodes2(1, :) = [3.0_dp, 0.0_dp]
    nodes2(2, :) = [1.5_dp, 1.0_dp]
    nodes2(3, :) = [0.0_dp, 0.0_dp]
    edges_second(2)%nodes = nodes2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = (status == Status_BAD_TANGENT)
    call print_status(name, case_id, case_success, success)

    ! CASE 17: Intersection at corner, but the corner "staddles" an edge.
    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_first(2)%nodes = nodes1

    nodes1(1, :) = [1.0_dp, 1.0_dp]
    nodes1(2, :) = [0.5_dp, 0.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(1, :) = [0.5_dp, 0.0_dp]
    nodes1(2, :) = [1.5_dp, -1.0_dp]
    edges_second(2)%nodes = nodes1

    intersection_%s = 0.5_dp
    intersection_%index_first = 2
    intersection_%t = 0.0_dp
    intersection_%index_second = 2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 18: Intersection at corner of **both** edges, the corners just
    !          "kiss", hence are ignored.
    nodes1(1, :) = [0.5_dp, 1.0_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes1
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [1.5_dp, 0.25_dp]
    edges_first(2)%nodes = nodes1

    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [0.0_dp, 1.0_dp]
    edges_second(2)%nodes = nodes1

    intersection_%s = 0.0_dp
    intersection_%index_first = 2
    intersection_%t = 0.0_dp
    intersection_%index_second = 2

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_IGNORED_CORNER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 19: Intersection at corner of **both** edges, "first" is interior.
    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes1
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [0.0_dp, 1.0_dp]
    edges_first(2)%nodes = nodes1

    nodes1(1, :) = [0.5_dp, 0.25_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 1.0_dp]
    edges_second(2)%nodes = nodes1

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 20: Intersection at corner of **both** edges, "second" is interior.
    nodes1(1, :) = [0.5_dp, 0.25_dp]
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_first(1)%nodes = nodes1
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 1.0_dp]
    edges_first(2)%nodes = nodes1

    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_second(1)%nodes = nodes1
    nodes1(1, :) = [1.0_dp, 0.0_dp]
    nodes1(2, :) = [0.0_dp, 1.0_dp]
    edges_second(2)%nodes = nodes1

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 21: Intersection at corner of **both** edges, one corner is
    !          completely contained in the area of the other corner.
    nodes1(1, :) = [0.5_dp, 1.0_dp]
    nodes1(2, :) = 0
    edges_first(1)%nodes = nodes1
    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.5_dp]
    edges_first(2)%nodes = nodes1

    nodes1(1, :) = [0.0_dp, 1.0_dp]
    nodes1(2, :) = 0
    edges_second(1)%nodes = nodes1
    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    edges_second(2)%nodes = nodes1

    call classify_intersection( &
         edges_first, edges_second, intersection_, enum_, status)
    case_success = ( &
         enum_ == IntersectionClassification_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_classify_intersection

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
    integer(c_int) :: first, second
    integer(c_int) :: all_types, status

    case_id = 1
    name = "add_st_vals"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND

    ! CASE 1: ``intersections`` must be allocated.
    num_intersections = 0
    all_types = 0
    allocate(st_vals(2, 1))
    st_vals(:, 1) = [0.25_dp, 0.75_dp]
    case_success = .NOT. allocated(intersections)
    ! Set the relevant edges as lines.
    allocate(edges_first(2)%nodes(2, 2))
    allocate(edges_second(3)%nodes(2, 2))
    edges_first(2)%nodes(1, :) = [0.0_dp, 0.0_dp]
    edges_first(2)%nodes(2, :) = [1.0_dp, 1.0_dp]
    edges_second(3)%nodes(1, :) = [1.0_dp, -0.5_dp]
    edges_second(3)%nodes(2, :) = [0.0_dp, 0.5_dp]
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, &
         2, 3, intersections, num_intersections, all_types, status)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         size(intersections) == 1 .AND. &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.25_dp, 0.75_dp, 2, 3, second) .AND. &
         all_types == 2**second .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: ``intersections`` must be re-allocated.
    all_types = 0
    st_vals(:, 1) = [0.0_dp, 0.5_dp]
    ! Set the relevant edges as lines.
    allocate(edges_first(3)%nodes(2, 2))
    allocate(edges_second(1)%nodes(2, 2))
    edges_first(3)%nodes(1, :) = [0.0_dp, 0.0_dp]
    edges_first(3)%nodes(2, :) = [1.0_dp, 1.0_dp]
    edges_second(1)%nodes(1, :) = [1.0_dp, -1.0_dp]
    edges_second(1)%nodes(2, :) = [-1.0_dp, 1.0_dp]
    case_success = ( &
         allocated(intersections) .AND. &
         num_intersections == 1)
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, &
         3, 1, intersections, num_intersections, all_types, status)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         size(intersections) == 2 .AND. &
         num_intersections == 2 .AND. &
         intersection_check( &
         intersections(2), 0.0_dp, 0.5_dp, 3, 1, second) .AND. &
         all_types == 2**second .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: ``intersections`` is large enough.
    num_intersections = 0
    all_types = 0
    st_vals(:, 1) = [0.875_dp, 0.125_dp]
    ! Set the relevant edges as lines.
    allocate(edges_second(2)%nodes(2, 2))
    edges_first(2)%nodes(1, :) = [0.0_dp, 0.0_dp]
    edges_first(2)%nodes(2, :) = [1.0_dp, 1.0_dp]
    edges_second(2)%nodes(1, :) = [0.75_dp, 0.625_dp]
    edges_second(2)%nodes(2, :) = [1.75_dp, 2.625_dp]
    case_success = allocated(intersections)
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, &
         2, 2, intersections, num_intersections, all_types, status)
    case_success = ( &
         case_success .AND. &
         allocated(intersections) .AND. &
         size(intersections) == 2 .AND. &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.875_dp, 0.125_dp, 2, 2, second) .AND. &
         all_types == 2**second .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Intersection causes error.
    all_types = 0
    st_vals(:, 1) = [0.5_dp, 0.5_dp]
    deallocate(edges_first(2)%nodes)
    allocate(edges_first(2)%nodes(3, 2))
    edges_first(2)%nodes(1, :) = [1.0_dp, 0.25_dp]
    edges_first(2)%nodes(2, :) = [-0.5_dp, -0.25_dp]
    edges_first(2)%nodes(3, :) = [0.0_dp, 0.25_dp]
    deallocate(edges_second(2)%nodes)
    allocate(edges_second(2)%nodes(3, 2))
    edges_second(2)%nodes(1, :) = [0.75_dp, 0.25_dp]
    edges_second(2)%nodes(2, :) = [-0.25_dp, -0.25_dp]
    edges_second(2)%nodes(3, :) = [-0.25_dp, 0.25_dp]
    call add_st_vals( &
         edges_first, edges_second, 1, st_vals, &
         2, 2, intersections, num_intersections, all_types, status)
    case_success = ( &
         all_types == 0 .AND. &
         status == Status_SAME_CURVATURE)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: ``intersections`` has edge ends, gets "over-allocated".
    deallocate(intersections)
    deallocate(st_vals)
    allocate(st_vals(2, 2))
    num_intersections = 0
    all_types = 0
    st_vals(:, 1) = [1.0_dp, 0.125_dp]
    st_vals(:, 2) = [0.5_dp, 0.5_dp]
    ! Set the relevant edges as lines.
    allocate(edges_first(1)%nodes(2, 2))
    edges_first(1)%nodes(1, :) = [0.0_dp, 0.0_dp]
    edges_first(1)%nodes(2, :) = [1.0_dp, 1.0_dp]
    edges_second(1)%nodes(1, :) = [0.0_dp, 1.0_dp]
    edges_second(1)%nodes(2, :) = [1.0_dp, 0.0_dp]
    call add_st_vals( &
         edges_first, edges_second, 2, st_vals, &
         1, 1, intersections, num_intersections, all_types, status)
    case_success = ( &
         allocated(intersections) .AND. &
         size(intersections) == 2 .AND. &
         num_intersections == 1 .AND. &
         intersection_check( &
         intersections(1), 0.5_dp, 0.5_dp, 1, 1, first) .AND. &
         all_types == 2**first .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_add_st_vals

  subroutine test_surfaces_intersection_points(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(28) :: name
    real(c_double) :: linear1(3, 2), linear2(3, 2)
    real(c_double) :: quadratic1(6, 2), quadratic2(6, 2)
    type(Intersection), allocatable :: intersections(:)
    integer(c_int) :: num_intersections
    integer(c_int) :: all_types, status
    integer(c_int) :: first, second

    case_id = 1
    name = "surfaces_intersection_points"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND

    ! CASE 1: Overlapping triangles (i.e. degree 1 surfaces).
    linear1(1, :) = 0
    linear1(2, :) = [2.0_dp, 0.0_dp]
    linear1(3, :) = [1.0_dp, 2.0_dp]
    linear2(1, :) = [0.0_dp, 1.0_dp]
    linear2(2, :) = [2.0_dp, 1.0_dp]
    linear2(3, :) = [1.0_dp, -1.0_dp]
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
         all_types == (2**first + 2**second) .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Tangent intersection (causes ``all_intersections()`` to error).
    quadratic1(1, :) = 0
    quadratic1(2, :) = [6.0_dp, 12.0_dp]
    quadratic1(3, :) = [12.0_dp, 6.0_dp]
    quadratic1(4, :) = [4.0_dp, -8.0_dp]
    quadratic1(5, :) = [10.0_dp, -5.0_dp]
    quadratic1(6, :) = [8.0_dp, -16.0_dp]
    quadratic2(1, :) = [4.0_dp, 10.0_dp]
    quadratic2(2, :) = [10.0_dp, 4.0_dp]
    quadratic2(3, :) = [16.0_dp, 16.0_dp]
    quadratic2(4, :) = [6.0_dp, 21.0_dp]
    quadratic2(5, :) = [12.0_dp, 24.0_dp]
    quadratic2(6, :) = [8.0_dp, 32.0_dp]
    call surfaces_intersection_points( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         intersections, num_intersections, all_types, status)
    case_success = ( &
         num_intersections == 0 .AND. &
         all_types == 0 .AND. &
         status == Status_PARALLEL)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Tangent intersection (causes ``add_st_vals()`` to error).
    num_intersections = 0
    quadratic1(1, :) = 0
    quadratic1(2, :) = [2.0_dp, 4.0_dp]
    quadratic1(3, :) = [4.0_dp, 0.0_dp]
    quadratic1(4, :) = [1.0_dp, 3.0_dp]
    quadratic1(5, :) = [3.0_dp, 3.0_dp]
    quadratic1(6, :) = [2.0_dp, 6.0_dp]
    quadratic2(1, :) = [4.0_dp, 4.0_dp]
    quadratic2(2, :) = [2.0_dp, 0.0_dp]
    quadratic2(3, :) = [0.0_dp, 4.0_dp]
    quadratic2(4, :) = [3.0_dp, 1.0_dp]
    quadratic2(5, :) = [1.0_dp, 1.0_dp]
    quadratic2(6, :) = [2.0_dp, -2.0_dp]
    call surfaces_intersection_points( &
         6, quadratic1, 2, 6, quadratic2, 2, &
         intersections, num_intersections, all_types, status)
    case_success = ( &
         num_intersections == 1 .AND. &
         all_types == 0 .AND. &
         status == Status_BAD_TANGENT)
    call print_status(name, case_id, case_success, success)

  end subroutine test_surfaces_intersection_points

  function intersection_check( &
       intersection_, s, t, index_first, index_second, &
       interior_curve) result(same)

    type(intersection), intent(in) :: intersection_
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

  subroutine test_remove_node(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(11) :: name
    integer(c_int) :: values(4)
    integer(c_int) :: remaining

    case_id = 1
    name = "remove_node"

    ! CASE 1: Not contained.
    values = [1, 3, 4, 8]
    remaining = 4
    call remove_node(2, values, remaining)
    case_success = ( &
         all(values == [1, 3, 4, 8]) .AND. &
         remaining == 4)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Is contained, not at end.
    values = [2, 3, 7, 8]
    remaining = 4
    call remove_node(3, values, remaining)
    case_success = ( &
         all(values == [2, 7, 8, 8]) .AND. &
         remaining == 3)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Is contained, from end.
    values = [3, 5, 8, 11]
    remaining = 4
    call remove_node(11, values, remaining)
    case_success = ( &
         all(values == [3, 5, 8, 11]) .AND. &
         remaining == 3)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: ``remaining`` is less than ``size(values)``, not contained.
    values = [11, 13, 18, -1]
    remaining = 3
    call remove_node(15, values, remaining)
    case_success = ( &
         all(values == [11, 13, 18, -1]) .AND. &
         remaining == 3)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: ``remaining`` is less than ``size(values)``, is contained.
    values = [1, 8, 9, -1]
    remaining = 3
    call remove_node(8, values, remaining)
    case_success = ( &
         all(values == [1, 9, 9, -1]) .AND. &
         remaining == 2)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: ``remaining`` is less than ``size(values)``, at end.
    values = [2, 3, 6, -1]
    remaining = 3
    call remove_node(6, values, remaining)
    case_success = ( &
         all(values == [2, 3, 6, -1]) .AND. &
         remaining == 2)
    call print_status(name, case_id, case_success, success)

  end subroutine test_remove_node

  subroutine test_get_next(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(8) :: name
    integer(c_int) :: first, second
    type(Intersection) :: intersections(5)
    type(SegmentNode) :: curr_node, next_node
    logical(c_bool) :: at_start
    integer(c_int) :: unused(4)
    integer(c_int) :: remaining

    case_id = 1
    name = "get_next"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND

    ! CASE 1: On an edge of first surface, no other intersections on that edge.
    curr_node%edge_index = 2
    curr_node%edge_param = 0.125_dp
    intersections(1)%index_first = 3
    intersections(1)%s = 0.5_dp
    intersections(1)%interior_curve = first
    remaining = 4
    call get_next( &
         1, intersections(:1), unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == curr_node%edge_index .AND. &
         next_node%edge_param == 1.0_dp .AND. &
         next_node%interior_curve == first .AND. &
         remaining == 4 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: On an edge of first surface, **multiple** other intersections
    !         on same edge.
    curr_node%edge_index = 3
    curr_node%edge_param = 0.25_dp
    ! An "acceptable" intersection that will be overtaken by the
    ! next since 0.25 < 0.5 < 0.875.
    intersections(1)%index_first = 3
    intersections(1)%s = 0.875_dp
    intersections(1)%interior_curve = second
    ! Closer to ``curr_node`` then previous.
    intersections(2)%index_first = 3
    intersections(2)%s = 0.5_dp
    intersections(2)%interior_curve = first
    ! On a different edge.
    intersections(3)%index_first = 2
    ! Same edge, but parameter comes **before** ``curr_node``.
    intersections(4)%index_first = 3
    intersections(4)%s = 0.125_dp
    intersections(4)%interior_curve = first
    ! Past the already accepted intersection: 0.25 < 0.5 < 0.625
    intersections(5)%index_first = 3
    intersections(5)%s = 0.625_dp
    intersections(5)%interior_curve = first
    ! Populate the list of indices to be removed.
    unused = [1, 2, 3, 5]
    remaining = 4
    call get_next( &
         5, intersections, unused, remaining, &
         2, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == intersections(2)%index_first .AND. &
         next_node%edge_param == intersections(2)%s .AND. &
         next_node%interior_curve == intersections(2)%interior_curve .AND. &
         remaining == 3 .AND. &
         all(unused == [1, 3, 5, 5]) .AND. &
         at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: On an edge of first surface, move to a corner that is **also** an
    !         intersection.
    curr_node%edge_index = 1
    curr_node%edge_param = 0.625_dp
    intersections(1)%index_first = 1
    intersections(1)%s = 1.0_dp
    intersections(1)%interior_curve = IntersectionClassification_TANGENT_FIRST
    unused(1) = 2  ! Won't contain ``intersection_index == 1``.
    remaining = 1
    call get_next( &
         1, intersections(:1), unused, remaining, &
         -1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == intersections(1)%index_first .AND. &
         next_node%edge_param == intersections(1)%s .AND. &
         next_node%interior_curve == intersections(1)%interior_curve .AND. &
         remaining == 1 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: On an edge of second surface, no other intersections on
    !         that edge.
    curr_node%edge_index = 6
    curr_node%edge_param = 0.625_dp
    intersections(1)%index_second = 1
    intersections(1)%s = 0.5_dp
    intersections(1)%interior_curve = first
    remaining = 2
    call get_next( &
         1, intersections(:1), unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == curr_node%edge_index .AND. &
         next_node%edge_param == 1.0_dp .AND. &
         next_node%interior_curve == second .AND. &
         remaining == 2 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: On an edge of second surface, **multiple** other intersections
    !         on same edge.
    curr_node%edge_index = 5
    curr_node%edge_param = 0.125_dp
    ! An "acceptable" intersection that will be overtaken by the
    ! next since 0.125 < 0.625 < 0.75.
    intersections(1)%index_second = 2
    intersections(1)%t = 0.75_dp
    intersections(1)%interior_curve = first
    ! On a different edge.
    intersections(2)%index_second = 3
    ! Closer to ``curr_node`` then previous accepted.
    intersections(3)%index_second = 2
    intersections(3)%t = 0.625_dp
    intersections(3)%interior_curve = second
    ! Same edge, but parameter comes **before** ``curr_node``.
    intersections(4)%index_second = 2
    intersections(4)%t = 0.0625_dp
    intersections(4)%interior_curve = first
    ! Past the already accepted intersection: 0.125 < 0.625 < 0.6875
    intersections(5)%index_second = 2
    intersections(5)%t = 0.6875_dp
    intersections(5)%interior_curve = second
    ! Populate the list of indices to be removed, missing 3.
    unused(:3) = [1, 2, 4]
    remaining = 3
    call get_next( &
         5, intersections, unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == intersections(3)%index_second + 3 .AND. &
         next_node%edge_param == intersections(3)%t .AND. &
         next_node%interior_curve == intersections(3)%interior_curve .AND. &
         remaining == 3 .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: On an edge of second surface, move to a corner that is **also**
    !         an intersection.
    curr_node%edge_index = 4
    curr_node%edge_param = 0.5_dp
    intersections(1)%index_second = 1
    intersections(1)%t = 1.0_dp
    intersections(1)%interior_curve = IntersectionClassification_TANGENT_FIRST
    unused = [1, -1, -1, -1]  ! Will be removed.
    remaining = 1
    call get_next( &
         1, intersections(:1), unused, remaining, &
         1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == intersections(1)%index_second + 3 .AND. &
         next_node%edge_param == intersections(1)%t .AND. &
         next_node%interior_curve == intersections(1)%interior_curve .AND. &
         remaining == 0 .AND. &
         all(unused == [1, -1, -1, -1]) .AND. &
         at_start)
    call print_status(name, case_id, case_success, success)

  end subroutine test_get_next

  subroutine test_to_front(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(8) :: name
    integer(c_int) :: first, second
    type(Intersection) :: intersections(1)
    type(SegmentNode) :: curr_node, next_node
    logical(c_bool) :: at_start

    case_id = 1
    name = "to_front"

    first = IntersectionClassification_FIRST
    second = IntersectionClassification_SECOND

    ! CASE 1: Unchanged node (i.e. not at the end).
    curr_node%edge_param = 0.5_dp
    curr_node%edge_index = 3
    curr_node%interior_curve = first
    call to_front( &
         0, intersections(:0), -1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == curr_node%edge_index .AND. &
         next_node%edge_param == curr_node%edge_param .AND. &
         next_node%interior_curve == curr_node%interior_curve .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Change the ``s`` parameter, no corresponding node.
    curr_node%edge_param = 1.0_dp
    curr_node%edge_index = 3
    curr_node%interior_curve = first
    call to_front( &
         0, intersections(:0), -1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == curr_node%edge_index - 2 .AND. &
         next_node%edge_param == 0.0_dp .AND. &
         next_node%interior_curve == curr_node%interior_curve .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Change the ``s`` parameter, new corner node exists.
    curr_node%edge_param = 1.0_dp
    curr_node%edge_index = 2
    curr_node%interior_curve = first
    intersections(1)%index_first = 3
    intersections(1)%s = 0.0
    intersections(1)%index_second = 1
    intersections(1)%t = 0.5
    intersections(1)%interior_curve = IntersectionClassification_TANGENT_FIRST
    call to_front( &
         1, intersections, -1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == intersections(1)%index_first .AND. &
         next_node%edge_param == intersections(1)%s .AND. &
         next_node%interior_curve == intersections(1)%interior_curve .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Change the ``t`` parameter, no corresponding node.
    curr_node%edge_param = 1.0_dp
    curr_node%edge_index = 5
    curr_node%interior_curve = second
    call to_front( &
         0, intersections(:0), -1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == curr_node%edge_index + 1 .AND. &
         next_node%edge_param == 0.0_dp .AND. &
         next_node%interior_curve == curr_node%interior_curve .AND. &
         .NOT. at_start)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Change the ``t`` parameter, new corner node exists.
    curr_node%edge_param = 1.0_dp
    curr_node%edge_index = 4
    curr_node%interior_curve = second
    intersections(1)%index_first = 3
    intersections(1)%s = 0.5_dp
    intersections(1)%index_second = 2
    intersections(1)%t = 0.0_dp
    intersections(1)%interior_curve = IntersectionClassification_TANGENT_SECOND
    call to_front( &
         1, intersections, 1, curr_node, next_node, at_start)
    case_success = ( &
         next_node%edge_index == intersections(1)%index_second + 3 .AND. &
         next_node%edge_param == intersections(1)%t .AND. &
         next_node%interior_curve == intersections(1)%interior_curve .AND. &
         at_start)
    call print_status(name, case_id, case_success, success)

  end subroutine test_to_front

  subroutine test_surfaces_intersect(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(18) :: name
    integer(c_int) :: contained, status
    real(c_double) :: linear1(3, 2), linear2(3, 2)
    real(c_double) :: quadratic1(6, 2), quadratic2(6, 2)

    case_id = 1
    name = "surfaces_intersect"

    ! CASE 1: Disjoint bounding boxes.
    linear1(1, :) = 0
    linear1(2, :) = [1.0_dp, 0.0_dp]
    linear1(3, :) = [0.0_dp, 1.0_dp]
    linear2(1, :) = 10
    linear2(2, :) = [11.0_dp, 10.0_dp]
    linear2(3, :) = [10.0_dp, 11.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, contained, status)
    case_success = ( &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Intersecting edges results in error.
    linear1(1, :) = 0
    linear1(2, :) = [5.0_dp, 0.0_dp]
    linear1(3, :) = [0.0_dp, 5.0_dp]
    linear2(1, :) = [2.0_dp, 0.0_dp]
    linear2(2, :) = [3.0_dp, 0.0_dp]
    linear2(3, :) = [2.0_dp, 1.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, contained, status)
    case_success = ( &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_PARALLEL)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Surface 2 contained in surface 1 (re-uses ``linear1``).
    linear2(1, :) = [2.0_dp, 1.0_dp]
    linear2(2, :) = [3.0_dp, 1.0_dp]
    linear2(3, :) = [2.0_dp, 2.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, contained, status)
    case_success = ( &
         contained == SurfaceContained_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Surface 1 contained in surface 2 (re-uses all data from CASE 3).
    call surfaces_intersect( &
         3, linear2, 1, 3, linear1, 1, contained, status)
    case_success = ( &
         contained == SurfaceContained_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Surfaces disjoint with overlapping bounding boxes (re-uses
    !         ``linear1`` from previous cases).
    linear2(1, :) = [4.0_dp, 2.0_dp]
    linear2(2, :) = [4.0_dp, 3.0_dp]
    linear2(3, :) = [3.0_dp, 3.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, contained, status)
    case_success = ( &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Surfaces intersect (re-uses ``linear1`` from previous cases).
    linear2(1, :) = [3.0_dp, -2.0_dp]
    linear2(2, :) = [3.0_dp, 3.0_dp]
    linear2(3, :) = [-2.0_dp, 3.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, contained, status)
    case_success = (status == -1234)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: The only intersection(s) are ``OPPOSED``.
    quadratic1(1, :) = [4.0_dp, 0.0_dp]
    quadratic1(2, :) = [2.0_dp, 4.0_dp]
    quadratic1(3, :) = [0.0_dp, 0.0_dp]
    quadratic1(4, :) = [3.0_dp, -2.0_dp]
    quadratic1(5, :) = [1.0_dp, -2.0_dp]
    quadratic1(6, :) = [2.0_dp, -4.0_dp]
    quadratic2(1, :) = [0.0_dp, 4.0_dp]
    quadratic2(2, :) = [2.0_dp, 0.0_dp]
    quadratic2(3, :) = [4.0_dp, 4.0_dp]
    quadratic2(4, :) = [1.0_dp, 6.0_dp]
    quadratic2(5, :) = [3.0_dp, 6.0_dp]
    quadratic2(6, :) = [2.0_dp, 8.0_dp]

    call surfaces_intersect( &
         6, quadratic1, 2, 6, quadratic2, 2, contained, status)
    case_success = ( &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: The only intersection(s) are ``IGNORED_CORNER``.
    linear1(1, :) = 0
    linear1(2, :) = [1.0_dp, 0.0_dp]
    linear1(3, :) = [0.0_dp, 1.0_dp]
    linear2(1, :) = 0
    linear2(2, :) = [-2.0_dp, 3.0_dp]
    linear2(3, :) = [1.0_dp, -3.0_dp]

    call surfaces_intersect( &
         3, linear1, 1, 3, linear2, 1, contained, status)
    case_success = ( &
         contained == SurfaceContained_NEITHER .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: The only intersection(s) are ``TANGENT_FIRST``.
    quadratic1(1, :) = [2.0_dp, 1.25_dp]
    quadratic1(2, :) = [4.0_dp, 0.75_dp]
    quadratic1(3, :) = [6.0_dp, 1.25_dp]
    quadratic1(4, :) = [3.0_dp, 3.0_dp]
    quadratic1(5, :) = [5.0_dp, 3.0_dp]
    quadratic1(6, :) = [4.0_dp, 5.0_dp]
    quadratic2(1, :) = [0.0_dp, 0.0_dp]
    quadratic2(2, :) = [4.0_dp, 2.0_dp]
    quadratic2(3, :) = [8.0_dp, 0.0_dp]
    quadratic2(4, :) = [1.5_dp, 7.5_dp]
    quadratic2(5, :) = [11.0_dp, 6.0_dp]
    quadratic2(6, :) = [8.0_dp, 8.0_dp]

    call surfaces_intersect( &
         6, quadratic1, 2, 6, quadratic2, 2, contained, status)
    case_success = ( &
         contained == SurfaceContained_FIRST .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: The only intersection(s) are ``TANGENT_SECOND`` (re-uses
    !          **all** data from CASE 9).
    call surfaces_intersect( &
         6, quadratic2, 2, 6, quadratic1, 2, contained, status)
    case_success = ( &
         contained == SurfaceContained_SECOND .AND. &
         status == Status_SUCCESS)
    call print_status(name, case_id, case_success, success)

  end subroutine test_surfaces_intersect

end module test_surface_intersection
