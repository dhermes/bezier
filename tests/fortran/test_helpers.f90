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

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_bool
  use helpers, only: &
       WIGGLE, cross_product, bbox, wiggle_interval, contains_nd, &
       vector_close, in_interval, convex_hull, polygon_collide, &
       solve2x2
  use types, only: dp
  use unit_test_helpers, only: MACHINE_EPS, print_status
  implicit none
  private &
       test_cross_product, test_bbox, test_wiggle_interval, &
       test_contains_nd, test_vector_close, test_in_interval, &
       test_convex_hull, test_polygon_collide, test_solve2x2
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
    call test_convex_hull(success)
    call test_polygon_collide(success)
    call test_solve2x2(success)

  end subroutine helpers_all_tests

  subroutine test_cross_product(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: vec0(2)
    real(c_double) :: vec1(2)
    real(c_double) :: result_
    integer :: case_id
    character(13) :: name

    case_id = 1
    name = "cross_product"

    ! CASE 1: Identical vector.
    vec0 = [1.0_dp, 7.0_dp] / 8
    vec1 = [-11.0_dp, 24.0_dp] / 32
    call cross_product(vec0, vec1, result_)
    case_success = (result_ == 101.0_dp / 256)
    call print_status(name, case_id, case_success, success)

  end subroutine test_cross_product

  subroutine test_bbox(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 2), nodes2(2, 6)
    real(c_double) :: left, right, bottom, top
    integer :: case_id
    character(4) :: name

    case_id = 1
    name = "bbox"

    ! CASE 1: Simple case (just two nodes).
    nodes1(:, 1) = [0.0_dp, 5.0_dp]
    nodes1(:, 2) = [1.0_dp, 3.0_dp]
    call bbox(2, nodes1, left, right, bottom, top)
    case_success = ( &
         left == 0.0_dp .AND. right == 1.0_dp &
         .AND. bottom == 3.0_dp .AND. top == 5.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: "Many" nodes.
    nodes2(:, 1) = [1.0_dp, 0.0_dp]
    nodes2(:, 2) = [2.0_dp, 1.0_dp]
    nodes2(:, 3) = [-1.0_dp, 2.0_dp]
    nodes2(:, 4) = [5.0_dp, -3.0_dp]
    nodes2(:, 5) = [4.0_dp, 4.0_dp]
    nodes2(:, 6) = [0.0_dp, 0.0_dp]
    call bbox(6, nodes2, left, right, bottom, top)
    case_success = ( &
         left == -1.0_dp .AND. right == 5.0_dp &
         .AND. bottom == -3.0_dp .AND. top == 4.0_dp)
    call print_status(name, case_id, case_success, success)

  end subroutine test_bbox

  subroutine test_wiggle_interval(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: value1, value2, value3, value4, value5
    real(c_double) :: result1, result2, result3, result4, result5
    logical(c_bool) :: &
         wiggle_success1, wiggle_success2, wiggle_success3, &
         wiggle_success4, wiggle_success5
    integer :: case_id
    character(15) :: name

    case_id = 1
    name = "wiggle_interval"

    ! CASE 1: **Exactly** at endpoint.
    value1 = 0.0_dp
    call wiggle_interval(value1, result1, wiggle_success1)
    case_success = (wiggle_success1 .AND. result1 == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Near endpoint.
    value1 = 1.0_dp + 0.5_dp**20
    call wiggle_interval(value1, result1, wiggle_success1)
    case_success = (.NOT. wiggle_success1)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: "Below" the interval.
    value1 = -0.25_dp
    call wiggle_interval(value1, result1, wiggle_success1)
    case_success = (.NOT. wiggle_success1)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: "Above" the interval.
    value1 = 1.5_dp
    call wiggle_interval(value1, result1, wiggle_success1)
    case_success = (.NOT. wiggle_success1)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Inside the interval.
    value1 = 0.25_dp
    call wiggle_interval(value1, result1, wiggle_success1)
    case_success = (wiggle_success1 .AND. result1 == value1)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Wiggle "below" the interval.
    value1 = -0.5_dp**60
    call wiggle_interval(value1, result1, wiggle_success1)
    case_success = (wiggle_success1 .AND. result1 == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Wiggle "above" the interval.
    value1 = 1.0_dp + MACHINE_EPS
    call wiggle_interval(value1, result1, wiggle_success1)
    case_success = (wiggle_success1 .AND. result1 == 1.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Test "lower outer" boundary, i.e. values in (-epsilon, 0).
    value1 = -WIGGLE + MACHINE_EPS * WIGGLE
    value2 = -WIGGLE + MACHINE_EPS * WIGGLE / 2
    value3 = -WIGGLE
    value4 = -WIGGLE - MACHINE_EPS * WIGGLE
    call wiggle_interval(value1, result1, wiggle_success1)
    call wiggle_interval(value2, result2, wiggle_success2)
    call wiggle_interval(value3, result3, wiggle_success3)
    call wiggle_interval(value4, result4, wiggle_success4)
    case_success = ( &
         wiggle_success1 .AND. result1 == 0.0_dp .AND. &
         wiggle_success2 .AND. result2 == 0.0_dp .AND. &
         .NOT. wiggle_success3 .AND. .NOT. wiggle_success4)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Test "upper outer" boundary, i.e. values in (1, 1 + epsilon).
    value1 = 1.0_dp + WIGGLE - 2 * MACHINE_EPS
    value2 = 1.0_dp + WIGGLE - MACHINE_EPS
    value3 = 1.0_dp + WIGGLE
    value4 = 1.0_dp + WIGGLE + MACHINE_EPS
    call wiggle_interval(value1, result1, wiggle_success1)
    call wiggle_interval(value2, result2, wiggle_success2)
    call wiggle_interval(value3, result3, wiggle_success3)
    call wiggle_interval(value4, result4, wiggle_success4)
    case_success = ( &
         wiggle_success1 .AND. result1 == 1.0_dp .AND. &
         wiggle_success2 .AND. result2 == 1.0_dp .AND. &
         .NOT. wiggle_success3 .AND. .NOT. wiggle_success4)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: Test "lower inner" boundary, i.e. values in (0, epsilon).
    value1 = WIGGLE - WIGGLE * MACHINE_EPS
    value2 = WIGGLE - WIGGLE * MACHINE_EPS / 2
    value3 = WIGGLE
    value4 = WIGGLE + WIGGLE * MACHINE_EPS
    call wiggle_interval(value1, result1, wiggle_success1)
    call wiggle_interval(value2, result2, wiggle_success2)
    call wiggle_interval(value3, result3, wiggle_success3)
    call wiggle_interval(value4, result4, wiggle_success4)
    case_success = ( &
         wiggle_success1 .AND. result1 == 0.0_dp .AND. &
         wiggle_success2 .AND. result2 == 0.0_dp .AND. &
         wiggle_success3 .AND. result3 == value3 .AND. &
         wiggle_success4 .AND. result4 == value4)
    call print_status(name, case_id, case_success, success)

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
    case_success = ( &
         wiggle_success1 .AND. result1 == value1 .AND. &
         wiggle_success2 .AND. result2 == value2 .AND. &
         wiggle_success3 .AND. result3 == value3 .AND. &
         wiggle_success4 .AND. result4 == 1.0_dp .AND. &
         wiggle_success5 .AND. result5 == 1.0_dp)
    call print_status(name, case_id, case_success, success)

  end subroutine test_wiggle_interval

  subroutine test_contains_nd(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 3), nodes2(3, 2), nodes3(4, 2)
    real(c_double) :: point1(2), point2(3), point3(4)
    logical(c_bool) :: is_contained
    integer :: case_id
    character(11) :: name

    case_id = 1
    name = "contains_nd"

    ! CASE 1: Below bounding box.
    nodes1(:, 1) = [0.0_dp, 1.0_dp]
    nodes1(:, 2) = [0.5_dp, 0.0_dp]
    nodes1(:, 3) = [1.0_dp, 2.0_dp]
    point1 = [-0.5_dp, 1.0_dp]
    call contains_nd(3, 2, nodes1, point1, is_contained)
    case_success = (.NOT. is_contained)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Above bounding box.
    nodes2(:, 1) = [0.0_dp, -4.0_dp, 2.0_dp]
    nodes2(:, 2) = [-1.0_dp, 1.0_dp, 3.0_dp]
    point2 = [-0.5_dp, 2.0_dp, 2.5_dp]
    call contains_nd(2, 3, nodes2, point2, is_contained)
    case_success = (.NOT. is_contained)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Inside bounding box.
    nodes3(:, 1) = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp]
    nodes3(:, 2) = [1.0_dp, -2.0_dp, -4.0_dp, 1.0_dp]
    point3 = [0.5_dp, 0.0_dp, 0.0_dp, 2.0_dp]
    call contains_nd(2, 4, nodes3, point3, is_contained)
    case_success = (is_contained)
    call print_status(name, case_id, case_success, success)

  end subroutine test_contains_nd

  subroutine test_vector_close(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: eps
    real(c_double) :: vec1(2)
    real(c_double) :: vec2(2)
    logical(c_bool) :: is_close
    integer :: case_id
    character(12) :: name

    eps = 0.5_dp**40
    case_id = 1
    name = "vector_close"

    ! CASE 1: Identical vector.
    vec1 = [0.5_dp, 4.0_dp]
    is_close = vector_close(2, vec1, vec1, eps)
    case_success = (is_close)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Far apart vectors.
    vec1 = [0.0_dp, 6.0_dp]
    vec2 = [1.0_dp, -4.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    case_success = (.NOT. is_close)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Close but different.
    vec1 = [2.25_dp, -3.5_dp]
    vec2 = vec1 + 0.5_dp**43 * [-5.0_dp, 12.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    case_success = (is_close)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Custom epsilon.
    vec1 = [3.0_dp, 4.0_dp]
    vec2 = [2.0_dp, 5.0_dp]
    is_close = vector_close(2, vec1, vec2, 0.5_dp)
    case_success = (is_close .AND. .NOT. vector_close(2, vec1, vec2, eps))
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Near zero.
    vec1 = [0.0_dp, 0.0_dp]
    vec2 = 0.5_dp**45 * [3.0_dp, 4.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    case_success = (is_close)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Near zero failure (i.e. not near enough).
    vec1 = 0.5_dp**20 * [1.0_dp, 0.0_dp]
    vec2 = [0.0_dp, 0.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    case_success = (.NOT. is_close)
    call print_status(name, case_id, case_success, success)

  end subroutine test_vector_close

  subroutine test_in_interval(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    logical(c_bool) :: is_inside
    integer :: case_id
    character(11) :: name

    case_id = 1
    name = "in_interval"

    ! CASE 1: Interior value.
    is_inside = in_interval(1.5_dp, 1.0_dp, 2.0_dp)
    case_success = (is_inside)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Barely inside.
    is_inside = in_interval(1.0_dp + MACHINE_EPS, 1.0_dp, 2.0_dp)
    case_success = (is_inside)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Barely outside.
    is_inside = in_interval(1.0_dp - MACHINE_EPS / 2, 1.0_dp, 2.0_dp)
    case_success = (.NOT. is_inside)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Exterior value.
    is_inside = in_interval(-1.0_dp, 1.0_dp, 2.0_dp)
    case_success = (.NOT. is_inside)
    call print_status(name, case_id, case_success, success)

  end subroutine test_in_interval

  subroutine test_convex_hull(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    real(c_double) :: points1(2, 5), points2(2, 2), points3(2, 100)
    integer(c_int) :: polygon_size, i, j
    real(c_double) :: polygon(2, 100)
    character(11) :: name

    case_id = 1
    name = "convex_hull"

    ! CASE 1: Triangle with centroid (i.e. inside) and a repeated corner.
    points1(1, :) = [0.0_dp, 0.0_dp, 1.0_dp, 3.0_dp, 0.0_dp]
    points1(2, :) = [0.0_dp, 3.0_dp, 1.0_dp, 0.0_dp, 3.0_dp]
    call convex_hull(5, points1, polygon_size, polygon(:, :5))
    case_success = ( &
         polygon_size == 3 .AND. &
         all(polygon(1, :3) == [0.0_dp, 3.0_dp, 0.0_dp]) .AND. &
         all(polygon(2, :3) == [0.0_dp, 0.0_dp, 3.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Degenerate case (2 points, i.e. the "polygon" is just a line).
    points2(1, :) = [0.0_dp, 1.0_dp]
    points2(2, :) = [0.0_dp, 0.0_dp]
    call convex_hull(2, points2, polygon_size, polygon(:, :2))
    case_success = ( &
         polygon_size == 2 .AND. &
         all(polygon(1, :2) == [0.0_dp, 1.0_dp]) .AND. &
         all(polygon(2, :2) == [0.0_dp, 0.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Same as CASE 2, but points are in reverse order.
    points2(1, :) = [1.0_dp, 0.0_dp]
    points2(2, :) = [0.0_dp, 0.0_dp]
    call convex_hull(2, points2, polygon_size, polygon(:, :2))
    case_success = ( &
         polygon_size == 2 .AND. &
         all(polygon(1, :2) == [0.0_dp, 1.0_dp]) .AND. &
         all(polygon(2, :2) == [0.0_dp, 0.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 4: 10-by-10 grid of points. Also covers case where ``polygon`` is
    !         allocated but not large enough.
    do i = 1, 10
       do j = 1, 10
          points3(1, 10 * (i - 1) + j) = j
          points3(2, 10 * (i - 1) + j) = i
       end do
    end do
    call convex_hull(100, points3, polygon_size, polygon)
    case_success = ( &
         polygon_size == 4 .AND. &
         all(polygon(1, :4) == [1.0_dp, 10.0_dp, 10.0_dp, 1.0_dp]) .AND. &
         all(polygon(2, :4) == [1.0_dp, 1.0_dp, 10.0_dp, 10.0_dp]))
    call print_status(name, case_id, case_success, success)

    ! CASE 5: No points (degenerate case).
    call convex_hull(0, points1(:, :0), polygon_size, polygon(:, :0))
    case_success = (polygon_size == 0)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: One point (degenerate case).
    points1(:, 1) = [2.0_dp, 3.0_dp]
    call convex_hull(1, points1(:, :1), polygon_size, polygon(:, :1))
    case_success = ( &
         polygon_size == 1 .AND. &
         all(polygon(1, :1) == [2.0_dp]) .AND. &
         all(polygon(2, :1) == [3.0_dp]))
    call print_status(name, case_id, case_success, success)

  end subroutine test_convex_hull

  subroutine test_polygon_collide(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    real(c_double) :: polygon1(2, 3), polygon2(2, 4), polygon3(2, 4)
    logical(c_bool) :: collision
    character(15) :: name

    case_id = 1
    name = "polygon_collide"

    ! CASE 1: Triangle and square that do not collide.
    polygon1(1, :) = [1.0_dp, 3.0_dp, -2.0_dp]
    polygon1(2, :) = [1.0_dp, 4.0_dp, 2.0_dp]
    polygon2(1, :) = [3.5_dp, 6.5_dp, 6.5_dp, 3.5_dp]
    polygon2(2, :) = [3.0_dp, 3.0_dp, 7.0_dp, 7.0_dp]
    call polygon_collide(3, polygon1, 4, polygon2, collision)
    case_success = (.NOT. collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Same as CASE 1 but with polygons swapped.
    call polygon_collide(4, polygon2, 3, polygon1, collision)
    case_success = (.NOT. collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Triangles that do collide.
    polygon1(1, :) = [0.0_dp, 3.0_dp, 0.0_dp]
    polygon1(2, :) = [0.0_dp, 0.0_dp, 3.0_dp]
    polygon2(1, :3) = [1.0_dp, 4.0_dp, 1.0_dp]
    polygon2(2, :3) = [1.0_dp, 1.0_dp, 4.0_dp]
    call polygon_collide(3, polygon1, 3, polygon2(:, :3), collision)
    case_success = (collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Two squares where an edge other than edge 1 in polygon 1 is
    !         the separating line.
    polygon2(1, :) = [1.0_dp, 1.0_dp, 0.0_dp, 0.0_dp]
    polygon2(2, :) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
    polygon3(1, :) = [2.0_dp, 3.0_dp, 3.0_dp, 2.0_dp]
    polygon3(2, :) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
    call polygon_collide(4, polygon2, 4, polygon3, collision)
    case_success = (.NOT. collision)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Square and triangle where none of the edges in polygon 1 are
    !         separating lines and an edge other than edge 1 in polygon 2
    !         **is** the separating line.
    polygon2(1, :) = [0.0_dp, 2.0_dp, 2.0_dp, 0.0_dp]
    polygon2(2, :) = [0.0_dp, 0.0_dp, 2.0_dp, 2.0_dp]
    polygon1(1, :) = [1.0_dp, 4.0_dp, 4.0_dp]
    polygon1(2, :) = [4.0_dp, 1.0_dp, 4.0_dp]
    call polygon_collide(4, polygon2, 3, polygon1, collision)
    case_success = (.NOT. collision)
    call print_status(name, case_id, case_success, success)

  end subroutine test_polygon_collide

  subroutine test_solve2x2(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    real(c_double) :: lhs(2, 2), rhs(2), x_val, y_val
    logical(c_bool) :: singular
    character(8) :: name

    case_id = 1
    name = "solve2x2"

    ! CASE 1: Solve without a row-swap.
    lhs(:, 1) = [2.0_dp, 1.0_dp]
    lhs(:, 2) = [3.0_dp, 2.0_dp]
    rhs = [31.0_dp, 19.0_dp]
    call solve2x2(lhs, rhs, singular, x_val, y_val)
    case_success = ( &
         .NOT. singular .AND. &
         x_val == 5.0_dp .AND. &
         y_val == 7.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Solve with a row-swap.
    lhs(:, 1) = [1.0_dp, 4.0_dp]
    lhs(:, 2) = [0.0_dp, 1.0_dp]
    rhs = [3.0_dp, 13.0_dp]
    call solve2x2(lhs, rhs, singular, x_val, y_val)
    case_success = ( &
         .NOT. singular .AND. &
         x_val == 3.0_dp .AND. &
         y_val == 1.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Singular due to first column all zero.
    lhs(:, 1) = 0
    call solve2x2(lhs, rhs, singular, x_val, y_val)
    case_success = singular
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Singular without a row-swap.
    lhs(:, 1) = [2.0_dp, 1.0_dp]
    lhs(:, 2) = [4.0_dp, 2.0_dp]
    call solve2x2(lhs, rhs, singular, x_val, y_val)
    case_success = singular
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Singular with a row-swap.
    lhs(:, 1) = [3.0_dp, 12.0_dp]
    lhs(:, 2) = [1.0_dp, 4.0_dp]
    call solve2x2(lhs, rhs, singular, x_val, y_val)
    case_success = singular
    call print_status(name, case_id, case_success, success)

  end subroutine test_solve2x2

end module test_helpers
