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

module helpers

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_bool
  use types, only: dp
  implicit none
  private min_index, sort_in_place, is_separating
  public &
       WIGGLE, VECTOR_CLOSE_EPS, cross_product, bbox, wiggle_interval, &
       contains_nd, vector_close, in_interval, in_sorted, convex_hull, &
       polygon_collide, solve2x2

  real(c_double), parameter :: WIGGLE = 0.5_dp**44
  ! NOTE: This is intended to be used as the default value for ``eps``
  !       in ``vector_close``.
  real(c_double), parameter :: VECTOR_CLOSE_EPS = 0.5_dp**40

contains

  pure subroutine cross_product( &
       vec0, vec1, result_) &
       bind(c, name='BEZ_cross_product')

    real(c_double), intent(in) :: vec0(2)
    real(c_double), intent(in) :: vec1(2)
    real(c_double), intent(out) :: result_

    result_ = vec0(1) * vec1(2) - vec0(2) * vec1(1)

  end subroutine cross_product

  subroutine bbox( &
       num_nodes, nodes, left, right, bottom, top) &
       bind(c, name='BEZ_bbox')

    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: nodes(2, num_nodes)
    real(c_double), intent(out) :: left, right, bottom, top
    ! Variables outside of signature.
    real(c_double) :: workspace(2)

    workspace = minval(nodes, 2)
    left = workspace(1)
    bottom = workspace(2)
    workspace = maxval(nodes, 2)
    right = workspace(1)
    top = workspace(2)

  end subroutine bbox

  subroutine wiggle_interval( &
       value_, result_, success) &
       bind(c, name='BEZ_wiggle_interval')

    real(c_double), intent(in) :: value_
    real(c_double), intent(out) :: result_
    logical(c_bool), intent(out) :: success

    success = .TRUE.
    if (-WIGGLE < value_ .AND. value_ < WIGGLE) then
       result_ = 0.0_dp
    else if (WIGGLE <= value_ .AND. value_ <= 1.0_dp - WIGGLE) then
       result_ = value_
    else if ( &
         1.0_dp - WIGGLE < value_ .AND. value_ < 1.0_dp + WIGGLE) then
       result_ = 1.0_dp
    else
       success = .FALSE.
    end if

  end subroutine wiggle_interval

  pure subroutine contains_nd( &
       num_nodes, dimension_, nodes, point, predicate) &
       bind(c, name='BEZ_contains_nd')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(in) :: point(dimension_)
    logical(c_bool), intent(out) :: predicate

    if (any(point < minval(nodes, 2))) then
       predicate = .FALSE.
    else if (any(maxval(nodes, 2) < point)) then
       predicate = .FALSE.
    else
       predicate = .TRUE.
    end if

  end subroutine contains_nd

  logical(c_bool) pure function vector_close( &
       num_values, vec1, vec2, eps) result(is_close) &
       bind(c, name='BEZ_vector_close')

    integer(c_int), intent(in) :: num_values
    real(c_double), intent(in) :: vec1(num_values)
    real(c_double), intent(in) :: vec2(num_values)
    real(c_double), intent(in) :: eps
    ! Variables outside of signature.
    real(c_double) :: size1, size2

    size1 = norm2(vec1)
    size2 = norm2(vec2)
    if (size1 == 0.0_dp) then
       is_close = size2 <= eps
    else if (size2 == 0.0_dp) then
       is_close = size1 <= eps
    else
       ! This is the most common path.
       is_close = norm2(vec1 - vec2) <= eps * min(size1, size2)
    end if

  end function vector_close

  pure function in_interval( &
       value_, start, end) result(predicate) &
       bind(c, name='BEZ_in_interval')

    real(c_double), intent(in) :: value_, start, end
    logical(c_bool) :: predicate

    predicate = (start <= value_) .AND. (value_ <= end)

  end function in_interval

  subroutine min_index(num_points, points, match)

    ! Finds the minimum point when sorting first by ``x`` then by
    ! ``y`` coordinate (in the case of a tie in ``x``).

    ! NOTE: This is a helper for ``sort_in_place``.
    ! NOTE: This assumes, but does not check, that ``num_points > 0``.

    integer(c_int), intent(in) :: num_points
    real(c_double), intent(inout) :: points(2, num_points)
    integer(c_int), intent(out) :: match
    ! Variables outside of signature.
    integer(c_int) :: i

    match = 1
    do i = 2, num_points
       if (points(1, i) < points(1, match)) then
          match = i
       else if (points(1, i) == points(1, match)) then
          if (points(2, i) < points(2, match)) then
             match = i
          end if
       end if
    end do

  end subroutine min_index

  subroutine sort_in_place(num_points, points, num_uniques)

    ! NOTE: This is a helper for ``convex_hull``. This does a simple quadratic
    !       sort (i.e. it is suboptimal) because it expects the number of
    !       points to be small. It **also** makes sure the values are unique
    !       to exact precision.

    integer(c_int), intent(in) :: num_points
    real(c_double), intent(inout) :: points(2, num_points)
    integer(c_int), intent(out) :: num_uniques
    ! Variables outside of signature.
    integer(c_int) :: i, match
    real(c_double) :: swap(2)

    num_uniques = num_points
    ! Place the "smallest" point first (since it can't be a duplicate).
    call min_index(num_points, points, match)
    if (match /= 1) then
       swap = points(:, match)
       points(:, match) = points(:, 1)
       points(:, 1) = swap
    end if

    i = 2
    do while (i <= num_uniques)
       call min_index(num_uniques + 1 - i, points(:, i:num_uniques), match)
       ! Shift ``match`` based on the slice ``points(:, i:...)``.
       match = match + i - 1
       if (match /= i) then
          swap = points(:, match)
          if (all(swap == points(:, i - 1))) then
             ! This means ``match`` is a duplicate.
             points(:, match) = points(:, num_uniques)
             points(:, num_uniques) = swap
             num_uniques = num_uniques - 1
          else
             points(:, match) = points(:, i)
             points(:, i) = swap
          end if
       end if
       ! Increment for next iteration.
       i = i + 1
    end do

  end subroutine sort_in_place

  logical(c_bool) pure function in_sorted( &
       num_values, values, value_) result(contained)

    ! NOTE: This assumes ``num_values >= 1``.
    ! NOTE: This subroutine is not part of the C ABI for this module,
    !       but it is (for now) public, so that it can be tested.

    integer(c_int), intent(in) :: num_values
    integer(c_int), intent(in) :: values(num_values)
    integer(c_int), intent(in) :: value_
    ! Variables outside of signature.
    integer(c_int) :: left, right, midpoint

    left = 1
    right = num_values

    do while (left < right)
       midpoint = (left + right) / 2
       if (value_ == values(midpoint)) then
          contained = .TRUE.
          return
       else if (value_ < values(midpoint)) then
          right = midpoint - 1
       else
          left = midpoint + 1
       end if
    end do

    contained = values(left) == value_

  end function in_sorted

  subroutine convex_hull( &
       num_points, points, polygon_size, polygon) &
       bind(c, name='BEZ_simple_convex_hull')

    ! NOTE: This uses Andrew's monotone chain convex hull algorithm and used
    !       (https://en.wikibooks.org/wiki/Algorithm_Implementation/
    !        Geometry/Convex_hull/Monotone_chain)
    !       as motivation. The code there is licensed CC BY-SA 3.0.

    ! NOTE: The "standard" implementation for this is in Qhull, but for now
    !       that is avoided as an extra dependency.

    integer(c_int), intent(in) :: num_points
    real(c_double), intent(in) :: points(2, num_points)
    integer(c_int), intent(out) :: polygon_size
    real(c_double), intent(out) :: polygon(2, num_points)
    ! Variables outside of signature.
    real(c_double) :: point1(2), point2(2), point3(2)
    integer(c_int) :: num_uniques
    real(c_double) :: uniques(2, num_points)
    integer(c_int) :: num_lower, num_upper
    integer(c_int) :: lower(num_points), upper(num_points)
    integer(c_int) :: i
    real(c_double) :: result_

    uniques = points
    call sort_in_place(num_points, uniques, num_uniques)

    ! In the "boring" case that we have fewer than 2 points, return.
    if (num_uniques == 0) then
       polygon_size = 0
       return
    else if (num_uniques < 3) then
       polygon_size = num_uniques
       polygon(:, :num_uniques) = uniques(:, :num_uniques)
       return
    end if

    ! First create a "lower" convex hull
    num_lower = 2  ! We know ``num_uniques >= 3``.
    lower(:2) = [1, 2]
    do i = 3, num_uniques
       point3 = uniques(:, i)
       result_ = -1.0_dp  ! Dummy value that is ``<= 0.0``.
       do while (num_lower > 1 .AND. result_ <= 0.0_dp)
          point1 = uniques(:, lower(num_lower - 1))
          point2 = uniques(:, lower(num_lower))
          ! If ``point3`` (the one we are considering) is more "inside"
          ! of ``point1`` than ``point2`` is, then we drop ``point2``.
          call cross_product(point2 - point1, point3 - point1, result_)
          if (result_ <= 0.0_dp) then
             num_lower = num_lower - 1
          end if
       end do
       num_lower = num_lower + 1
       lower(num_lower) = i
    end do

    ! Then create an "upper" convex hull
    num_upper = 1
    upper(1) = num_uniques
    upper_loop: do i = num_uniques - 1, 1, -1
       ! Don't consider indices from the lower hull (other than the ends).
       if (i > 1 .AND. in_sorted(num_lower, lower(:num_lower), i)) then
          cycle upper_loop
       end if
       point3 = uniques(:, i)
       result_ = -1.0_dp  ! Dummy value that is ``<= 0.0``.
       do while (num_upper > 1 .AND. result_ <= 0.0_dp)
          point1 = uniques(:, upper(num_upper - 1))
          point2 = uniques(:, upper(num_upper))
          ! If ``point3`` (the one we are considering) is more "inside"
          ! of ``point1`` than ``point2`` is, then we drop ``point2``.
          call cross_product(point2 - point1, point3 - point1, result_)
          if (result_ <= 0.0_dp) then
             num_upper = num_upper - 1
          end if
       end do
       num_upper = num_upper + 1
       upper(num_upper) = i
    end do upper_loop

    ! The endpoints are **both** double counted, so we skip the "end"
    ! of each hull (lower and upper).
    polygon_size = 0
    do i = 1, num_lower - 1
       polygon_size = polygon_size + 1
       polygon(:, polygon_size) = uniques(:, lower(i))
    end do
    do i = 1, num_upper - 1
       polygon_size = polygon_size + 1
       polygon(:, polygon_size) = uniques(:, upper(i))
    end do

  end subroutine convex_hull

  logical(c_bool) pure function is_separating( &
       edge_direction, polygon_size1, polygon1, polygon_size2, polygon2) &
       result(predicate)

    ! NOTE: This assumes, but does not check, that ``polygon_sizeX`` is
    !       at least two.
    real(c_double), intent(in) :: edge_direction(2)
    integer(c_int), intent(in) :: polygon_size1
    real(c_double), intent(in) :: polygon1(2, polygon_size1)
    integer(c_int), intent(in) :: polygon_size2
    real(c_double), intent(in) :: polygon2(2, polygon_size2)
    ! Variables outside of signature.
    real(c_double) :: vertex(2), norm_squared, cp_result, param
    real(c_double) :: min_param1, max_param1
    real(c_double) :: min_param2, max_param2
    integer(c_int) :: i

    ! NOTE: We assume throughout that ``norm_squared != 0``. If it **were**
    !       zero that would mean the ``edge_direction`` corresponds to an
    !       invalid edge.
    norm_squared = dot_product(edge_direction, edge_direction)

    ! Compute the parameters along the "separating axis" for the vertices
    ! of the first polygon.
    vertex = polygon1(:, 1)
    call cross_product(edge_direction, vertex, cp_result)
    param = cp_result / norm_squared
    min_param1 = param
    max_param1 = param
    do i = 2, polygon_size1
       vertex = polygon1(:, i)
       call cross_product(edge_direction, vertex, cp_result)
       param = cp_result / norm_squared
       min_param1 = min(min_param1, param)
       max_param1 = max(max_param1, param)
    end do

    ! Compute the parameters along the "separating axis" for the vertices
    ! of the first polygon.
    vertex = polygon2(:, 1)
    call cross_product(edge_direction, vertex, cp_result)
    param = cp_result / norm_squared
    min_param2 = param
    max_param2 = param
    do i = 2, polygon_size2
       vertex = polygon2(:, i)
       call cross_product(edge_direction, vertex, cp_result)
       param = cp_result / norm_squared
       min_param2 = min(min_param2, param)
       max_param2 = max(max_param2, param)
    end do

    ! Determine if [min_param1, max_param1] and [min_param2, max_param2]
    ! are disjoint intervals.
    predicate = (min_param1 > max_param2 .OR. max_param1 < min_param2)

  end function is_separating

  subroutine polygon_collide( &
       polygon_size1, polygon1, polygon_size2, polygon2, collision) &
       bind(c, name='BEZ_polygon_collide')

    ! This determines if two **convex** polygons collide. The polygons
    ! are given as ``2 x N`` arrays of ``x-y`` points (one per column)
    ! in order they appear (the first and last node are **not** the same
    ! as is sometimes expected for closed polygons).

    ! This code uses the Separating axis theorem (SAT) [1] [2] to quickly
    ! determine if the polygons intersect.
    ! [1]: https://en.wikipedia.org/wiki/Hyperplane_separation_theorem
    ! [2]: https://hackmd.io/s/ryFmIZrsl

    ! A "standard" implementation for this is in GEOS [3], but for now
    ! that is avoided as an extra dependency.
    ! [3]: https://trac.osgeo.org/geos

    integer(c_int), intent(in) :: polygon_size1
    real(c_double), intent(in) :: polygon1(2, polygon_size1)
    integer(c_int), intent(in) :: polygon_size2
    real(c_double), intent(in) :: polygon2(2, polygon_size2)
    logical(c_bool), intent(out) :: collision
    ! Variables outside of signature.
    real(c_double) :: edge_direction(2)
    integer(c_int) :: i

    collision = .TRUE.

    ! First handle the "wrap-around" edge from polygon1.
    edge_direction = polygon1(:, 1) - polygon1(:, polygon_size1)
    if (is_separating( &
         edge_direction, polygon_size1, polygon1, &
         polygon_size2, polygon2)) then
       collision = .FALSE.
       return
    end if

    ! Then, check all other edges from polygon1.
    do i = 2, polygon_size1
       edge_direction = polygon1(:, i) - polygon1(:, i - 1)
       if (is_separating( &
            edge_direction, polygon_size1, polygon1, &
            polygon_size2, polygon2)) then
          collision = .FALSE.
          return
       end if
    end do

    ! Then, handle the "wrap-around" edge from polygon2.
    edge_direction = polygon2(:, 1) - polygon2(:, polygon_size2)
    if (is_separating( &
         edge_direction, polygon_size1, polygon1, &
         polygon_size2, polygon2)) then
       collision = .FALSE.
       return
    end if

    ! Then, check all other edges from polygon1.
    do i = 2, polygon_size2
       edge_direction = polygon2(:, i) - polygon2(:, i - 1)
       if (is_separating( &
            edge_direction, polygon_size1, polygon1, &
            polygon_size2, polygon2)) then
          collision = .FALSE.
          return
       end if
    end do

  end subroutine polygon_collide

  subroutine solve2x2(lhs, rhs, singular, x_val, y_val)

    real(c_double), intent(in) :: lhs(2, 2)
    real(c_double), intent(in) :: rhs(2)
    logical(c_bool), intent(out) :: singular
    real(c_double), intent(out) :: x_val
    real(c_double), intent(out) :: y_val
    ! Variables outside of signature.
    real(c_double) :: ratio, denominator

    ! A <--> lhs(1, 1)
    ! B <--> lhs(1, 2)
    ! C <--> lhs(2, 1)
    ! D <--> lhs(2, 2)
    ! E <--> rhs(1)
    ! F <--> rhs(2)
    if (abs(lhs(2, 1)) > abs(lhs(1, 1))) then
       ! NOTE: We know there is no division by zero here since ``C``
       !       is **strictly** bigger than **some** value (in magnitude).
       ! [A | B][x] = [E]
       ! [C | D][y]   [F]
       ratio = lhs(1, 1) / lhs(2, 1)
       ! r = A / C
       ! [A - rC | B - rD][x]   [E - rF]
       ! [C      | D     ][y] = [F     ]
       ! ==> 0x + (B - rD) y = E - rF
       denominator = lhs(1, 2) - ratio * lhs(2, 2)
       if (denominator == 0.0_dp) then
          singular = .TRUE.
          return
       end if

       y_val = (rhs(1) - ratio * rhs(2)) / denominator
       ! Cx + Dy = F ==> x = (F - Dy) / C
       x_val = (rhs(2) - lhs(2, 2) * y_val) / lhs(2, 1)
       singular = .FALSE.
    else
       if (lhs(1, 1) == 0.0_dp) then
          singular = .TRUE.
          return
       end if
       ! [A | B][x] = [E]
       ! [C | D][y]   [F]
       ratio = lhs(2, 1) / lhs(1, 1)
       ! r = C / A
       ! [A      | B     ][x] = [E     ]
       ! [C - rA | D - rB][y]   [F - rE]
       ! ==> 0x + (D - rB) y = F - rE
       denominator = lhs(2, 2) - ratio * lhs(1, 2)
       if (denominator == 0.0_dp) then
          singular = .TRUE.
          return
       end if

       y_val = (rhs(2) - ratio * rhs(1)) / denominator
       ! Ax + By = E ==> x = (E - B y) / A
       x_val = (rhs(1) - lhs(1, 2) * y_val) / lhs(1, 1)
       singular = .FALSE.
    end if

  end subroutine solve2x2

end module helpers
