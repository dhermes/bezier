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
  private min_index, sort_in_place
  public &
       WIGGLE, VECTOR_CLOSE_EPS, cross_product, bbox, wiggle_interval, &
       contains_nd, vector_close, in_interval, ulps_away, convex_hull

  real(c_double), parameter :: WIGGLE = 0.5_dp**45
  ! NOTE: This is intended to be used as the default value for ``eps``
  !       in ``vector_close``.
  real(c_double), parameter :: VECTOR_CLOSE_EPS = 0.5_dp**40

contains

  subroutine cross_product( &
       vec0, vec1, result_) &
       bind(c, name='cross_product')

    real(c_double), intent(in) :: vec0(1, 2)
    real(c_double), intent(in) :: vec1(1, 2)
    real(c_double), intent(out) :: result_

    result_ = vec0(1, 1) * vec1(1, 2) - vec0(1, 2) * vec1(1, 1)

  end subroutine cross_product

  subroutine bbox( &
       num_nodes, nodes, left, right, bottom, top) &
       bind(c, name='bbox')

    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: nodes(num_nodes, 2)
    real(c_double), intent(out) :: left, right, bottom, top
    ! Variables outside of signature.
    real(c_double) :: workspace(2)

    workspace = minval(nodes, 1)
    left = workspace(1)
    bottom = workspace(2)
    workspace = maxval(nodes, 1)
    right = workspace(1)
    top = workspace(2)

  end subroutine bbox

  subroutine wiggle_interval( &
       value_, result_, success) &
       bind(c, name='wiggle_interval')

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
       bind(c, name='contains_nd')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    real(c_double), intent(in) :: point(dimension_)
    logical(c_bool), intent(out) :: predicate

    if (any(point < minval(nodes, 1))) then
       predicate = .FALSE.
    else if (any(maxval(nodes, 1) < point)) then
       predicate = .FALSE.
    else
       predicate = .TRUE.
    end if

  end subroutine contains_nd

  logical(c_bool) pure function vector_close( &
       num_values, vec1, vec2, eps) result(is_close) &
       bind(c, name='vector_close')

    integer(c_int), intent(in) :: num_values
    real(c_double), intent(in) :: vec1(1, num_values)
    real(c_double), intent(in) :: vec2(1, num_values)
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
       bind(c, name='in_interval')

    real(c_double), intent(in) :: value_, start, end
    logical(c_bool) :: predicate

    predicate = (start <= value_) .AND. (value_ <= end)

  end function in_interval

  pure function ulps_away( &
       value1, value2, num_bits, eps) result(predicate) &
       bind(c, name='ulps_away')

    real(c_double), intent(in) :: value1, value2
    integer(c_int), intent(in) :: num_bits
    real(c_double), intent(in) :: eps
    logical(c_bool) :: predicate

    if (value1 == 0.0_dp) then
       predicate = abs(value2) < eps
    else if (value2 == 0.0_dp) then
       predicate = abs(value1) < eps
    else
       ! NOTE: This assumes `spacing` always returns a positive (which is
       !       not the case with ``numpy.spacing``, which keeps the sign
       !       of the input).
       ! NOTE: This will be the most common path through this function.
       predicate = abs(value1 - value2) <= num_bits * spacing(value1)
    end if

  end function ulps_away

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

  subroutine sort_in_place(num_points, points)

    ! NOTE: This is a helper for ``convex_hull``. This does a simple quadratic
    !       sort (i.e. it is suboptimal) because it expects the number of
    !       points to be small.

    integer(c_int), intent(in) :: num_points
    real(c_double), intent(inout) :: points(2, num_points)
    ! Variables outside of signature.
    integer(c_int) :: i, match
    real(c_double) :: swap(2)

    do i = 1, num_points - 1
       call min_index(num_points + 1 - i, points(:, i:), match)
       ! Shift ``match`` based on the slice ``points(:, i:)``.
       match = match + i - 1
       if (match /= i) then
          swap = points(:, i)
          points(:, i) = points(:, match)
          points(:, match) = swap
       end if
    end do

  end subroutine sort_in_place

  subroutine convex_hull(num_points, points, polygon_size, polygon)

    ! NOTE: This uses Andrew's monotone chain convex hull algorithm and used
    !       (https://en.wikibooks.org/wiki/Algorithm_Implementation/
    !        Geometry/Convex_hull/Monotone_chain)
    !       as motivation. The code there is licensed CC BY-SA 3.0.

    integer(c_int), intent(in) :: num_points
    real(c_double), intent(in) :: points(2, num_points)
    integer(c_int), intent(out) :: polygon_size
    real(c_double), intent(out) :: polygon(2, num_points)
    ! Variables outside of signature.
    real(c_double) :: point1(1, 2), point2(1, 2), point3(1, 2)
    integer(c_int) :: num_uniques
    real(c_double) :: uniques(2, num_points)
    integer(c_int) :: num_lower, num_upper
    real(c_double) :: lower(2, num_points), upper(2, num_points)
    integer(c_int) :: i, j
    real(c_double) :: result_

    ! First make sure all the points are unique to exact precision.
    num_uniques = 0
    unique_loop: do i = 1, num_points
       point1(1, :) = points(:, i)
       do j = 1, num_uniques
          if (all(point1(1, :) == uniques(:, j))) then
             cycle unique_loop
          end if
       end do
       ! If we haven't cycled, then we know it is unique.
       num_uniques = num_uniques + 1
       uniques(:, num_uniques) = points(:, i)
    end do unique_loop

    ! In the "boring" case that we have fewer than 2 points, return.
    if (num_uniques == 0) then
       polygon_size = 0
       return
    else if (num_uniques == 1) then
       polygon_size = 1
       polygon(:, 1) = uniques(:, 1)
       return
    end if

    call sort_in_place(num_uniques, uniques)

    ! First create a "lower" convex hull
    num_lower = 0
    do i = 1, num_uniques
       point3(1, :) = uniques(:, i)
       result_ = -1.0_dp  ! Dummy value that is ``<= 0.0``.
       do while (num_lower > 1 .AND. result_ <= 0.0_dp)
          point1(1, :) = lower(:, num_lower - 1)
          point2(1, :) = lower(:, num_lower)
          ! If ``point3`` (the one we are considering) is more "inside"
          ! of ``point1`` than ``point2`` is, then we drop ``point2``.
          call cross_product(point2 - point1, point3 - point1, result_)
          if (result_ <= 0.0_dp) then
             num_lower = num_lower - 1
          end if
       end do
       num_lower = num_lower + 1
       lower(:, num_lower) = point3(1, :)
    end do

    ! Then create an "upper" convex hull
    num_upper = 0
    do i = num_uniques, 1, -1
       point3(1, :) = uniques(:, i)
       result_ = -1.0_dp  ! Dummy value that is ``<= 0.0``.
       do while (num_upper > 1 .AND. result_ <= 0.0_dp)
          point1(1, :) = upper(:, num_upper - 1)
          point2(1, :) = upper(:, num_upper)
          ! If ``point3`` (the one we are considering) is more "inside"
          ! of ``point1`` than ``point2`` is, then we drop ``point2``.
          call cross_product(point2 - point1, point3 - point1, result_)
          if (result_ <= 0.0_dp) then
             num_upper = num_upper - 1
          end if
       end do
       num_upper = num_upper + 1
       upper(:, num_upper) = point3(1, :)
    end do

    ! The endpoints are **both** double counted, so we skip the "end"
    ! of each hull (lower and upper).
    polygon_size = 0
    do i = 1, num_lower - 1
       polygon_size = polygon_size + 1
       polygon(:, polygon_size) = lower(:, i)
    end do
    do i = 1, num_upper - 1
       polygon_size = polygon_size + 1
       polygon(:, polygon_size) = upper(:, i)
    end do

  end subroutine convex_hull

end module helpers
