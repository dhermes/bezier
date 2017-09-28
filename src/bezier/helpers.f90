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
  private
  public &
       WIGGLE, cross_product, bbox, wiggle_interval, contains_nd, &
       vector_close, in_interval, ulps_away

  real(c_double), parameter :: WIGGLE = 0.5_dp**45

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

end module helpers
