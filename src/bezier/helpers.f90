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

  use types, only: dp
  implicit none
  private
  public cross_product, bbox, wiggle_interval

contains

  subroutine cross_product(vec0, vec1, result_)

    real(dp), intent(in) :: vec0(1, 2)
    real(dp), intent(in) :: vec1(1, 2)
    real(dp), intent(out) :: result_

    result_ = vec0(1, 1) * vec1(1, 2) - vec0(1, 2) * vec1(1, 1)

  end subroutine cross_product

  subroutine bbox(num_nodes, nodes, left, right, bottom, top)

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    integer :: num_nodes
    real(dp), intent(in) :: nodes(num_nodes, 2)
    real(dp), intent(out) :: left, right, bottom, top
    ! Variables outside of signature.
    real(dp) :: workspace(2)

    workspace = minval(nodes, 1)
    left = workspace(1)
    bottom = workspace(2)
    workspace = maxval(nodes, 1)
    right = workspace(1)
    top = workspace(2)

  end subroutine bbox

  subroutine wiggle_interval(value_, result_, success)

    real(dp), intent(in) :: value_
    real(dp), intent(out) :: result_
    logical(1), intent(out) :: success

    success = .TRUE.
    if (-0.5_dp**45 < value_ .AND. value_ < 0.5_dp**45) then
       result_ = 0.0_dp
    else if (0.5_dp**45 <= value_ .AND. value_ <= 1.0_dp - 0.5_dp**45) then
       result_ = value_
    else if ( &
         1.0_dp - 0.5_dp**45 < value_ .AND. value_ < 1.0_dp + 0.5_dp**45) then
       result_ = 1.0_dp
    else
       success = .FALSE.
    end if

  end subroutine wiggle_interval

end module helpers
