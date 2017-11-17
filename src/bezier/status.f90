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

module status

  use, intrinsic :: iso_c_binding, only: c_int
  implicit none
  private
  public &
       Status_SUCCESS, Status_PARALLEL, Status_WIGGLE_FAIL, &
       Status_NO_CONVERGE, Status_TOO_SMALL, Status_SAME_CURVATURE, &
       Status_BAD_TANGENT, Status_EDGE_END

  ! Values of Status enum:
  integer(c_int), parameter :: Status_SUCCESS = 0
  ! PARALLEL: Corresponds to ``NotImplementedError('Line segments parallel.')``
  integer(c_int), parameter :: Status_PARALLEL = 1
  ! WIGGLE_FAIL: Indicates that ``wiggle_interval`` failed.
  integer(c_int), parameter :: Status_WIGGLE_FAIL = 2
  integer(c_int), parameter :: Status_NO_CONVERGE = 3
  integer(c_int), parameter :: Status_TOO_SMALL = 4
  integer(c_int), parameter :: Status_SAME_CURVATURE = 5
  integer(c_int), parameter :: Status_BAD_TANGENT = 6
  integer(c_int), parameter :: Status_EDGE_END = 7

end module status
