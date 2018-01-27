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
       Status_SUCCESS, Status_PARALLEL, Status_NO_CONVERGE, &
       Status_INSUFFICIENT_SPACE, Status_SAME_CURVATURE, Status_BAD_TANGENT, &
       Status_EDGE_END, Status_UNKNOWN

  ! Values of Status enum:
  ! SUCCESS: Procedure exited with no error.
  integer(c_int), parameter :: Status_SUCCESS = 0
  ! PARALLEL: Corresponds to ``NotImplementedError('Line segments parallel.')``
  integer(c_int), parameter :: Status_PARALLEL = 1
  ! NO_CONVERGE: An iterative algorithm has failed to converge. Used by
  !              ``curve_intersection.all_intersections()``.
  integer(c_int), parameter :: Status_NO_CONVERGE = 2
  ! INSUFFICIENT_SPACE: Intended to be used by ABI versions of procedures. Will
  !                     be used when the caller has not allocated enough space
  !                     for the output values. (On the Fortran side, an
  !                     ``allocatable`` array can be resized, but in the ABI,
  !                     we can only use what the caller has passed.)
  integer(c_int), parameter :: Status_INSUFFICIENT_SPACE = 3
  ! SAME_CURVATURE: Can't classify curve-curve intersection when the curves
  !                 have identical curvature at the point of intersection.
  integer(c_int), parameter :: Status_SAME_CURVATURE = 4
  ! BAD_TANGENT: Two curves are tangent where they intersect, but they move
  !              in an opposite direction while defining overlapping arcs.
  integer(c_int), parameter :: Status_BAD_TANGENT = 5
  ! EDGE_END: A surface-surface intersection point occurs at the **end** of
  !           an edge (only intersections at the beginning of an edge should
  !           be used).
  integer(c_int), parameter :: Status_EDGE_END = 6
  ! UNKNOWN: Signifies a block of code reached an "impossible" state. Either
  !          the code has violated some mathematical invariant or the author
  !          misunderstood the possible states of the system.
  integer(c_int), parameter :: Status_UNKNOWN = 999

end module status
