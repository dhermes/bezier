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

program unit_test

  use, intrinsic :: iso_c_binding, only: c_bool
  use curve_intersection, only: free_curve_intersections_workspace
  use triangle_intersection, only: free_triangle_intersections_workspace
  use test_helpers, only: helpers_all_tests
  use test_curve, only: curve_all_tests
  use test_triangle, only: triangle_all_tests
  use test_curve_intersection, only: curve_intersection_all_tests
  use test_triangle_intersection, only: triangle_intersection_all_tests
  implicit none

  logical(c_bool) :: success

  success = .TRUE.
  call helpers_all_tests(success)
  call curve_all_tests(success)
  call triangle_all_tests(success)
  call curve_intersection_all_tests(success)
  call triangle_intersection_all_tests(success)
  ! Wrap up allocated globals.
  call free_curve_intersections_workspace()
  call free_triangle_intersections_workspace()
  if (.NOT. success) then
     stop 1  ! LCOV_EXCL_LINE
  end if

end program unit_test
