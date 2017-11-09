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

program test

  use, intrinsic :: iso_c_binding, only: c_bool
  use test_helpers, only: helpers_all_tests
  use test_curve, only: curve_all_tests
  use test_surface, only: surface_all_tests
  use test_curve_intersection, only: curve_intersection_all_tests
  use test_surface_intersection, only: surface_intersection_all_tests
  implicit none

  logical(c_bool) :: success

  success = .TRUE.
  call helpers_all_tests(success)
  call curve_all_tests(success)
  call surface_all_tests(success)
  call curve_intersection_all_tests(success)
  call surface_intersection_all_tests(success)
  if (.NOT. success) then
     stop 1  ! LCOV_EXCL_LINE
  end if

end program test
