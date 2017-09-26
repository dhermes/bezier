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

module unit_test_helpers

  implicit none
  private
  public print_status

contains

  subroutine print_status(name, case_id, success)
    character(len=*), intent(in) :: name
    integer, intent(inout) :: case_id
    logical, intent(in) :: success
    ! Variables outside of signature.
    character(len=8) :: status_msg

    if (success) then
       status_msg = " success"
    else
       status_msg = " failure"
    end if

    write (*, "(A15, A, I2, A)") name, ": Case ", case_id, status_msg

    ! Increment case ID for next case.
    case_id = case_id + 1

  end subroutine print_status

end module unit_test_helpers
