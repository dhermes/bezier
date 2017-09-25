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

  use iso_c_binding, only: c_double, c_bool
  use helpers, only: vector_close
  use types, only: dp
  implicit none
  private
  public test_vector_close

contains

  subroutine test_vector_close(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical(c_bool) :: is_close
    real(c_double) :: eps
    real(c_double) :: vec1(1, 2)
    real(c_double) :: vec2(1, 2)

    eps = 0.5_dp**40

    ! CASE 1: Identical vector.
    vec1(1, :) = [0.5_dp, 4.0_dp]
    is_close = vector_close(2, vec1, vec1, eps)
    if (is_close) then
       write (*, "(A)") "vector_close: Case 1 success"
    else
       write (*, "(A)") "vector_close: Case 1 failure"
       success = .FALSE.
    end if

    ! CASE 2: Far apart vectors.
    vec1(1, :) = [0.0_dp, 6.0_dp]
    vec2(1, :) = [1.0_dp, -4.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (.NOT. is_close) then
       write (*, "(A)") "vector_close: Case 2 success"
    else
       write (*, "(A)") "vector_close: Case 2 failure"
       success = .FALSE.
    end if

    ! CASE 3: Close but different.
    vec1(1, :) = [2.25_dp, -3.5_dp]
    vec2(1, :) = vec1(1, :) + 0.5_dp**43 * [-5.0_dp, 12.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (is_close) then
       write (*, "(A)") "vector_close: Case 3 success"
    else
       write (*, "(A)") "vector_close: Case 3 failure"
       success = .FALSE.
    end if

    ! CASE 4: Custom epsilon.
    vec1(1, :) = [3.0_dp, 4.0_dp]
    vec2(1, :) = [2.0_dp, 5.0_dp]
    is_close = vector_close(2, vec1, vec2, 0.5_dp)
    if (is_close .AND. .NOT. vector_close(2, vec1, vec2, eps)) then
       write (*, "(A)") "vector_close: Case 4 success"
    else
       write (*, "(A)") "vector_close: Case 4 failure"
       success = .FALSE.
    end if

    ! CASE 5: Near zero.
    vec1(1, :) = [0.0_dp, 0.0_dp]
    vec2(1, :) = 0.5_dp**45 * [3.0_dp, 4.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (is_close) then
       write (*, "(A)") "vector_close: Case 5 success"
    else
       write (*, "(A)") "vector_close: Case 5 failure"
       success = .FALSE.
    end if

    ! CASE 6: Near zero failure (i.e. not near enough).
    vec1(1, :) = 0.5_dp**20 * [1.0_dp, 0.0_dp]
    vec2(1, :) = [0.0_dp, 0.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (.NOT. is_close) then
       write (*, "(A)") "vector_close: Case 6 success"
    else
       write (*, "(A)") "vector_close: Case 6 failure"
       success = .FALSE.
    end if

  end subroutine test_vector_close

end module test_helpers
