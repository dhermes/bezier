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
  use helpers, only: cross_product, bbox, vector_close, in_interval
  use types, only: dp
  implicit none
  private test_cross_product, test_vector_close, test_in_interval
  public helpers_all_tests

contains

  subroutine helpers_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_cross_product(success)
    call test_bbox(success)
    call test_vector_close(success)
    call test_in_interval(success)

  end subroutine helpers_all_tests

  subroutine test_cross_product(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: vec0(1, 2)
    real(c_double) :: vec1(1, 2)
    real(c_double) :: result_

    ! CASE 1: Identical vector.
    vec0(1, :) = [1.0_dp, 7.0_dp] / 8
    vec1(1, :) = [-11.0_dp, 24.0_dp] / 32
    call cross_product(vec0, vec1, result_)
    if (result_ == 101.0_dp / 256) then
       write (*, "(A)") "cross_product: Case 1 success"
    else
       write (*, "(A)") "cross_product: Case 1 failure"
       success = .FALSE.
    end if

  end subroutine test_cross_product

  subroutine test_bbox(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes(6, 2)
    real(c_double) :: left, right, bottom, top

    ! CASE 1: Simple case (just two nodes).
    nodes(1, :) = [0.0_dp, 5.0_dp]
    nodes(2, :) = [1.0_dp, 3.0_dp]
    call bbox(2, nodes(:2, :), left, right, bottom, top)
    if (left == 0.0_dp .AND. right == 1.0_dp &
         .AND. bottom == 3.0_dp .AND. top == 5.0_dp) then
       write (*, "(A)") "         bbox: Case 1 success"
    else
       write (*, "(A)") "         bbox: Case 1 failure"
       success = .FALSE.
    end if

    ! CASE 2: "Many" nodes.
    nodes(1, :) = [1.0_dp, 0.0_dp]
    nodes(2, :) = [2.0_dp, 1.0_dp]
    nodes(3, :) = [-1.0_dp, 2.0_dp]
    nodes(4, :) = [5.0_dp, -3.0_dp]
    nodes(5, :) = [4.0_dp, 4.0_dp]
    nodes(6, :) = [0.0_dp, 0.0_dp]
    call bbox(6, nodes, left, right, bottom, top)
    if (left == -1.0_dp .AND. right == 5.0_dp &
         .AND. bottom == -3.0_dp .AND. top == 4.0_dp) then
       write (*, "(A)") "         bbox: Case 2 success"
    else
       write (*, "(A)") "         bbox: Case 2 failure"
       success = .FALSE.
    end if

  end subroutine test_bbox

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
       write (*, "(A)") " vector_close: Case 1 success"
    else
       write (*, "(A)") " vector_close: Case 1 failure"
       success = .FALSE.
    end if

    ! CASE 2: Far apart vectors.
    vec1(1, :) = [0.0_dp, 6.0_dp]
    vec2(1, :) = [1.0_dp, -4.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (.NOT. is_close) then
       write (*, "(A)") " vector_close: Case 2 success"
    else
       write (*, "(A)") " vector_close: Case 2 failure"
       success = .FALSE.
    end if

    ! CASE 3: Close but different.
    vec1(1, :) = [2.25_dp, -3.5_dp]
    vec2(1, :) = vec1(1, :) + 0.5_dp**43 * [-5.0_dp, 12.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (is_close) then
       write (*, "(A)") " vector_close: Case 3 success"
    else
       write (*, "(A)") " vector_close: Case 3 failure"
       success = .FALSE.
    end if

    ! CASE 4: Custom epsilon.
    vec1(1, :) = [3.0_dp, 4.0_dp]
    vec2(1, :) = [2.0_dp, 5.0_dp]
    is_close = vector_close(2, vec1, vec2, 0.5_dp)
    if (is_close .AND. .NOT. vector_close(2, vec1, vec2, eps)) then
       write (*, "(A)") " vector_close: Case 4 success"
    else
       write (*, "(A)") " vector_close: Case 4 failure"
       success = .FALSE.
    end if

    ! CASE 5: Near zero.
    vec1(1, :) = [0.0_dp, 0.0_dp]
    vec2(1, :) = 0.5_dp**45 * [3.0_dp, 4.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (is_close) then
       write (*, "(A)") " vector_close: Case 5 success"
    else
       write (*, "(A)") " vector_close: Case 5 failure"
       success = .FALSE.
    end if

    ! CASE 6: Near zero failure (i.e. not near enough).
    vec1(1, :) = 0.5_dp**20 * [1.0_dp, 0.0_dp]
    vec2(1, :) = [0.0_dp, 0.0_dp]
    is_close = vector_close(2, vec1, vec2, eps)
    if (.NOT. is_close) then
       write (*, "(A)") " vector_close: Case 6 success"
    else
       write (*, "(A)") " vector_close: Case 6 failure"
       success = .FALSE.
    end if

  end subroutine test_vector_close

  subroutine test_in_interval(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical(c_bool) :: is_inside

    ! CASE 1: Interior value.
    is_inside = in_interval(1.5_dp, 1.0_dp, 2.0_dp)
    if (is_inside) then
       write (*, "(A)") "  in_interval: Case 1 success"
    else
       write (*, "(A)") "  in_interval: Case 1 failure"
       success = .FALSE.
    end if

    ! CASE 2: Barely inside.
    is_inside = in_interval(1.0_dp + 0.5_dp**52, 1.0_dp, 2.0_dp)
    if (is_inside) then
       write (*, "(A)") "  in_interval: Case 2 success"
    else
       write (*, "(A)") "  in_interval: Case 2 failure"
       success = .FALSE.
    end if

    ! CASE 3: Barely outside.
    is_inside = in_interval(1.0_dp - 0.5_dp**53, 1.0_dp, 2.0_dp)
    if (.NOT. is_inside) then
       write (*, "(A)") "  in_interval: Case 3 success"
    else
       write (*, "(A)") "  in_interval: Case 3 failure"
       success = .FALSE.
    end if

    ! CASE 4: Exterior value.
    is_inside = in_interval(-1.0_dp, 1.0_dp, 2.0_dp)
    if (.NOT. is_inside) then
       write (*, "(A)") "  in_interval: Case 4 success"
    else
       write (*, "(A)") "  in_interval: Case 4 failure"
       success = .FALSE.
    end if

  end subroutine test_in_interval

end module test_helpers
