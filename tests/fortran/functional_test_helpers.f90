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

module functional_test_helpers

  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use types, only: dp
  use status, only: &
       Status_SUCCESS, Status_PARALLEL, Status_NO_CONVERGE, &
       Status_SAME_CURVATURE, Status_BAD_TANGENT, Status_EDGE_END, &
       Status_UNKNOWN
  use curve_intersection, only: &
       all_intersections, free_curve_intersections_workspace
  implicit none
  private CURVE_INTERSECTIONS, max_ulp, show_curve_intersect_result, sort_vals
  public intersect_and_check, clean_up

  real(c_double), allocatable :: CURVE_INTERSECTIONS(:, :)

contains

  integer(c_int) function max_ulp(m, n, guess, known) result(ulp_err)
    integer(c_int), intent(in) :: m, n
    real(c_double), intent(in) :: guess(m, n)
    real(c_double), intent(in) :: known(m, n)
    ! Variables outside of signature.
    real(c_double) :: rel_errs(m, n)

    if (m == 0 .OR. n == 0) then
       ulp_err = 0
       return
    end if

    rel_errs = guess - known
    where (known == 0.0_dp)
       where (guess == 0.0_dp)
          rel_errs = 0.0_dp
       elsewhere (abs(guess) < 1.0_dp)
          rel_errs = exponent(guess) - 1
       elsewhere
          rel_errs = HUGE(ulp_err)
       end where
    elsewhere
       rel_errs = abs(rel_errs / spacing(known))
    end where

    ulp_err = minval(rel_errs)
    if (ulp_err >= 0) then
       ulp_err = maxval(rel_errs)
    end if

  end function max_ulp

  subroutine show_curve_intersect_result( &
       case_id, status, actual_n, expected_n, known)

    integer(c_int), intent(in) :: case_id, status, actual_n, expected_n
    real(c_double), intent(in) :: known(2, expected_n)
    ! Variables outside of signature.
    integer(c_int) :: ulp_err

    if (status == Status_SUCCESS) then
       if (actual_n == expected_n) then
          ulp_err = max_ulp( &
               2, expected_n, CURVE_INTERSECTIONS(:, :expected_n), known)
          write (*, '(A, I2, A, I0, A, I0)') &
               "Case ", &
               case_id, &
               " (success): num_intersections = ", &
               expected_n, &
               ", Max ULP error = ", &
               ulp_err
       else
          write (*, '(A, I2, A, I0, A, I0, A)') &
               "Case ", &
               case_id, &
               " (failure): ", &
               actual_n, &
               " intersections (expected ", &
               expected_n, &
               ")"
       end if
    else if (status == Status_PARALLEL) then
       write (*, '(A, I2, A)') &
            "Case ", &
            case_id, &
            " (failure): PARALLEL"
    else if (status == Status_NO_CONVERGE) then
       write (*, '(A, I2, A)') &
            "Case ", &
            case_id, &
            " (failure): NO_CONVERGE"
    else if (status == Status_SAME_CURVATURE) then
       write (*, '(A, I2, A)') &
            "Case ", &
            case_id, &
            " (failure): SAME_CURVATURE"
    else if (status == Status_BAD_TANGENT) then
       write (*, '(A, I2, A)') &
            "Case ", &
            case_id, &
            " (failure): BAD_TANGENT"
    else if (status == Status_EDGE_END) then
       write (*, '(A, I2, A)') &
            "Case ", &
            case_id, &
            " (failure): EDGE_END"
    else if (status == Status_UNKNOWN) then
       write (*, '(A, I2, A)') &
            "Case ", &
            case_id, &
            " (failure): UNKNOWN"
    else
       write (*, '(A, I2, A, I0, A)') &
            "Case ", &
            case_id, &
            " (failure): ", &
            status, &
            " intersection candidates"
    end if

  end subroutine show_curve_intersect_result

  subroutine sort_vals(n, status)
    ! Sort by the first row.

    integer(c_int), intent(in) :: n
    integer(c_int), intent(in) :: status
    ! Variables outside of signature.
    real(c_double) :: swap(2)
    integer(c_int) :: i, w

    if (status /= Status_SUCCESS) then
       return
    end if

    do i = 1, n
       w = minloc(CURVE_INTERSECTIONS(1, i:n), 1)
       if (w /= 1) then
          swap = CURVE_INTERSECTIONS(:, i + w - 1)
          CURVE_INTERSECTIONS(:, i + w - 1) = CURVE_INTERSECTIONS(:, i)
          CURVE_INTERSECTIONS(:, i) = swap
       end if
    end do

  end subroutine sort_vals

  subroutine intersect_and_check( &
       case_id, num_nodes1, nodes1, num_nodes2, nodes2, &
       expected_n, expected_params)

    integer(c_int), intent(in) :: case_id
    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(num_nodes1, 2)
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(num_nodes2, 2)
    integer(c_int), intent(in) :: expected_n
    real(c_double), intent(in) :: expected_params(2, expected_n)
    ! Variables outside of signature.
    integer(c_int) :: num_intersections, status

    call all_intersections( &
         num_nodes1, nodes1, &
         num_nodes2, nodes2, &
         CURVE_INTERSECTIONS, num_intersections, &
         status)
    call sort_vals(num_intersections, status)
    call show_curve_intersect_result( &
         case_id, status, num_intersections, expected_n, expected_params)

  end subroutine intersect_and_check

  subroutine clean_up()
    if (allocated(CURVE_INTERSECTIONS)) then
       deallocate(CURVE_INTERSECTIONS)
    end if
    call free_curve_intersections_workspace()
  end subroutine clean_up

end module functional_test_helpers
