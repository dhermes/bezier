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

  use, intrinsic :: iso_c_binding, only: c_double
  use types, only: dp
  implicit none
  private
  public MACHINE_EPS, print_status, get_random_nodes, binary_round, get_id_mat

  ! NOTE: Should probably use ``d1mach`` to determine this.
  real(c_double), parameter :: MACHINE_EPS = 0.5_dp**52

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

    write (*, "(A26, A, I2, A)") name, ": Case ", case_id, status_msg

    ! Increment case ID for next case.
    case_id = case_id + 1

  end subroutine print_status

  subroutine get_random_nodes(nodes, multiplier, modulus, num_bits)
    real(c_double), intent(inout) :: nodes(:, :)
    integer, intent(in) :: multiplier, modulus
    integer, optional, intent(in) :: num_bits
    ! Variables outside of signature.
    integer, allocatable :: old_seed(:)
    integer, allocatable :: new_seed(:)
    integer :: size_, i

    ! Get the size of the random seed.
    call random_seed(size=size_)
    ! Store the old seed before replacing it.
    allocate(old_seed(size_))
    call random_seed(get=old_seed)

    ! Set a temporary seed based on ``multiplier`` and ``modulus``.
    new_seed = [(mod(multiplier * i, modulus), i = 0, size_ - 1)]
    call random_seed(put=new_seed)
    ! Populate the matrix with random values from our seed.
    call random_number(nodes)

    ! If ``num_bits`` is specified, then binary round the ``nodes``.
    if (present(num_bits)) then
       call binary_round(nodes, num_bits)
    end if

    ! Put back the original seed.
    call random_seed(put=old_seed)

  end subroutine get_random_nodes

  subroutine binary_round(value_, num_bits)
    real(c_double), intent(inout) :: value_(:, :)
    integer, intent(in) :: num_bits
    ! Variables outside of signature.
    real(c_double), allocatable :: work(:, :)

    ! E.g. The first 4 bits of -6691/512 = -13.068359375 = -0x1.a23p+3
    ! are -0x1.ap+3 = -13. To get this, we use `fraction` to compute the
    ! unsigned significand (i.e. `1 + mantissa`) 6691/4096.
    work = 2 * fraction(abs(value_))
    ! Then we multiply by 2^(num_bits) and truncate to an integer
    ! 26 = 0x1.ap+4.
    work = floor(work * 2.0_dp**num_bits)
    ! Then we remove the added exponents (1 from doubling the output
    ! of `fraction` and `num_bits` from the scaling).
    work = work * 0.5_dp**(1 + num_bits)
    ! Finally we restore the original exponent and sign.
    value_ = sign(work * 2.0_dp**exponent(value_), value_)

  end subroutine binary_round

  function get_id_mat(n) result(id_mat)
    integer, intent(in) :: n
    real(c_double) :: id_mat(n, n)
    ! Variables outside of signature.
    integer :: i

    ! Populate both as the identity matrix.
    id_mat = 0
    forall (i = 1:n)
       id_mat(i, i) = 1
    end forall

  end function get_id_mat

end module unit_test_helpers
