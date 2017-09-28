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

module test_surface

  use, intrinsic :: iso_c_binding, only: c_bool, c_double
  use surface, only: &
       de_casteljau_one_round
  use types, only: dp
  use unit_test_helpers, only: print_status, get_random_nodes
  implicit none
  private &
       test_de_casteljau_one_round
  public surface_all_tests

contains

  subroutine surface_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_de_casteljau_one_round(success)

  end subroutine surface_all_tests

  subroutine test_de_casteljau_one_round(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes1(3, 2), nodes2(6, 2), nodes3(10, 2)
    real(c_double) :: new_nodes1(1, 2), new_nodes2(3, 2), new_nodes3(6, 2)
    real(c_double) :: expected1(1, 2), expected2(3, 2), expected3(6, 2)
    real(c_double) :: lambda1, lambda2, lambda3
    integer :: case_id
    character(:), allocatable :: name

    case_id = 1
    name = "de_casteljau_one_round"

    ! CASE 1: Linear reference triangle (just goes to point).
    nodes1(1, :) = 0
    nodes1(2, :) = [1.0_dp, 0.0_dp]
    nodes1(3, :) = [0.0_dp, 1.0_dp]
    lambda1 = 0.125_dp
    lambda2 = 0.5_dp
    lambda3 = 0.375_dp
    expected1(1, :) = [lambda2, lambda3]
    call de_casteljau_one_round( &
         3, 2, nodes1, 1, &
         lambda1, lambda2, lambda3, new_nodes1)
    if (all(new_nodes1 == expected1)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 2: Quadratic surface.
    call get_random_nodes(nodes2, 790931, 1483, num_bits=8)
    ! NOTE: Use a fixed seed so the test is deterministic and round
    !       the nodes to 8 bits of precision to avoid round-off.
    lambda1 = 0.625_dp
    lambda2 = 0.25_dp
    lambda3 = 0.125_dp
    expected2(1, :) = ( &
         lambda1 * nodes2(1, :) + &  ! p200
         lambda2 * nodes2(2, :) + &  ! p110
         lambda3 * nodes2(4, :))  ! p101
    expected2(2, :) = ( &
         lambda1 * nodes2(2, :) + &  ! p110
         lambda2 * nodes2(3, :) + &  ! p020
         lambda3 * nodes2(5, :))  ! p011
    expected2(3, :) = ( &
         lambda1 * nodes2(4, :) + &  ! p101
         lambda2 * nodes2(5, :) + &  ! p011
         lambda3 * nodes2(6, :))  ! p002
    call de_casteljau_one_round( &
         6, 2, nodes2, 2, &
         lambda1, lambda2, lambda3, new_nodes2)
    if (all(new_nodes2 == expected2)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

    ! CASE 3: Cubic surface.
    nodes3(1, :) = [0.0_dp, 0.0_dp]
    nodes3(2, :) = [3.25_dp, 1.5_dp]
    nodes3(3, :) = [6.5_dp, 1.5_dp]
    nodes3(4, :) = [10.0_dp, 0.0_dp]
    nodes3(5, :) = [1.5_dp, 3.25_dp]
    nodes3(6, :) = [5.0_dp, 5.0_dp]
    nodes3(7, :) = [10.0_dp, 5.25_dp]
    nodes3(8, :) = [1.5_dp, 6.5_dp]
    nodes3(9, :) = [5.25_dp, 10.0_dp]
    nodes3(10, :) = [0.0_dp, 10.0_dp]
    lambda1 = 0.375_dp
    lambda2 = 0.25_dp
    lambda3 = 0.375_dp
    expected3(1, :) = ( &
         lambda1 * nodes3(1, :) + &  ! p300
         lambda2 * nodes3(2, :) + &  ! p210
         lambda3 * nodes3(5, :))  ! p201
    expected3(2, :) = ( &
         lambda1 * nodes3(2, :) + &  ! p210
         lambda2 * nodes3(3, :) + &  ! p120
         lambda3 * nodes3(6, :))  ! p111
    expected3(3, :) = ( &
         lambda1 * nodes3(3, :) + &  ! p120
         lambda2 * nodes3(4, :) + &  ! p030
         lambda3 * nodes3(7, :))  ! p021
    expected3(4, :) = ( &
         lambda1 * nodes3(5, :) + &  ! p201
         lambda2 * nodes3(6, :) + &  ! p111
         lambda3 * nodes3(8, :))  ! p102
    expected3(5, :) = ( &
         lambda1 * nodes3(6, :) + &  ! p111
         lambda2 * nodes3(7, :) + &  ! p021
         lambda3 * nodes3(9, :))  ! p012
    expected3(6, :) = ( &
         lambda1 * nodes3(8, :) + &  ! p102
         lambda2 * nodes3(9, :) + &  ! p012
         lambda3 * nodes3(10, :))  ! p003

    call de_casteljau_one_round( &
         10, 2, nodes3, 3, &
         lambda1, lambda2, lambda3, new_nodes3)
    if (all(new_nodes3 == expected3)) then
       call print_status(name, case_id, .TRUE.)
    else
       call print_status(name, case_id, .FALSE.)
       success = .FALSE.
    end if

  end subroutine test_de_casteljau_one_round

end module test_surface
