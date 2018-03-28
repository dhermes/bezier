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

  use, intrinsic :: iso_c_binding, only: c_bool, c_double, c_int, c_ptr, c_loc
  use surface, only: &
       de_casteljau_one_round, evaluate_barycentric, &
       evaluate_barycentric_multi, evaluate_cartesian_multi, jacobian_both, &
       jacobian_det, specialize_surface, subdivide_nodes, compute_edge_nodes, &
       shoelace_for_area, compute_area
  use types, only: dp
  use unit_test_helpers, only: &
       print_status, get_random_nodes, get_id_mat, ref_triangle_uniform_params
  implicit none
  private &
       test_de_casteljau_one_round, test_evaluate_barycentric, &
       test_evaluate_barycentric_multi, test_evaluate_cartesian_multi, &
       test_jacobian_both, test_jacobian_det, &
       test_specialize_surface, test_subdivide_nodes, &
       subdivide_points_check, test_compute_edge_nodes, &
       test_shoelace_for_area, test_compute_area
  public surface_all_tests

contains

  subroutine surface_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_de_casteljau_one_round(success)
    call test_evaluate_barycentric(success)
    call test_evaluate_barycentric_multi(success)
    call test_evaluate_cartesian_multi(success)
    call test_jacobian_both(success)
    call test_jacobian_det(success)
    call test_specialize_surface(success)
    call test_subdivide_nodes(success)
    call test_compute_edge_nodes(success)
    call test_shoelace_for_area(success)
    call test_compute_area(success)

  end subroutine surface_all_tests

  subroutine test_de_casteljau_one_round(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 3), nodes2(2, 6), nodes3(2, 10)
    real(c_double) :: new_nodes1(2, 1), new_nodes2(2, 3), new_nodes3(2, 6)
    real(c_double) :: expected1(2, 1), expected2(2, 3), expected3(2, 6)
    real(c_double) :: lambda1, lambda2, lambda3
    integer :: case_id
    character(22) :: name

    case_id = 1
    name = "de_casteljau_one_round"

    ! CASE 1: Linear reference triangle (just goes to point).
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    nodes1(:, 3) = [0.0_dp, 1.0_dp]
    lambda1 = 0.125_dp
    lambda2 = 0.5_dp
    lambda3 = 0.375_dp
    expected1(:, 1) = [lambda2, lambda3]
    call de_casteljau_one_round( &
         3, 2, nodes1, 1, &
         lambda1, lambda2, lambda3, new_nodes1)
    case_success = all(new_nodes1 == expected1)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic surface.
    call get_random_nodes(nodes2, 790931, 1483, num_bits=8)
    ! NOTE: Use a fixed seed so the test is deterministic and round
    !       the nodes to 8 bits of precision to avoid round-off.
    lambda1 = 0.625_dp
    lambda2 = 0.25_dp
    lambda3 = 0.125_dp
    expected2(:, 1) = ( &
         lambda1 * nodes2(:, 1) + &  ! p200
         lambda2 * nodes2(:, 2) + &  ! p110
         lambda3 * nodes2(:, 4))  ! p101
    expected2(:, 2) = ( &
         lambda1 * nodes2(:, 2) + &  ! p110
         lambda2 * nodes2(:, 3) + &  ! p020
         lambda3 * nodes2(:, 5))  ! p011
    expected2(:, 3) = ( &
         lambda1 * nodes2(:, 4) + &  ! p101
         lambda2 * nodes2(:, 5) + &  ! p011
         lambda3 * nodes2(:, 6))  ! p002
    call de_casteljau_one_round( &
         6, 2, nodes2, 2, &
         lambda1, lambda2, lambda3, new_nodes2)
    case_success = all(new_nodes2 == expected2)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Cubic surface.
    nodes3(:, 1) = [0.0_dp, 0.0_dp]
    nodes3(:, 2) = [3.25_dp, 1.5_dp]
    nodes3(:, 3) = [6.5_dp, 1.5_dp]
    nodes3(:, 4) = [10.0_dp, 0.0_dp]
    nodes3(:, 5) = [1.5_dp, 3.25_dp]
    nodes3(:, 6) = [5.0_dp, 5.0_dp]
    nodes3(:, 7) = [10.0_dp, 5.25_dp]
    nodes3(:, 8) = [1.5_dp, 6.5_dp]
    nodes3(:, 9) = [5.25_dp, 10.0_dp]
    nodes3(:, 10) = [0.0_dp, 10.0_dp]
    lambda1 = 0.375_dp
    lambda2 = 0.25_dp
    lambda3 = 0.375_dp
    expected3(:, 1) = ( &
         lambda1 * nodes3(:, 1) + &  ! p300
         lambda2 * nodes3(:, 2) + &  ! p210
         lambda3 * nodes3(:, 5))  ! p201
    expected3(:, 2) = ( &
         lambda1 * nodes3(:, 2) + &  ! p210
         lambda2 * nodes3(:, 3) + &  ! p120
         lambda3 * nodes3(:, 6))  ! p111
    expected3(:, 3) = ( &
         lambda1 * nodes3(:, 3) + &  ! p120
         lambda2 * nodes3(:, 4) + &  ! p030
         lambda3 * nodes3(:, 7))  ! p021
    expected3(:, 4) = ( &
         lambda1 * nodes3(:, 5) + &  ! p201
         lambda2 * nodes3(:, 6) + &  ! p111
         lambda3 * nodes3(:, 8))  ! p102
    expected3(:, 5) = ( &
         lambda1 * nodes3(:, 6) + &  ! p111
         lambda2 * nodes3(:, 7) + &  ! p021
         lambda3 * nodes3(:, 9))  ! p012
    expected3(:, 6) = ( &
         lambda1 * nodes3(:, 8) + &  ! p102
         lambda2 * nodes3(:, 9) + &  ! p012
         lambda3 * nodes3(:, 10))  ! p003

    call de_casteljau_one_round( &
         10, 2, nodes3, 3, &
         lambda1, lambda2, lambda3, new_nodes3)
    case_success = all(new_nodes3 == expected3)
    call print_status(name, case_id, case_success, success)

  end subroutine test_de_casteljau_one_round

  subroutine test_evaluate_barycentric(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 3), nodes2(2, 6), nodes3(3, 6)
    real(c_double) :: nodes4(2, 10), nodes5(2, 15)
    real(c_double) :: point1(2, 1), point2(3, 1)
    real(c_double) :: expected1(2, 1), expected2(3, 1)
    real(c_double) :: lambda1, lambda2, lambda3
    integer :: index_, i, j, k, trinomial
    integer :: case_id
    character(20) :: name

    case_id = 1
    name = "evaluate_barycentric"

    ! CASE 1: Linear surface.
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.5_dp]
    nodes1(:, 3) = [0.0_dp, 1.25_dp]
    lambda1 = 0.25_dp
    lambda2 = 0.5_dp
    lambda3 = 0.25_dp
    expected1(:, 1) = [0.5_dp, 0.5625_dp]
    call evaluate_barycentric( &
         3, 2, nodes1, 1, lambda1, lambda2, lambda3, point1)
    case_success = all(point1 == expected1)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic surface.
    lambda1 = 0.0_dp
    lambda2 = 0.25_dp
    lambda3 = 0.75_dp
    nodes2(:, 1) = [0.0_dp, 0.0_dp]
    nodes2(:, 2) = [0.5_dp, 0.0_dp]
    nodes2(:, 3) = [1.0_dp, 0.5_dp]
    nodes2(:, 4) = [0.5_dp, 1.25_dp]
    nodes2(:, 5) = [0.0_dp, 1.25_dp]
    nodes2(:, 6) = [0.0_dp, 0.5_dp]
    expected1(:, 1) = [0.0625_dp, 0.78125_dp]
    call evaluate_barycentric( &
         6, 2, nodes2, 2, lambda1, lambda2, lambda3, point1)
    case_success = all(point1 == expected1)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Quadratic surface in 3D.
    nodes3(:, 1) = [0.0_dp, 0.0_dp, 1.0_dp]
    nodes3(:, 2) = [0.5_dp, 0.0_dp, 0.25_dp]
    nodes3(:, 3) = [1.0_dp, 0.5_dp, 0.0_dp]
    nodes3(:, 4) = [0.5_dp, 1.25_dp, 1.25_dp]
    nodes3(:, 5) = [0.0_dp, 1.25_dp, 0.5_dp]
    nodes3(:, 6) = [0.0_dp, 0.5_dp, -1.0_dp]
    lambda1 = 0.125_dp
    lambda2 = 0.375_dp
    lambda3 = 0.5_dp
    expected2(:, 1) = [0.25_dp, 0.8203125_dp, 0.1328125_dp]
    call evaluate_barycentric( &
         6, 3, nodes3, 2, lambda1, lambda2, lambda3, point2)
    case_success = all(point2 == expected2)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Cubic surface.
    nodes4(:, 1) = [0.0_dp, 0.0_dp]
    nodes4(:, 2) = [0.25_dp, 0.0_dp]
    nodes4(:, 3) = [0.75_dp, 0.25_dp]
    nodes4(:, 4) = [1.0_dp, 0.0_dp]
    nodes4(:, 5) = [0.0_dp, 0.25_dp]
    nodes4(:, 6) = [0.375_dp, 0.25_dp]
    nodes4(:, 7) = [0.5_dp, 0.25_dp]
    nodes4(:, 8) = [0.0_dp, 0.5_dp]
    nodes4(:, 9) = [0.25_dp, 0.75_dp]
    nodes4(:, 10) = [0.0_dp, 1.0_dp]
    lambda1 = 0.125_dp
    lambda2 = 0.5_dp
    lambda3 = 0.375_dp
    expected1(:, 1) = [0.447265625_dp, 0.37060546875_dp]
    call evaluate_barycentric( &
         10, 2, nodes4, 3, lambda1, lambda2, lambda3, point1)
    case_success = all(point1 == expected1)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Quartic (random) surface.
    call get_random_nodes(nodes5, 64441, 222, num_bits=8)
    ! NOTE: Use a fixed seed so the test is deterministic and round
    !       the nodes to 8 bits of precision to avoid round-off.
    lambda1 = 0.125_dp
    lambda2 = 0.375_dp
    lambda3 = 0.5_dp
    ! We will manually compute the expected output.
    expected1 = 0
    index_ = 1
    do k = 0, 4
       do j = 0, 4 - k
          i = 4 - j - k
          if (maxval([i, j, k]) == 4) then
             trinomial = 1
          else if (maxval([i, j, k]) == 3) then
             trinomial = 4
          else if (minval([i, j, k]) == 0) then
             trinomial = 6
          else
             trinomial = 12
          end if
          expected1(:, 1) = ( &
               expected1(:, 1) + &
               ( &
               trinomial * lambda1**i * lambda2**j * lambda3**k * &
               nodes5(:, index_)))
          ! Update the index.
          index_ = index_ + 1
       end do
    end do
    call evaluate_barycentric( &
         15, 2, nodes5, 4, lambda1, lambda2, lambda3, point1)
    case_success = all(point1 == expected1)
    call print_status(name, case_id, case_success, success)

  end subroutine test_evaluate_barycentric

  subroutine test_evaluate_barycentric_multi(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes(2, 3), param_vals(3, 3)
    real(c_double) :: evaluated(2, 3), expected(2, 3)
    real(c_double) :: nodes_deg0(2, 1)
    integer :: i
    integer :: case_id
    character(26) :: name

    case_id = 1
    name = "evaluate_barycentric_multi"

    ! CASE 1: Basic evaluation.
    nodes(:, 1) = [0.0_dp, 0.0_dp]
    nodes(:, 2) = [2.0_dp, 1.0_dp]
    nodes(:, 3) = [-3.0_dp, 2.0_dp]
    param_vals(1, :) = [1.0_dp, 0.0_dp, 0.0_dp]
    param_vals(2, :) = [0.0_dp, 1.0_dp, 0.0_dp]
    param_vals(3, :) = [0.0_dp, 0.5_dp, 0.5_dp]
    expected(:, 1) = [0.0_dp, 0.0_dp]
    expected(:, 2) = [2.0_dp, 1.0_dp]
    expected(:, 3) = [-0.5_dp, 1.5_dp]
    call evaluate_barycentric_multi( &
         3, 2, nodes, 1, 3, param_vals, evaluated)
    case_success = all(evaluated == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Values outside the domain.
    nodes(:, 1) = [0.0_dp, 0.0_dp]
    nodes(:, 2) = [3.0_dp, -1.0_dp]
    nodes(:, 3) = [1.0_dp, 0.0_dp]
    param_vals(1, :) = [0.25_dp, 0.25_dp, 0.25_dp]
    param_vals(2, :) = [-1.0_dp, -1.0_dp, 3.0_dp]
    param_vals(3, :) = [0.125_dp, 0.75_dp, 0.125_dp]
    expected(:, 1) = [1.0_dp, -0.25_dp]
    expected(:, 2) = [0.0_dp, 1.0_dp]
    expected(:, 3) = [2.375_dp, -0.75_dp]
    call evaluate_barycentric_multi( &
         3, 2, nodes, 1, 3, param_vals, evaluated)
    case_success = all(evaluated == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Degree zero "surface"
    nodes_deg0(:, 1) = [1.5_dp, 7.75_dp]
    param_vals(1, :) = [0.5_dp, 0.25_dp, 0.25_dp]
    param_vals(2, :) = [0.25_dp, 0.125_dp, 0.625_dp]
    param_vals(3, :) = [0.75_dp, 0.25_dp, 0.0_dp]
    forall (i = 1:3)
       expected(:, i) = nodes_deg0(:, 1)
    end forall
    call evaluate_barycentric_multi( &
         1, 2, nodes_deg0, 0, 3, param_vals, evaluated)
    case_success = all(evaluated == expected)
    call print_status(name, case_id, case_success, success)

  end subroutine test_evaluate_barycentric_multi

  subroutine test_evaluate_cartesian_multi(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 6), param_vals1(4, 2)
    real(c_double) :: evaluated1(2, 4), expected1(2, 4)
    real(c_double) :: nodes2(2, 3), param_vals2(3, 2)
    real(c_double) :: evaluated2(2, 3), expected2(2, 3)
    real(c_double) :: nodes3(3, 1), param_vals3(16, 2)
    real(c_double) :: evaluated3(3, 16), expected3(3, 16)
    integer :: i
    integer :: case_id
    character(24) :: name

    case_id = 1
    name = "evaluate_cartesian_multi"

    ! CASE 1: Basic evaluation
    nodes1(:, 1) = [0.0_dp, 0.0_dp]
    nodes1(:, 2) = [1.0_dp, 0.75_dp]
    nodes1(:, 3) = [2.0_dp, 1.0_dp]
    nodes1(:, 4) = [-1.5_dp, 1.0_dp]
    nodes1(:, 5) = [-0.5_dp, 1.5_dp]
    nodes1(:, 6) = [-3.0_dp, 2.0_dp]
    param_vals1(1, :) = [0.25_dp, 0.75_dp]
    param_vals1(2, :) = [0.0_dp, 0.0_dp]
    param_vals1(3, :) = [0.5_dp, 0.25_dp]
    param_vals1(4, :) = [0.25_dp, 0.375_dp]
    expected1(:, 1) = [-1.75_dp, 1.75_dp]
    expected1(:, 2) = [0.0_dp, 0.0_dp]
    expected1(:, 3) = [0.25_dp, 1.0625_dp]
    expected1(:, 4) = [-0.625_dp, 1.046875_dp]
    call evaluate_cartesian_multi( &
         6, 2, nodes1, 2, 4, param_vals1, evaluated1)
    case_success = all(evaluated1 == expected1)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Values outside the domain.
    nodes2(:, 1) = [0.0_dp, 2.0_dp]
    nodes2(:, 2) = [1.0_dp, 1.0_dp]
    nodes2(:, 3) = [1.0_dp, 0.0_dp]
    param_vals2(1, :) = [0.25_dp, 0.25_dp]
    param_vals2(2, :) = [-1.0_dp, 3.0_dp]
    param_vals2(3, :) = [0.75_dp, 0.125_dp]
    expected2(:, 1) = [0.5_dp, 1.25_dp]
    expected2(:, 2) = [2.0_dp, -3.0_dp]
    expected2(:, 3) = [0.875_dp, 1.0_dp]
    call evaluate_cartesian_multi( &
         3, 2, nodes2, 1, 3, param_vals2, evaluated2)
    case_success = all(evaluated2 == expected2)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Degree zero "surface"
    nodes3(:, 1) = [-1.5_dp, 0.75_dp, 5.0_dp]
    forall (i = 1:16)
       expected3(:, i) = nodes3(:, 1)
    end forall
    call evaluate_cartesian_multi( &
         1, 3, nodes3, 0, 16, param_vals3, evaluated3)
    case_success = all(evaluated3 == expected3)
    call print_status(name, case_id, case_success, success)

  end subroutine test_evaluate_cartesian_multi

  subroutine test_jacobian_both(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(1, 3), expected1(2, 1), new_nodes1(2, 1)
    real(c_double) :: nodes2(3, 6), expected2(6, 3), new_nodes2(6, 3)
    real(c_double) :: nodes3(2, 10), expected3(4, 6), new_nodes3(4, 6)
    integer :: case_id
    character(13) :: name

    case_id = 1
    name = "jacobian_both"

    ! CASE 1: Linear surface.
    ! B(s, t) = -2s + 2t + 3
    nodes1(1, :) = [3.0_dp, 1.0_dp, 5.0_dp]
    ! B_s = -2
    ! B_t = 2
    expected1(:, 1) = [-2.0_dp, 2.0_dp]
    call jacobian_both( &
         3, 1, nodes1, 1, new_nodes1)
    case_success = all(new_nodes1 == expected1)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic surface.
    ! B(s, t) = [
    !     4 s t - 2 s + 5 t^2 - 6 t + 3,
    !     -s (s - 2 t),
    !     8 s^2 - 10 s t - 2 s - 13 t^2 + 12 t + 1,
    ! ]
    !
    nodes2(:, 1) = [3.0_dp, 0.0_dp, 1.0_dp]
    nodes2(:, 2) = [2.0_dp, 0.0_dp, 0.0_dp]
    nodes2(:, 3) = [1.0_dp, -1.0_dp, 7.0_dp]
    nodes2(:, 4) = [0.0_dp, 0.0_dp, 7.0_dp]
    nodes2(:, 5) = [1.0_dp, 1.0_dp, 1.0_dp]
    nodes2(:, 6) = [2.0_dp, 0.0_dp, 0.0_dp]
    ! B_s = [
    !     4 t - 2,
    !     -2 s + 2 t,
    !     16 s - 10 t - 2,
    ! ]
    ! B_t = [
    !     4 s + 10 t - 6,
    !     2 s,
    !    -10 s - 26 t + 12,
    ! ]
    expected2(:, 1) = [-2.0_dp, 0.0_dp, -2.0_dp, -6.0_dp, 0.0_dp, 12.0_dp]
    expected2(:, 2) = [-2.0_dp, -2.0_dp, 14.0_dp, -2.0_dp, 2.0_dp, 2.0_dp]
    expected2(:, 3) = [2.0_dp, 2.0_dp, -12.0_dp, 4.0_dp, 0.0_dp, -14.0_dp]
    call jacobian_both( &
         6, 3, nodes2, 2, new_nodes2)
    case_success = all(new_nodes2 == expected2)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Cubic surface.
    ! B(s, t) = [
    !     -2s^3 + 9s^2t + 12st^2 - 12st + 3s - 2t^3 + 6t,
    !     (-10s^3 - 30s^2t + 30s^2 - 36st^2 + 42st -
    !          18s - 15t^3 + 30t^2 - 18t + 7),
    ! ]
    nodes3(:, 1) = [0.0_dp, 7.0_dp]
    nodes3(:, 2) = [1.0_dp, 1.0_dp]
    nodes3(:, 3) = [2.0_dp, 5.0_dp]
    nodes3(:, 4) = [1.0_dp, 9.0_dp]
    nodes3(:, 5) = [2.0_dp, 1.0_dp]
    nodes3(:, 6) = [1.0_dp, 2.0_dp]
    nodes3(:, 7) = [3.0_dp, 3.0_dp]
    nodes3(:, 8) = [4.0_dp, 5.0_dp]
    nodes3(:, 9) = [5.0_dp, 1.0_dp]
    nodes3(:, 10) = [4.0_dp, 4.0_dp]
    ! B_s = [
    !     -6s^2 + 18st + 12t^2 - 12t + 3,
    !     -30s^2 - 60st + 60s - 36t^2 + 42t - 18,
    ! ]
    ! B_t = [
    !     9s^2 + 24st - 12s - 6t^2 + 6,
    !     -30s^2 - 72st + 42s - 45t^2 + 60t - 18,
    ! ]
    expected3(:, 1) = [3.0_dp, -18.0_dp, 6.0_dp, -18.0_dp]
    expected3(:, 2) = [3.0_dp, 12.0_dp, 0.0_dp, 3.0_dp]
    expected3(:, 3) = [-3.0_dp, 12.0_dp, 3.0_dp, -6.0_dp]
    expected3(:, 4) = [-3.0_dp, 3.0_dp, 6.0_dp, 12.0_dp]
    expected3(:, 5) = [6.0_dp, 3.0_dp, 12.0_dp, -3.0_dp]
    expected3(:, 6) = [3.0_dp, -12.0_dp, 0.0_dp, -3.0_dp]
    call jacobian_both( &
         10, 2, nodes3, 3, new_nodes3)
    case_success = all(new_nodes3 == expected3)
    call print_status(name, case_id, case_success, success)

  end subroutine test_jacobian_both

  subroutine test_jacobian_det(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 3), param_vals1(13, 2), evaluated1(13)
    real(c_double) :: nodes2(2, 6), param_vals2(4, 2)
    real(c_double) :: expected2(4), evaluated2(4)
    integer :: case_id
    character(12) :: name

    case_id = 1
    name = "jacobian_det"

    ! CASE 1: Linear surface (determinant is area).
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    nodes1(:, 3) = [0.0_dp, 2.0_dp]
    call get_random_nodes(param_vals1, 308975611, 101301, num_bits=8)
    call jacobian_det( &
         3, nodes1, 1, 13, param_vals1, evaluated1)
    case_success = all(evaluated1 == 2.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic (i.e. non-linear) surface.
    nodes2(:, 1) = 0
    nodes2(:, 2) = [0.5_dp, 0.0_dp]
    nodes2(:, 3) = [1.0_dp, 0.0_dp]
    nodes2(:, 4) = [0.0_dp, 0.5_dp]
    nodes2(:, 5) = [1.0_dp, 1.0_dp]
    nodes2(:, 6) = [0.0_dp, 1.0_dp]
    ! B(s, t) = [s(t + 1), t(s + 1)]
    param_vals2(1, :) = [0.125_dp, 0.125_dp]
    param_vals2(2, :) = [0.5_dp, 0.375_dp]
    param_vals2(3, :) = [0.25_dp, 0.75_dp]
    param_vals2(4, :) = [1.0_dp, 0.0_dp]
    ! det(DB) = s + t + 1
    expected2 = param_vals2(:, 1) + param_vals2(:, 2) + 1.0_dp
    call jacobian_det( &
         6, nodes2, 2, 4, param_vals2, evaluated2)
    case_success = all(evaluated2 == expected2)
    call print_status(name, case_id, case_success, success)

  end subroutine test_jacobian_det

  subroutine test_specialize_surface(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    real(c_double), allocatable :: id_mat(:, :), specialized(:, :)
    real(c_double), allocatable :: expected(:, :)
    real(c_double) :: &
         weights0(3), weights1(3), weights2(3), weights3(3), &
         weights4(3), weights5(3)
    character(18) :: name

    case_id = 1
    name = "specialize_surface"

    weights0 = [1.0_dp, 0.0_dp, 0.0_dp]
    weights1 = [0.5_dp, 0.5_dp, 0.0_dp]
    weights2 = [0.0_dp, 1.0_dp, 0.0_dp]
    weights3 = [0.5_dp, 0.0_dp, 0.5_dp]
    weights4 = [0.0_dp, 0.5_dp, 0.5_dp]
    weights5 = [0.0_dp, 0.0_dp, 1.0_dp]

    ! CASE 1: degree == 1, subdivision A
    allocate(id_mat(3, 3))
    id_mat = get_id_mat(3)
    allocate(specialized(3, 3))
    allocate(expected(3, 3))
    call specialize_surface( &
         3, 3, id_mat, 1, &
         weights0, weights1, weights3, &
         specialized)
    expected(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp]
    expected(:, 2) = [0.5_dp, 0.5_dp, 0.0_dp]
    expected(:, 3) = [0.5_dp, 0.0_dp, 0.5_dp]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: degree == 1, subdivision B
    call specialize_surface( &
         3, 3, id_mat, 1, &
         weights4, weights3, weights1, &
         specialized)
    expected(:, 1) = [0.0_dp, 0.5_dp, 0.5_dp]
    expected(:, 2) = [0.5_dp, 0.0_dp, 0.5_dp]
    expected(:, 3) = [0.5_dp, 0.5_dp, 0.0_dp]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: degree == 1, subdivision C
    call specialize_surface( &
         3, 3, id_mat, 1, &
         weights1, weights2, weights4, &
         specialized)
    expected(:, 1) = [0.5_dp, 0.5_dp, 0.0_dp]
    expected(:, 2) = [0.0_dp, 1.0_dp, 0.0_dp]
    expected(:, 3) = [0.0_dp, 0.5_dp, 0.5_dp]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: degree == 1, subdivision D
    call specialize_surface( &
         3, 3, id_mat, 1, &
         weights3, weights4, weights5, &
         specialized)
    expected(:, 1) = [0.5_dp, 0.0_dp, 0.5_dp]
    expected(:, 2) = [0.0_dp, 0.5_dp, 0.5_dp]
    expected(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! De-allocate and prepare for degree == 2
    deallocate(id_mat)
    deallocate(specialized)
    deallocate(expected)

    ! CASE 5: degree == 2, subdivision A
    id_mat = get_id_mat(6)
    allocate(specialized(6, 6))
    allocate(expected(6, 6))
    call specialize_surface( &
         6, 6, id_mat, 2, &
         weights0, weights1, weights3, &
         specialized)
    expected(:, 1) = [1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    expected(:, 2) = [0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    expected(:, 3) = [0.25_dp, 0.5_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    expected(:, 4) = [0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp]
    expected(:, 5) = [0.25_dp, 0.25_dp, 0.0_dp, 0.25_dp, 0.25_dp, 0.0_dp]
    expected(:, 6) = [0.25_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.25_dp]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: degree == 2, subdivision B
    call specialize_surface( &
         6, 6, id_mat, 2, &
         weights4, weights3, weights1, &
         specialized)
    expected(:, 1) = [0.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, 0.5_dp, 0.25_dp]
    expected(:, 2) = [0.0_dp, 0.25_dp, 0.0_dp, 0.25_dp, 0.25_dp, 0.25_dp]
    expected(:, 3) = [0.25_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.25_dp]
    expected(:, 4) = [0.0_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.0_dp]
    expected(:, 5) = [0.25_dp, 0.25_dp, 0.0_dp, 0.25_dp, 0.25_dp, 0.0_dp]
    expected(:, 6) = [0.25_dp, 0.5_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: degree == 2, subdivision C
    call specialize_surface( &
         6, 6, id_mat, 2, &
         weights1, weights2, weights4, &
         specialized)
    expected(:, 1) = [0.25_dp, 0.5_dp, 0.25_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    expected(:, 2) = [0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    expected(:, 3) = [0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    expected(:, 4) = [0.0_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp, 0.0_dp]
    expected(:, 5) = [0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 0.0_dp]
    expected(:, 6) = [0.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, 0.5_dp, 0.25_dp]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: degree == 2, subdivision D
    call specialize_surface( &
         6, 6, id_mat, 2, &
         weights3, weights4, weights5, &
         specialized)
    expected(:, 1) = [0.25_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.25_dp]
    expected(:, 2) = [0.0_dp, 0.25_dp, 0.0_dp, 0.25_dp, 0.25_dp, 0.25_dp]
    expected(:, 3) = [0.0_dp, 0.0_dp, 0.25_dp, 0.0_dp, 0.5_dp, 0.25_dp]
    expected(:, 4) = [0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp]
    expected(:, 5) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.5_dp]
    expected(:, 6) = [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! De-allocate and prepare for degree == 3
    deallocate(id_mat)
    deallocate(specialized)
    deallocate(expected)

    ! CASE 9: degree == 3, subdivision A
    id_mat = get_id_mat(10)
    allocate(specialized(10, 10))
    allocate(expected(10, 10))
    call specialize_surface( &
         10, 10, id_mat, 3, &
         weights0, weights1, weights3, &
         specialized)
    expected(:, 1) = 0.125_dp * [8, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 2) = 0.125_dp * [4, 4, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 3) = 0.125_dp * [2, 4, 2, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 4) = 0.125_dp * [1, 3, 3, 1, 0, 0, 0, 0, 0, 0]
    expected(:, 5) = 0.125_dp * [4, 0, 0, 0, 4, 0, 0, 0, 0, 0]
    expected(:, 6) = 0.125_dp * [2, 2, 0, 0, 2, 2, 0, 0, 0, 0]
    expected(:, 7) = 0.125_dp * [1, 2, 1, 0, 1, 2, 1, 0, 0, 0]
    expected(:, 8) = 0.125_dp * [2, 0, 0, 0, 4, 0, 0, 2, 0, 0]
    expected(:, 9) = 0.125_dp * [1, 1, 0, 0, 2, 2, 0, 1, 1, 0]
    expected(:, 10) = 0.125_dp * [1, 0, 0, 0, 3, 0, 0, 3, 0, 1]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 10: degree == 3, subdivision B
    call specialize_surface( &
         10, 10, id_mat, 3, &
         weights4, weights3, weights1, &
         specialized)
    expected(:, 1) = 0.125_dp * [0, 0, 0, 1, 0, 0, 3, 0, 3, 1]
    expected(:, 2) = 0.125_dp * [0, 0, 1, 0, 0, 2, 1, 1, 2, 1]
    expected(:, 3) = 0.125_dp * [0, 1, 0, 0, 1, 2, 0, 2, 1, 1]
    expected(:, 4) = 0.125_dp * [1, 0, 0, 0, 3, 0, 0, 3, 0, 1]
    expected(:, 5) = 0.125_dp * [0, 0, 1, 1, 0, 2, 2, 1, 1, 0]
    expected(:, 6) = 0.125_dp * [0, 1, 1, 0, 1, 2, 1, 1, 1, 0]
    expected(:, 7) = 0.125_dp * [1, 1, 0, 0, 2, 2, 0, 1, 1, 0]
    expected(:, 8) = 0.125_dp * [0, 1, 2, 1, 1, 2, 1, 0, 0, 0]
    expected(:, 9) = 0.125_dp * [1, 2, 1, 0, 1, 2, 1, 0, 0, 0]
    expected(:, 10) = 0.125_dp * [1, 3, 3, 1, 0, 0, 0, 0, 0, 0]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 11: degree == 3, subdivision C
    call specialize_surface( &
         10, 10, id_mat, 3, &
         weights1, weights2, weights4, &
         specialized)
    expected(:, 1) = 0.125_dp * [1, 3, 3, 1, 0, 0, 0, 0, 0, 0]
    expected(:, 2) = 0.125_dp * [0, 2, 4, 2, 0, 0, 0, 0, 0, 0]
    expected(:, 3) = 0.125_dp * [0, 0, 4, 4, 0, 0, 0, 0, 0, 0]
    expected(:, 4) = 0.125_dp * [0, 0, 0, 8, 0, 0, 0, 0, 0, 0]
    expected(:, 5) = 0.125_dp * [0, 1, 2, 1, 1, 2, 1, 0, 0, 0]
    expected(:, 6) = 0.125_dp * [0, 0, 2, 2, 0, 2, 2, 0, 0, 0]
    expected(:, 7) = 0.125_dp * [0, 0, 0, 4, 0, 0, 4, 0, 0, 0]
    expected(:, 8) = 0.125_dp * [0, 0, 1, 1, 0, 2, 2, 1, 1, 0]
    expected(:, 9) = 0.125_dp * [0, 0, 0, 2, 0, 0, 4, 0, 2, 0]
    expected(:, 10) = 0.125_dp * [0, 0, 0, 1, 0, 0, 3, 0, 3, 1]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 12: degree == 3, subdivision D
    call specialize_surface( &
         10, 10, id_mat, 3, &
         weights3, weights4, weights5, &
         specialized)
    expected(:, 1) = 0.125_dp * [1, 0, 0, 0, 3, 0, 0, 3, 0, 1]
    expected(:, 2) = 0.125_dp * [0, 1, 0, 0, 1, 2, 0, 2, 1, 1]
    expected(:, 3) = 0.125_dp * [0, 0, 1, 0, 0, 2, 1, 1, 2, 1]
    expected(:, 4) = 0.125_dp * [0, 0, 0, 1, 0, 0, 3, 0, 3, 1]
    expected(:, 5) = 0.125_dp * [0, 0, 0, 0, 2, 0, 0, 4, 0, 2]
    expected(:, 6) = 0.125_dp * [0, 0, 0, 0, 0, 2, 0, 2, 2, 2]
    expected(:, 7) = 0.125_dp * [0, 0, 0, 0, 0, 0, 2, 0, 4, 2]
    expected(:, 8) = 0.125_dp * [0, 0, 0, 0, 0, 0, 0, 4, 0, 4]
    expected(:, 9) = 0.125_dp * [0, 0, 0, 0, 0, 0, 0, 0, 4, 4]
    expected(:, 10) = 0.125_dp * [0, 0, 0, 0, 0, 0, 0, 0, 0, 8]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! De-allocate and prepare for degree == 4
    deallocate(id_mat)
    deallocate(specialized)
    deallocate(expected)

    ! CASE 13: degree == 4, subdivision A
    id_mat = get_id_mat(15)
    allocate(specialized(15, 15))
    allocate(expected(15, 15))
    call specialize_surface( &
         15, 15, id_mat, 4, &
         weights0, weights1, weights3, &
         specialized)
    expected(:, 1) = 0.0625_dp * [16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 2) = 0.0625_dp * [8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 3) = 0.0625_dp * [4, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 4) = 0.0625_dp * [2, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 5) = 0.0625_dp * [1, 4, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 6) = 0.0625_dp * [8, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 7) = 0.0625_dp * [4, 4, 0, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 8) = 0.0625_dp * [2, 4, 2, 0, 0, 2, 4, 2, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 9) = 0.0625_dp * [1, 3, 3, 1, 0, 1, 3, 3, 1, 0, 0, 0, 0, 0, 0]
    expected(:, 10) = 0.0625_dp * [4, 0, 0, 0, 0, 8, 0, 0, 0, 4, 0, 0, 0, 0, 0]
    expected(:, 11) = 0.0625_dp * [2, 2, 0, 0, 0, 4, 4, 0, 0, 2, 2, 0, 0, 0, 0]
    expected(:, 12) = 0.0625_dp * [1, 2, 1, 0, 0, 2, 4, 2, 0, 1, 2, 1, 0, 0, 0]
    expected(:, 13) = 0.0625_dp * [2, 0, 0, 0, 0, 6, 0, 0, 0, 6, 0, 0, 2, 0, 0]
    expected(:, 14) = 0.0625_dp * [1, 1, 0, 0, 0, 3, 3, 0, 0, 3, 3, 0, 1, 1, 0]
    expected(:, 15) = 0.0625_dp * [1, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0, 4, 0, 1]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 14: degree == 4, subdivision B
    call specialize_surface( &
         15, 15, id_mat, 4, &
         weights4, weights3, weights1, &
         specialized)
    expected(:, 1) = 0.0625_dp * [0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 6, 0, 4, 1]
    expected(:, 2) = 0.0625_dp * [0, 0, 0, 1, 0, 0, 0, 3, 1, 0, 3, 3, 1, 3, 1]
    expected(:, 3) = 0.0625_dp * [0, 0, 1, 0, 0, 0, 2, 2, 0, 1, 4, 1, 2, 2, 1]
    expected(:, 4) = 0.0625_dp * [0, 1, 0, 0, 0, 1, 3, 0, 0, 3, 3, 0, 3, 1, 1]
    expected(:, 5) = 0.0625_dp * [1, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0, 4, 0, 1]
    expected(:, 6) = 0.0625_dp * [0, 0, 0, 1, 1, 0, 0, 3, 3, 0, 3, 3, 1, 1, 0]
    expected(:, 7) = 0.0625_dp * [0, 0, 1, 1, 0, 0, 2, 3, 1, 1, 3, 2, 1, 1, 0]
    expected(:, 8) = 0.0625_dp * [0, 1, 1, 0, 0, 1, 3, 2, 0, 2, 3, 1, 1, 1, 0]
    expected(:, 9) = 0.0625_dp * [1, 1, 0, 0, 0, 3, 3, 0, 0, 3, 3, 0, 1, 1, 0]
    expected(:, 10) = 0.0625_dp * [0, 0, 1, 2, 1, 0, 2, 4, 2, 1, 2, 1, 0, 0, 0]
    expected(:, 11) = 0.0625_dp * [0, 1, 2, 1, 0, 1, 3, 3, 1, 1, 2, 1, 0, 0, 0]
    expected(:, 12) = 0.0625_dp * [1, 2, 1, 0, 0, 2, 4, 2, 0, 1, 2, 1, 0, 0, 0]
    expected(:, 13) = 0.0625_dp * [0, 1, 3, 3, 1, 1, 3, 3, 1, 0, 0, 0, 0, 0, 0]
    expected(:, 14) = 0.0625_dp * [1, 3, 3, 1, 0, 1, 3, 3, 1, 0, 0, 0, 0, 0, 0]
    expected(:, 15) = 0.0625_dp * [1, 4, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 15: degree == 4, subdivision C
    call specialize_surface( &
         15, 15, id_mat, 4, &
         weights1, weights2, weights4, &
         specialized)
    expected(:, 1) = 0.0625_dp * [1, 4, 6, 4, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 2) = 0.0625_dp * [0, 2, 6, 6, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 3) = 0.0625_dp * [0, 0, 4, 8, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 4) = 0.0625_dp * [0, 0, 0, 8, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 5) = 0.0625_dp * [0, 0, 0, 0, 16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    expected(:, 6) = 0.0625_dp * [0, 1, 3, 3, 1, 1, 3, 3, 1, 0, 0, 0, 0, 0, 0]
    expected(:, 7) = 0.0625_dp * [0, 0, 2, 4, 2, 0, 2, 4, 2, 0, 0, 0, 0, 0, 0]
    expected(:, 8) = 0.0625_dp * [0, 0, 0, 4, 4, 0, 0, 4, 4, 0, 0, 0, 0, 0, 0]
    expected(:, 9) = 0.0625_dp * [0, 0, 0, 0, 8, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0]
    expected(:, 10) = 0.0625_dp * [0, 0, 1, 2, 1, 0, 2, 4, 2, 1, 2, 1, 0, 0, 0]
    expected(:, 11) = 0.0625_dp * [0, 0, 0, 2, 2, 0, 0, 4, 4, 0, 2, 2, 0, 0, 0]
    expected(:, 12) = 0.0625_dp * [0, 0, 0, 0, 4, 0, 0, 0, 8, 0, 0, 4, 0, 0, 0]
    expected(:, 13) = 0.0625_dp * [0, 0, 0, 1, 1, 0, 0, 3, 3, 0, 3, 3, 1, 1, 0]
    expected(:, 14) = 0.0625_dp * [0, 0, 0, 0, 2, 0, 0, 0, 6, 0, 0, 6, 0, 2, 0]
    expected(:, 15) = 0.0625_dp * [0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 6, 0, 4, 1]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

    ! CASE 16: degree == 4, subdivision D
    call specialize_surface( &
         15, 15, id_mat, 4, &
         weights3, weights4, weights5, &
         specialized)
    expected(:, 1) = 0.0625_dp * [1, 0, 0, 0, 0, 4, 0, 0, 0, 6, 0, 0, 4, 0, 1]
    expected(:, 2) = 0.0625_dp * [0, 1, 0, 0, 0, 1, 3, 0, 0, 3, 3, 0, 3, 1, 1]
    expected(:, 3) = 0.0625_dp * [0, 0, 1, 0, 0, 0, 2, 2, 0, 1, 4, 1, 2, 2, 1]
    expected(:, 4) = 0.0625_dp * [0, 0, 0, 1, 0, 0, 0, 3, 1, 0, 3, 3, 1, 3, 1]
    expected(:, 5) = 0.0625_dp * [0, 0, 0, 0, 1, 0, 0, 0, 4, 0, 0, 6, 0, 4, 1]
    expected(:, 6) = 0.0625_dp * [0, 0, 0, 0, 0, 2, 0, 0, 0, 6, 0, 0, 6, 0, 2]
    expected(:, 7) = 0.0625_dp * [0, 0, 0, 0, 0, 0, 2, 0, 0, 2, 4, 0, 4, 2, 2]
    expected(:, 8) = 0.0625_dp * [0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 4, 2, 2, 4, 2]
    expected(:, 9) = 0.0625_dp * [0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 6, 0, 6, 2]
    expected(:, 10) = 0.0625_dp * [0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 8, 0, 4]
    expected(:, 11) = 0.0625_dp * [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 4, 4, 4]
    expected(:, 12) = 0.0625_dp * [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 8, 4]
    expected(:, 13) = 0.0625_dp * [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 0, 8]
    expected(:, 14) = 0.0625_dp * [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8, 8]
    expected(:, 15) = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]
    case_success = all(specialized == expected)
    call print_status(name, case_id, case_success, success)

  end subroutine test_specialize_surface

  subroutine test_subdivide_nodes(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 3)
    real(c_double) :: nodes_a1(2, 3), nodes_b1(2, 3)
    real(c_double) :: nodes_c1(2, 3), nodes_d1(2, 3)
    real(c_double) :: expected_a1(2, 3), expected_b1(2, 3)
    real(c_double) :: expected_c1(2, 3), expected_d1(2, 3)
    real(c_double) :: nodes2(2, 6)
    real(c_double) :: nodes_a2(2, 6), nodes_b2(2, 6)
    real(c_double) :: nodes_c2(2, 6), nodes_d2(2, 6)
    real(c_double) :: expected_a2(2, 6), expected_b2(2, 6)
    real(c_double) :: expected_c2(2, 6), expected_d2(2, 6)
    real(c_double) :: nodes3(2, 10)
    real(c_double) :: nodes_a3(2, 10), nodes_b3(2, 10)
    real(c_double) :: nodes_c3(2, 10), nodes_d3(2, 10)
    real(c_double) :: expected_a3(2, 10), expected_b3(2, 10)
    real(c_double) :: expected_c3(2, 10), expected_d3(2, 10)
    integer :: case_id
    character(25) :: name

    case_id = 1
    name = "subdivide_nodes (Surface)"

    ! CASE 1: Linear surface.
    nodes1(:, 1) = 0
    nodes1(:, 2) = [1.0_dp, 0.0_dp]
    nodes1(:, 3) = [0.0_dp, 1.0_dp]
    expected_a1(:, 1) = 0
    expected_a1(:, 2) = [0.5_dp, 0.0_dp]
    expected_a1(:, 3) = [0.0_dp, 0.5_dp]
    expected_b1(:, 1) = 0.5_dp
    expected_b1(:, 2) = [0.0_dp, 0.5_dp]
    expected_b1(:, 3) = [0.5_dp, 0.0_dp]
    expected_c1(:, 1) = [0.5_dp, 0.0_dp]
    expected_c1(:, 2) = [1.0_dp, 0.0_dp]
    expected_c1(:, 3) = 0.5_dp
    expected_d1(:, 1) = [0.0_dp, 0.5_dp]
    expected_d1(:, 2) = 0.5_dp
    expected_d1(:, 3) = [0.0_dp, 1.0_dp]
    call subdivide_nodes( &
         3, 2, nodes1, 1, &
         nodes_a1, nodes_b1, nodes_c1, nodes_d1)
    case_success = ( &
         all(nodes_a1 == expected_a1) .AND. &
         all(nodes_b1 == expected_b1) .AND. &
         all(nodes_c1 == expected_c1) .AND. &
         all(nodes_d1 == expected_d1))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Evaluate subdivided parts of a linear surface.
    call subdivide_points_check( &
         3, 2, 1, 290289, 111228, case_success)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Quadratic surface.
    nodes2(:, 1) = [0.0_dp, 0.0_dp]
    nodes2(:, 2) = [0.5_dp, 0.25_dp]
    nodes2(:, 3) = [1.0_dp, 0.0_dp]
    nodes2(:, 4) = [0.5_dp, 0.75_dp]
    nodes2(:, 5) = [0.0_dp, 1.0_dp]
    nodes2(:, 6) = [0.0_dp, 0.5_dp]
    expected_a2(:, 1) = [0.0_dp, 0.0_dp]
    expected_a2(:, 2) = [0.25_dp, 0.125_dp]
    expected_a2(:, 3) = [0.5_dp, 0.125_dp]
    expected_a2(:, 4) = [0.25_dp, 0.375_dp]
    expected_a2(:, 5) = [0.25_dp, 0.5_dp]
    expected_a2(:, 6) = [0.25_dp, 0.5_dp]
    expected_b2(:, 1) = [0.25_dp, 0.625_dp]
    expected_b2(:, 2) = [0.25_dp, 0.625_dp]
    expected_b2(:, 3) = [0.25_dp, 0.5_dp]
    expected_b2(:, 4) = [0.5_dp, 0.5_dp]
    expected_b2(:, 5) = [0.25_dp, 0.5_dp]
    expected_b2(:, 6) = [0.5_dp, 0.125_dp]
    expected_c2(:, 1) = [0.5_dp, 0.125_dp]
    expected_c2(:, 2) = [0.75_dp, 0.125_dp]
    expected_c2(:, 3) = [1.0_dp, 0.0_dp]
    expected_c2(:, 4) = [0.5_dp, 0.5_dp]
    expected_c2(:, 5) = [0.5_dp, 0.5_dp]
    expected_c2(:, 6) = [0.25_dp, 0.625_dp]
    expected_d2(:, 1) = [0.25_dp, 0.5_dp]
    expected_d2(:, 2) = [0.25_dp, 0.625_dp]
    expected_d2(:, 3) = [0.25_dp, 0.625_dp]
    expected_d2(:, 4) = [0.25_dp, 0.625_dp]
    expected_d2(:, 5) = [0.0_dp, 0.75_dp]
    expected_d2(:, 6) = [0.0_dp, 0.5_dp]
    call subdivide_nodes( &
         6, 2, nodes2, 2, &
         nodes_a2, nodes_b2, nodes_c2, nodes_d2)
    case_success = ( &
         all(nodes_a2 == expected_a2) .AND. &
         all(nodes_b2 == expected_b2) .AND. &
         all(nodes_c2 == expected_c2) .AND. &
         all(nodes_d2 == expected_d2))
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Evaluate subdivided parts of a quadratic surface.
    call subdivide_points_check( &
         6, 2, 2, 219803, 7086, case_success)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Cubic surface.
    nodes3(:, 1) = [0.0_dp, 0.0_dp]
    nodes3(:, 2) = [3.25_dp, 1.5_dp]
    nodes3(:, 3) = [6.5_dp, 1.5_dp]
    nodes3(:, 4) = [10.0_dp, 0.0_dp]
    nodes3(:, 5) = [1.5_dp, 3.25_dp]
    nodes3(:, 6) = [5.0_dp, 5.0_dp]
    nodes3(:, 7) = [10.0_dp, 5.25_dp]
    nodes3(:, 8) = [1.5_dp, 6.5_dp]
    nodes3(:, 9) = [5.25_dp, 10.0_dp]
    nodes3(:, 10) = [0.0_dp, 10.0_dp]
    expected_a3(:, 1) = [0.0_dp, 0.0_dp]
    expected_a3(:, 2) = [1.625_dp, 0.75_dp]
    expected_a3(:, 3) = [3.25_dp, 1.125_dp]
    expected_a3(:, 4) = [4.90625_dp, 1.125_dp]
    expected_a3(:, 5) = [0.75_dp, 1.625_dp]
    expected_a3(:, 6) = [2.4375_dp, 2.4375_dp]
    expected_a3(:, 7) = [4.3125_dp, 2.875_dp]
    expected_a3(:, 8) = [1.125_dp, 3.25_dp]
    expected_a3(:, 9) = [2.875_dp, 4.3125_dp]
    expected_a3(:, 10) = [1.125_dp, 4.90625_dp]
    expected_b3(:, 1) = [6.96875_dp, 6.96875_dp]
    expected_b3(:, 2) = [4.8125_dp, 6.65625_dp]
    expected_b3(:, 3) = [2.875_dp, 5.96875_dp]
    expected_b3(:, 4) = [1.125_dp, 4.90625_dp]
    expected_b3(:, 5) = [6.65625_dp, 4.8125_dp]
    expected_b3(:, 6) = [4.75_dp, 4.75_dp]
    expected_b3(:, 7) = [2.875_dp, 4.3125_dp]
    expected_b3(:, 8) = [5.96875_dp, 2.875_dp]
    expected_b3(:, 9) = [4.3125_dp, 2.875_dp]
    expected_b3(:, 10) = [4.90625_dp, 1.125_dp]
    expected_c3(:, 1) = [4.90625_dp, 1.125_dp]
    expected_c3(:, 2) = [6.5625_dp, 1.125_dp]
    expected_c3(:, 3) = [8.25_dp, 0.75_dp]
    expected_c3(:, 4) = [10.0_dp, 0.0_dp]
    expected_c3(:, 5) = [5.96875_dp, 2.875_dp]
    expected_c3(:, 6) = [7.875_dp, 2.9375_dp]
    expected_c3(:, 7) = [10.0_dp, 2.625_dp]
    expected_c3(:, 8) = [6.65625_dp, 4.8125_dp]
    expected_c3(:, 9) = [8.8125_dp, 5.125_dp]
    expected_c3(:, 10) = [6.96875_dp, 6.96875_dp]
    expected_d3(:, 1) = [1.125_dp, 4.90625_dp]
    expected_d3(:, 2) = [2.875_dp, 5.96875_dp]
    expected_d3(:, 3) = [4.8125_dp, 6.65625_dp]
    expected_d3(:, 4) = [6.96875_dp, 6.96875_dp]
    expected_d3(:, 5) = [1.125_dp, 6.5625_dp]
    expected_d3(:, 6) = [2.9375_dp, 7.875_dp]
    expected_d3(:, 7) = [5.125_dp, 8.8125_dp]
    expected_d3(:, 8) = [0.75_dp, 8.25_dp]
    expected_d3(:, 9) = [2.625_dp, 10.0_dp]
    expected_d3(:, 10) = [0.0_dp, 10.0_dp]
    call subdivide_nodes( &
         10, 2, nodes3, 3, &
         nodes_a3, nodes_b3, nodes_c3, nodes_d3)
    case_success = ( &
         all(nodes_a3 == expected_a3) .AND. &
         all(nodes_b3 == expected_b3) .AND. &
         all(nodes_c3 == expected_c3) .AND. &
         all(nodes_d3 == expected_d3))
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Evaluate subdivided parts of a cubic surface.
    call subdivide_points_check( &
         10, 2, 3, 439028340, 2184938, case_success)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Evaluate subdivided parts of a quartic surface.
    call subdivide_points_check( &
         15, 2, 4, 1029038, 890012, case_success)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Evaluate subdivided parts of a quantic surface. This
    !         will use ``specialize_surface`` rather than a hard-coded
    !         process.
    call subdivide_points_check( &
         21, 2, 5, 91889, 1303497, case_success)
    call print_status(name, case_id, case_success, success)

  end subroutine test_subdivide_nodes

  subroutine subdivide_points_check( &
       num_nodes, dimension_, degree, multiplier, modulus, success)
    integer, intent(in) :: num_nodes, dimension_, degree
    integer, intent(in) :: multiplier, modulus
    logical, intent(out) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes(dimension_, num_nodes)
    real(c_double) :: nodes_a(dimension_, num_nodes)
    real(c_double) :: nodes_b(dimension_, num_nodes)
    real(c_double) :: nodes_c(dimension_, num_nodes)
    real(c_double) :: nodes_d(dimension_, num_nodes)
    real(c_double) :: ref_triangle(561, 2), local_vals(561, 2)
    real(c_double) :: evaluated1(dimension_, 561)
    real(c_double) :: evaluated2(dimension_, 561)

    success = .TRUE.

    call get_random_nodes(nodes, multiplier, modulus, num_bits=8)
    call subdivide_nodes( &
         num_nodes, dimension_, nodes, degree, &
         nodes_a, nodes_b, nodes_c, nodes_d)

    ref_triangle = ref_triangle_uniform_params(5)

    ! Evaluate the **original** at ``local_vals`` and the subdivided version
    ! on the full unit triangle.
    call evaluate_cartesian_multi( &
         num_nodes, dimension_, nodes_a, degree, &
         561, ref_triangle, evaluated1)
    local_vals = 0.5_dp * ref_triangle
    call evaluate_cartesian_multi( &
         num_nodes, dimension_, nodes, degree, &
         561, local_vals, evaluated2)
    success = ( &
         success .AND. &
         all(evaluated1 == evaluated2))

    call evaluate_cartesian_multi( &
         num_nodes, dimension_, nodes_b, degree, &
         561, ref_triangle, evaluated1)
    local_vals = 0.5_dp - 0.5_dp * ref_triangle
    call evaluate_cartesian_multi( &
         num_nodes, dimension_, nodes, degree, &
         561, local_vals, evaluated2)
    success = ( &
         success .AND. &
         all(evaluated1 == evaluated2))

    call evaluate_cartesian_multi( &
         num_nodes, dimension_, nodes_c, degree, &
         561, ref_triangle, evaluated1)
    local_vals = 0.5_dp * ref_triangle
    local_vals(:, 1) = local_vals(:, 1) + 0.5_dp
    call evaluate_cartesian_multi( &
         num_nodes, dimension_, nodes, degree, &
         561, local_vals, evaluated2)
    success = ( &
         success .AND. &
         all(evaluated1 == evaluated2))

    call evaluate_cartesian_multi( &
         num_nodes, dimension_, nodes_d, degree, &
         561, ref_triangle, evaluated1)
    local_vals = 0.5_dp * ref_triangle
    local_vals(:, 2) = local_vals(:, 2) + 0.5_dp
    call evaluate_cartesian_multi( &
         num_nodes, dimension_, nodes, degree, &
         561, local_vals, evaluated2)
    success = ( &
         success .AND. &
         all(evaluated1 == evaluated2))

  end subroutine subdivide_points_check

  subroutine test_compute_edge_nodes(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(18) :: name
    real(c_double), allocatable :: nodes(:, :)
    real(c_double), allocatable :: nodes1(:, :), nodes2(:, :), nodes3(:, :)

    case_id = 1
    name = "compute_edge_nodes"

    ! CASE 1: Linear surface.
    allocate(nodes(2, 3))
    allocate(nodes1(2, 2))
    allocate(nodes2(2, 2))
    allocate(nodes3(2, 2))
    nodes(:, 1) = [1.0_dp, 2.0_dp]
    nodes(:, 2) = [4.0_dp, 2.5_dp]
    nodes(:, 3) = [0.0_dp, 4.0_dp]

    call compute_edge_nodes( &
         3, 2, nodes, 1, nodes1, nodes2, nodes3)
    case_success = ( &
         all(nodes1 == nodes(:, [1, 2])) .AND. &
         all(nodes2 == nodes(:, [2, 3])) .AND. &
         all(nodes3 == nodes(:, [3, 1])))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic surface.
    deallocate(nodes)
    deallocate(nodes1)
    deallocate(nodes2)
    deallocate(nodes3)
    allocate(nodes(2, 6))
    allocate(nodes1(2, 3))
    allocate(nodes2(2, 3))
    allocate(nodes3(2, 3))
    nodes(:, 1) = 0
    nodes(:, 2) = [1.25_dp, 0.5_dp]
    nodes(:, 3) = [2.0_dp, 1.0_dp]
    nodes(:, 4) = [-1.5_dp, 0.75_dp]
    nodes(:, 5) = [0.0_dp, 2.0_dp]
    nodes(:, 6) = [-3.0_dp, 3.0_dp]

    call compute_edge_nodes( &
         6, 2, nodes, 2, nodes1, nodes2, nodes3)
    case_success = ( &
         all(nodes1 == nodes(:, [1, 2, 3])) .AND. &
         all(nodes2 == nodes(:, [3, 5, 6])) .AND. &
         all(nodes3 == nodes(:, [6, 4, 1])))
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Cubic surface.
    deallocate(nodes)
    deallocate(nodes1)
    deallocate(nodes2)
    deallocate(nodes3)
    allocate(nodes(2, 10))
    allocate(nodes1(2, 4))
    allocate(nodes2(2, 4))
    allocate(nodes3(2, 4))
    nodes(:, 1) = [0.0_dp, 0.0_dp]
    nodes(:, 2) = [0.328125_dp, 0.1484375_dp]
    nodes(:, 3) = [0.65625_dp, 0.1484375_dp]
    nodes(:, 4) = [1.0_dp, 0.0_dp]
    nodes(:, 5) = [0.1484375_dp, 0.328125_dp]
    nodes(:, 6) = [0.5_dp, 0.5_dp]
    nodes(:, 7) = [1.0_dp, 0.53125_dp]
    nodes(:, 8) = [0.1484375_dp, 0.65625_dp]
    nodes(:, 9) = [0.53125_dp, 1.0_dp]
    nodes(:, 10) = [0.0_dp, 1.0_dp]

    call compute_edge_nodes( &
         10, 2, nodes, 3, nodes1, nodes2, nodes3)
    case_success = ( &
         all(nodes1 == nodes(:, [1, 2, 3, 4])) .AND. &
         all(nodes2 == nodes(:, [4, 7, 9, 10])) .AND. &
         all(nodes3 == nodes(:, [10, 8, 5, 1])))
    call print_status(name, case_id, case_success, success)

  end subroutine test_compute_edge_nodes

  subroutine test_shoelace_for_area(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(17) :: name
    real(c_double) :: nodes1(2, 2), nodes2(2, 3), nodes3(2, 4)
    real(c_double) :: nodes4(2, 5), nodes5(2, 6)
    real(c_double) :: shoelace
    logical(c_bool) :: not_implemented

    case_id = 1
    name = "shoelace_for_area"

    ! CASE 1: Linear edge.
    nodes1(:, 1) = 1
    nodes1(:, 2) = [2.0_dp, 4.0_dp]
    call shoelace_for_area(2, nodes1, shoelace, not_implemented)
    case_success = ( &
         shoelace == 1.0_dp .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic edge.
    nodes2(:, 1) = 0
    nodes2(:, 2) = [1.0_dp, 3.0_dp]
    nodes2(:, 3) = [4.0_dp, 3.0_dp]
    call shoelace_for_area(3, nodes2, shoelace, not_implemented)
    case_success = ( &
         shoelace == -3.0_dp .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Cubic edge.
    nodes3(:, 1) = 0
    nodes3(:, 2) = [1.0_dp, 9.0_dp]
    nodes3(:, 3) = [2.0_dp, -6.0_dp]
    nodes3(:, 4) = [3.0_dp, 3.0_dp]
    call shoelace_for_area(4, nodes3, shoelace, not_implemented)
    case_success = ( &
         shoelace == 0.0_dp .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Quartic edge.
    nodes4(:, 1) = 0
    nodes4(:, 2) = [1.0_dp, -4.0_dp]
    nodes4(:, 3) = [2.0_dp, 3.0_dp]
    nodes4(:, 4) = [3.0_dp, -4.0_dp]
    nodes4(:, 5) = [4.0_dp, 0.0_dp]
    call shoelace_for_area(5, nodes4, shoelace, not_implemented)
    case_success = ( &
         shoelace == 4.0_dp .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Unsupported degree.
    nodes5 = 1
    call shoelace_for_area(6, nodes5, shoelace, not_implemented)
    case_success = not_implemented
    call print_status(name, case_id, case_success, success)

  end subroutine test_shoelace_for_area

  subroutine test_compute_area(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    character(12) :: name
    type(c_ptr) :: nodes_pointers(4)
    real(c_double) :: area
    logical(c_bool) :: not_implemented
    real(c_double), target :: linear1(2, 2), linear2(2, 2), linear3(2, 2)
    real(c_double), target :: &
         quadratic1(2, 3), quadratic2(2, 3), quadratic3(2, 3)
    real(c_double), target :: quintic1(2, 6)

    case_id = 1
    name = "compute_area"

    ! CASE 1: Curved clam shell (edges are mixed degree).
    linear1(:, 1) = 0
    linear1(:, 2) = [2.0_dp, 0.0_dp]
    nodes_pointers(1) = c_loc(linear1)
    quadratic1(:, 1) = [2.0_dp, 0.0_dp]
    quadratic1(:, 2) = [1.0_dp, 3.0_dp]
    quadratic1(:, 3) = 0
    nodes_pointers(2) = c_loc(quadratic1)
    call compute_area( &
         2, [2, 3], nodes_pointers(:2), area, not_implemented)
    case_success = ( &
         area == 2.0_dp .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Curved triangle (all edges quadratic).
    quadratic1(:, 1) = 0
    quadratic1(:, 2) = [1.0_dp, -1.0_dp]
    quadratic1(:, 3) = [2.0_dp, 0.0_dp]
    nodes_pointers(1) = c_loc(quadratic1)
    quadratic2(:, 1) = [2.0_dp, 0.0_dp]
    quadratic2(:, 2) = [1.5_dp, 1.5_dp]
    quadratic2(:, 3) = [0.0_dp, 2.0_dp]
    nodes_pointers(2) = c_loc(quadratic2)
    quadratic3(:, 1) = [0.0_dp, 2.0_dp]
    quadratic3(:, 2) = 1
    quadratic3(:, 3) = 0
    nodes_pointers(3) = c_loc(quadratic3)
    call compute_area( &
         3, [3, 3, 3], nodes_pointers(:3), area, not_implemented)
    case_success = ( &
         3.0_dp * area == 8.0_dp .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Curved quadrilateral (edges are mixed degree).
    linear1(:, 1) = 0
    linear1(:, 2) = [1.0_dp, 0.0_dp]
    nodes_pointers(1) = c_loc(linear1)
    linear2(:, 1) = [1.0_dp, 0.0_dp]
    linear2(:, 2) = 1
    nodes_pointers(2) = c_loc(linear2)
    quadratic1(:, 1) = 1
    quadratic1(:, 2) = [0.5_dp, 2.0_dp]
    quadratic1(:, 3) = [0.0_dp, 1.0_dp]
    nodes_pointers(3) = c_loc(quadratic1)
    linear3(:, 1) = [0.0_dp, 1.0_dp]
    linear3(:, 2) = 0
    nodes_pointers(4) = c_loc(linear3)
    call compute_area( &
         4, [2, 2, 3, 2], nodes_pointers, area, not_implemented)
    case_success = ( &
         3.0_dp * area == 4.0_dp .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Unsupported degree.
    quintic1 = 1
    nodes_pointers(1) = c_loc(quintic1)
    call compute_area( &
         1, [6], nodes_pointers(:1), area, not_implemented)
    case_success = not_implemented
    call print_status(name, case_id, case_success, success)

  end subroutine test_compute_area

end module test_surface
