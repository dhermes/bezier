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

module test_curve

  use, intrinsic :: iso_c_binding, only: c_bool, c_double, c_int
  use curve, only: &
       CurveData, evaluate_curve_barycentric, evaluate_multi, &
       specialize_curve, evaluate_hodograph, subdivide_nodes, newton_refine, &
       LOCATE_MISS, LOCATE_INVALID, locate_point, elevate_nodes, &
       get_curvature, reduce_pseudo_inverse, full_reduce, compute_length, &
       curves_equal, subdivide_curve
  use types, only: dp
  use unit_test_helpers, only: &
       MACHINE_EPS, print_status, get_random_nodes, get_id_mat
  implicit none
  private &
       test_evaluate_curve_barycentric, test_evaluate_multi, &
       test_specialize_curve, test_evaluate_hodograph, test_subdivide_nodes, &
       subdivide_points_check, test_newton_refine, test_locate_point, &
       test_elevate_nodes, test_get_curvature, test_reduce_pseudo_inverse, &
       pseudo_inverse_helper, test_full_reduce, test_compute_length, &
       test_curves_equal, test_subdivide_curve
  public curve_all_tests

contains

  subroutine curve_all_tests(success)
    logical(c_bool), intent(inout) :: success

    call test_evaluate_curve_barycentric(success)
    call test_evaluate_multi(success)
    call test_specialize_curve(success)
    call test_evaluate_hodograph(success)
    call test_subdivide_nodes(success)
    call test_newton_refine(success)
    call test_locate_point(success)
    call test_elevate_nodes(success)
    call test_get_curvature(success)
    call test_reduce_pseudo_inverse(success)
    call test_full_reduce(success)
    call test_compute_length(success)
    call test_curves_equal(success)
    call test_subdivide_curve(success)

  end subroutine curve_all_tests

  subroutine test_evaluate_curve_barycentric(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes(4, 3)
    real(c_double) :: lambda1(3), lambda2(3)
    real(c_double) :: evaluated(3, 3), expected(3, 3)
    integer :: case_id
    character(26) :: name

    case_id = 1
    name = "evaluate_curve_barycentric"

    ! CASE 1: Slightly complex example.
    nodes(1, :) = [0.0_dp, 0.0_dp, 0.0_dp]
    nodes(2, :) = [0.5_dp, 3.0_dp, 1.0_dp]
    nodes(3, :) = [1.5_dp, 4.0_dp, 1.0_dp]
    nodes(4, :) = [2.0_dp, 8.0_dp, 1.0_dp]
    lambda1 = [0.25_dp, 0.5_dp, 0.75_dp]
    lambda2 = [0.25_dp, 0.125_dp, -0.75_dp]
    expected(1, :) = [0.125_dp, 0.453125_dp, 0.109375_dp]
    expected(2, :) = [0.0859375_dp, 0.390625_dp, 0.119140625_dp]
    expected(3, :) = [0.421875_dp, -2.109375_dp, -0.421875_dp]
    call evaluate_curve_barycentric( &
         4, 3, nodes, 3, lambda1, lambda2, evaluated)

    case_success = all(evaluated == expected)
    call print_status(name, case_id, case_success, success)

  end subroutine test_evaluate_curve_barycentric

  subroutine test_evaluate_multi(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 3), nodes2(3, 2)
    real(c_double) :: s_vals1(129), s_vals2(65)
    real(c_double) :: evaluated1(129, 3), expected1(129, 3)
    real(c_double) :: evaluated2(65, 2), expected2(65, 2)
    integer :: i
    integer :: case_id
    character(14) :: name

    case_id = 1
    name = "evaluate_multi"

    ! CASE 1: Linear curve.
    do i = 1, 129
       s_vals1(i) = (i - 1) / 128.0_dp
    end do
    nodes1(1, :) = [1.0_dp, 1.0_dp, -7.0_dp]
    nodes1(2, :) = [2.0_dp, -1.0_dp, -4.0_dp]
    ! B(s) = [s + 1, 1 - 2 s, 3 s - 7]
    expected1(:, 1) = 1.0_dp + s_vals1
    expected1(:, 2) = 1.0_dp - 2.0_dp * s_vals1
    expected1(:, 3) = -7.0_dp + 3.0_dp * s_vals1
    call evaluate_multi( &
         2, 3, nodes1, 129, s_vals1, evaluated1)
    case_success = all(evaluated1 == expected1)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic curve.
    do i = 1, 65
       s_vals2(i) = (i - 1) / 64.0_dp
    end do
    nodes2(1, :) = [0.0_dp, 0.0_dp]
    nodes2(2, :) = [2.0_dp, -1.0_dp]
    nodes2(3, :) = [3.0_dp, 2.0_dp]
    ! B(s) = [s(4 - s), 2s(2s - 1)]
    expected2(:, 1) = s_vals2 * (4.0_dp - s_vals2)
    expected2(:, 2) = 2.0_dp * s_vals2 * (2.0_dp * s_vals2 - 1.0_dp)
    call evaluate_multi( &
         3, 2, nodes2, 65, s_vals2, evaluated2)
    case_success = all(evaluated2 == expected2)
    call print_status(name, case_id, case_success, success)

  end subroutine test_evaluate_multi

  subroutine test_specialize_curve(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: true_start, true_end
    real(c_double) :: nodes1(2, 2)
    real(c_double) :: new_nodes1(2, 2), expected1(2, 2)
    real(c_double) :: nodes2(3, 2), left_nodes(3, 2), right_nodes(3, 2)
    real(c_double) :: new_nodes2(3, 2)
    real(c_double) :: nodes3(4, 2), new_nodes3(4, 2), expected3(4, 2)
    real(c_double) :: nodes4(5, 2), new_nodes4(5, 2), expected4(5, 2)
    integer :: case_id
    character(16) :: name

    case_id = 1
    name = "specialize_curve"

    ! CASE 1: Linear curve.
    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 1.0_dp]
    expected1(1, :) = [0.25_dp, 0.25_dp]
    expected1(2, :) = [0.75_dp, 0.75_dp]
    call specialize_curve( &
         2, 2, nodes1, 0.25_dp, 0.75_dp, 0.125_dp, 0.25_dp, &
         new_nodes1, true_start, true_end)
    case_success = ( &
         all(new_nodes1 == expected1) .AND. &
         true_start == 0.15625_dp .AND. true_end == 0.21875_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic curve, after subdivision.
    nodes2(1, :) = [0.0_dp, 1.0_dp]
    nodes2(2, :) = [1.0_dp, 6.0_dp]
    nodes2(3, :) = [3.0_dp, 5.0_dp]
    call subdivide_nodes( &
         3, 2, nodes2, left_nodes, right_nodes)
    call specialize_curve( &
         3, 2, nodes2, 0.0_dp, 0.5_dp, 0.0_dp, 1.0_dp, &
         new_nodes2, true_start, true_end)
    case_success = ( &
         all(new_nodes2 == left_nodes) .AND. &
         true_start == 0.0_dp .AND. true_end == 0.5_dp)
    ! Do a "second" check for the right-hand nodes.
    call specialize_curve( &
         3, 2, nodes2, 0.5_dp, 1.0_dp, 0.0_dp, 1.0_dp, &
         new_nodes2, true_start, true_end)
    case_success = ( &
         case_success .AND. &
         all(new_nodes2 == right_nodes) .AND. &
         true_start == 0.5_dp .AND. true_end == 1.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Cubic curve.
    nodes3(1, :) = [0.0_dp, 0.0_dp]
    nodes3(2, :) = [1.0_dp, -1.0_dp]
    nodes3(3, :) = [1.0_dp, -2.0_dp]
    nodes3(4, :) = [3.0_dp, 2.0_dp]
    expected3(1, :) = [171.0_dp, -187.0_dp] / 512.0_dp
    expected3(2, :) = [375.0_dp, -423.0_dp] / 512.0_dp
    expected3(3, :) = [499.0_dp, -579.0_dp] / 512.0_dp
    expected3(4, :) = [735.0_dp, -335.0_dp] / 512.0_dp
    call specialize_curve( &
         4, 2, nodes3, 0.125_dp, 0.625_dp, 0.0_dp, 1.0_dp, &
         new_nodes3, true_start, true_end)
    case_success = ( &
         all(new_nodes3 == expected3) .AND. &
         true_start == 0.125_dp .AND. true_end == 0.625_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Quartic curve.
    nodes4(1, :) = [0.0_dp, 5.0_dp]
    nodes4(2, :) = [1.0_dp, 6.0_dp]
    nodes4(3, :) = [1.0_dp, 7.0_dp]
    nodes4(4, :) = [3.0_dp, 6.0_dp]
    nodes4(5, :) = [3.0_dp, 7.0_dp]
    expected4(1, :) = [1.5625_dp, 6.375_dp]
    expected4(2, :) = [1.78125_dp, 6.4375_dp]
    expected4(3, :) = [2.015625_dp, 6.46875_dp]
    expected4(4, :) = [2.2578125_dp, 6.484375_dp]
    expected4(5, :) = [2.47265625_dp, 6.5234375_dp]
    call specialize_curve( &
         5, 2, nodes4, 0.5_dp, 0.75_dp, 0.0_dp, 1.0_dp, &
         new_nodes4, true_start, true_end)

    case_success = ( &
         all(new_nodes4 == expected4) .AND. &
         true_start == 0.5_dp .AND. true_end == 0.75_dp)
    call print_status(name, case_id, case_success, success)

  end subroutine test_specialize_curve

  subroutine test_evaluate_hodograph(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: s_vals(5), s_val
    real(c_double) :: hodograph(1, 2), expected(1, 2)
    real(c_double) :: nodes1(2, 2), nodes2(3, 2), nodes3(4, 2)
    integer :: i
    integer :: case_id
    character(18) :: name

    case_id = 1
    name = "evaluate_hodograph"

    ! CASE 1: Linear curve.
    nodes1(1, :) = [0.0_dp, 0.0_dp]
    nodes1(2, :) = [1.0_dp, 1.0_dp]
    expected(1, :) = nodes1(2, :) - nodes1(1, :)
    s_vals(:2) = [0.25_dp, 0.75_dp]
    case_success = .TRUE.
    do i = 1, 2
       call evaluate_hodograph( &
            s_vals(i), 2, 2, nodes1, hodograph)
       case_success = ( &
            case_success .AND. &
            all(hodograph == expected))
    end do
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic curve.
    nodes2(1, :) = [0.0_dp, 0.0_dp]
    nodes2(2, :) = [0.5_dp, 1.0_dp]
    nodes2(3, :) = [1.25_dp, 0.25_dp]
    !  B(s) = [s(s + 4)/4, s(8 - 7s)/4]
    ! B'(s) = [(2 + s)/2, (4 - 7s)/2]
    s_vals = [0.0_dp, 0.25_dp, 0.5_dp, 0.625_dp, 0.875_dp]
    case_success = .TRUE.
    do i = 1, 5
       s_val = s_vals(i)
       call evaluate_hodograph( &
            s_val, 3, 2, nodes2, hodograph)
       expected(1, :) = [ &
            (2.0_dp + s_val) / 2.0_dp, &
            (4.0_dp - 7.0_dp * s_val) / 2.0_dp]
       case_success = ( &
            case_success .AND. &
            all(hodograph == expected))
    end do
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Cubic curve.
    nodes3(1, :) = [0.0_dp, 0.0_dp]
    nodes3(2, :) = [0.25_dp, 1.0_dp]
    nodes3(3, :) = [0.75_dp, 0.5_dp]
    nodes3(4, :) = [1.25_dp, 1.0_dp]
    !  B(s) = [s(3 + 3s - s^2)/4, s(5s^2 - 9s + 6)/2]
    ! B'(s) = [3(1 + 2s - s^2)/4, 3(5s^2 - 6s + 2)/2]
    s_vals = [0.125_dp, 0.5_dp, 0.75_dp, 1.0_dp, 1.125_dp]
    case_success = .TRUE.
    do i = 1, 5
       s_val = s_vals(i)
       call evaluate_hodograph( &
            s_val, 4, 2, nodes3, hodograph)
       expected(1, 1) = ( &
            3.0_dp * (1.0_dp + 2.0_dp * s_val - s_val * s_val) / 4.0_dp)
       expected(1, 2) = ( &
            1.5_dp * (5.0_dp * s_val * s_val - 6.0_dp * s_val + 2.0_dp))
       case_success = ( &
            case_success .AND. &
            all(hodograph == expected))
    end do
    call print_status(name, case_id, case_success, success)

  end subroutine test_evaluate_hodograph

  subroutine test_subdivide_nodes(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 2), l_expected1(2, 2), r_expected1(2, 2)
    real(c_double) :: left_nodes1(2, 2), right_nodes1(2, 2)
    real(c_double) :: nodes2(3, 2), l_expected2(3, 2), r_expected2(3, 2)
    real(c_double) :: left_nodes2(3, 2), right_nodes2(3, 2)
    real(c_double) :: nodes3(4, 2), l_expected3(4, 2), r_expected3(4, 2)
    real(c_double) :: left_nodes3(4, 2), right_nodes3(4, 2)
    integer :: case_id
    character(23) :: name

    case_id = 1
    name = "subdivide_nodes (Curve)"

    ! CASE 1: Linear curve.
    nodes1(1, :) = [0.0_dp, 1.0_dp]
    nodes1(2, :) = [4.0_dp, 6.0_dp]
    l_expected1(1, :) = [0.0_dp, 1.0_dp]
    l_expected1(2, :) = [2.0_dp, 3.5_dp]
    r_expected1(1, :) = [2.0_dp, 3.5_dp]
    r_expected1(2, :) = [4.0_dp, 6.0_dp]
    call subdivide_nodes( &
         2, 2, nodes1, left_nodes1, right_nodes1)

    case_success = ( &
         all(left_nodes1 == l_expected1) .AND. &
         all(right_nodes1 == r_expected1))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Evaluate subdivided parts of a line.
    call subdivide_points_check( &
         2, 2, 909203, 10489, case_success)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Quadratic curve.
    nodes2(1, :) = [0.0_dp, 1.0_dp]
    nodes2(2, :) = [4.0_dp, 6.0_dp]
    nodes2(3, :) = [7.0_dp, 3.0_dp]
    l_expected2(1, :) = [0.0_dp, 1.0_dp]
    l_expected2(2, :) = [2.0_dp, 3.5_dp]
    l_expected2(3, :) = [3.75_dp, 4.0_dp]
    r_expected2(1, :) = [3.75_dp, 4.0_dp]
    r_expected2(2, :) = [5.5_dp, 4.5_dp]
    r_expected2(3, :) = [7.0_dp, 3.0_dp]
    call subdivide_nodes( &
         3, 2, nodes2, left_nodes2, right_nodes2)

    case_success = ( &
         all(left_nodes2 == l_expected2) .AND. &
         all(right_nodes2 == r_expected2))
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Evaluate subdivided parts of a quadratic curve.
    call subdivide_points_check( &
         3, 2, 10920388, 7756, case_success)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Cubic curve.
    nodes3(1, :) = [0.0_dp, 1.0_dp]
    nodes3(2, :) = [4.0_dp, 6.0_dp]
    nodes3(3, :) = [7.0_dp, 3.0_dp]
    nodes3(4, :) = [6.0_dp, 5.0_dp]
    l_expected3(1, :) = [0.0_dp, 1.0_dp]
    l_expected3(2, :) = [2.0_dp, 3.5_dp]
    l_expected3(3, :) = [3.75_dp, 4.0_dp]
    l_expected3(4, :) = [4.875_dp, 4.125_dp]
    r_expected3(1, :) = [4.875_dp, 4.125_dp]
    r_expected3(2, :) = [6.0_dp, 4.25_dp]
    r_expected3(3, :) = [6.5_dp, 4.0_dp]
    r_expected3(4, :) = [6.0_dp, 5.0_dp]
    call subdivide_nodes( &
         4, 2, nodes3, left_nodes3, right_nodes3)

    case_success = ( &
         all(left_nodes3 == l_expected3) .AND. &
         all(right_nodes3 == r_expected3))
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Evaluate subdivided parts of a cubic curve.
    call subdivide_points_check( &
         4, 2, 33981020, 81902302, case_success)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Quartic curve (not hardcoded).
    call subdivide_points_check( &
         5, 2, 1190, 881, case_success)
    call print_status(name, case_id, case_success, success)

  end subroutine test_subdivide_nodes

  subroutine subdivide_points_check( &
       num_nodes, dimension_, multiplier, modulus, success)
    integer, intent(in) :: num_nodes, dimension_, multiplier, modulus
    logical, intent(out) :: success
    ! Variables outside of signature.
    real(c_double) :: nodes(num_nodes, dimension_)
    real(c_double) :: left(num_nodes, dimension_)
    real(c_double) :: right(num_nodes, dimension_)
    real(c_double) :: unit_interval(33), left_s(33), right_s(33)
    real(c_double) :: evaluated1(33, dimension_)
    real(c_double) :: evaluated2(33, dimension_)
    integer :: i

    success = .TRUE.

    call get_random_nodes(nodes, multiplier, modulus, num_bits=8)
    call subdivide_nodes( &
         num_nodes, dimension_, nodes, left, right)

    do i = 1, 33
       unit_interval(i) = (i - 1) / 32.0_dp
       left_s(i) = (i - 1) / 64.0_dp
       right_s(i) = (i + 31) / 64.0_dp
    end do

    call evaluate_multi( &
         num_nodes, dimension_, nodes, 33, left_s, evaluated1)
    call evaluate_multi( &
         num_nodes, dimension_, left, 33, unit_interval, evaluated2)
    success = ( &
         success .AND. &
         all(evaluated1 == evaluated2))

    call evaluate_multi( &
         num_nodes, dimension_, nodes, 33, right_s, evaluated1)
    call evaluate_multi( &
         num_nodes, dimension_, right, 33, unit_interval, evaluated2)
    success = ( &
         success .AND. &
         all(evaluated1 == evaluated2))

  end subroutine subdivide_points_check

  subroutine test_newton_refine(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes(4, 3)
    real(c_double) :: point(1, 3)
    real(c_double) :: new_s
    integer :: case_id
    character(13) :: name

    case_id = 1
    name = "newton_refine"

    ! CASE 1: Cubic in 3D.
    nodes(1, :) = [0.0_dp, 0.0_dp, 0.0_dp]
    nodes(2, :) = [1.0_dp, -1.0_dp, 1.0_dp]
    nodes(3, :) = [3.0_dp, 2.0_dp, 2.0_dp]
    nodes(4, :) = [2.0_dp, 2.0_dp, 4.0_dp]
    ! curve(1/2) = p
    point(1, :) = [1.75_dp, 0.625_dp, 1.625_dp]
    call newton_refine( &
         4, 3, nodes, point, 0.25_dp, new_s)
    case_success = (110.0_dp * new_s == 57.0_dp)
    call print_status(name, case_id, case_success, success)

  end subroutine test_newton_refine

  subroutine test_locate_point(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(3, 3), point1(1, 3)
    real(c_double) :: nodes2(3, 2), point2(1, 2)
    real(c_double) :: nodes3(4, 2), point3(1, 2)
    real(c_double) :: s_val
    integer :: case_id
    character(20) :: name

    case_id = 1
    name = "locate_point (Curve)"

    ! CASE 1: Quadratic in 3D (with match).
    nodes1(1, :) = 0
    nodes1(2, :) = [3.0_dp, 0.0_dp, -1.0_dp]
    nodes1(3, :) = [1.0_dp, 1.0_dp, 3.0_dp]
    point1(1, :) = [43.0_dp, 1.0_dp, -11.0_dp] / 64.0_dp
    call locate_point( &
         3, 3, nodes1, point1, s_val)

    case_success = (s_val == 0.125_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic in 2D (without match).
    nodes2(1, :) = 0
    nodes2(2, :) = [0.5_dp, 1.0_dp]
    nodes2(3, :) = [1.0_dp, 0.0_dp]
    point2(1, :) = [0.5_dp, 2.0_dp]
    call locate_point( &
         3, 2, nodes2, point2, s_val)

    case_success = (s_val == LOCATE_MISS)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Self-intersecting cubic in 3D (will cause failure).
    nodes3(1, :) = [0.0_dp, 2.0_dp]
    nodes3(2, :) = [-1.0_dp, 0.0_dp]
    nodes3(3, :) = [1.0_dp, 1.0_dp]
    nodes3(4, :) = [-0.75_dp, 1.625_dp]
    point3(1, :) = [-0.25_dp, 1.375_dp]
    call locate_point( &
         4, 2, nodes3, point3, s_val)

    case_success = (s_val == LOCATE_INVALID)
    call print_status(name, case_id, case_success, success)

  end subroutine test_locate_point

  subroutine test_elevate_nodes(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 2), elevated1(3, 2), expected1(3, 2)
    real(c_double) :: nodes2(3, 3), elevated2(4, 3), expected2(4, 3)
    integer :: case_id
    character(13) :: name

    case_id = 1
    name = "elevate_nodes"

    ! CASE 1: Linear curve.
    nodes1(1, :) = 0
    nodes1(2, :) = [2.0_dp, 4.0_dp]
    expected1(1, :) = [0.0_dp, 0.0_dp]
    expected1(2, :) = [1.0_dp, 2.0_dp]
    expected1(3, :) = [2.0_dp, 4.0_dp]
    call elevate_nodes( &
         2, 2, nodes1, elevated1)
    case_success = all(elevated1 == expected1)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Quadratic curve.
    nodes2(1, :) = [0.0_dp, 0.5_dp, 0.75_dp]
    nodes2(2, :) = [3.0_dp, 0.5_dp, 3.0_dp]
    nodes2(3, :) = [6.0_dp, 0.5_dp, 2.25_dp]
    expected2(1, :) = [0.0_dp, 0.5_dp, 0.75_dp]
    expected2(2, :) = [2.0_dp, 0.5_dp, 2.25_dp]
    expected2(3, :) = [4.0_dp, 0.5_dp, 2.75_dp]
    expected2(4, :) = [6.0_dp, 0.5_dp, 2.25_dp]
    call elevate_nodes( &
         3, 3, nodes2, elevated2)
    case_success = all(elevated2 == expected2)
    call print_status(name, case_id, case_success, success)

  end subroutine test_elevate_nodes

  subroutine test_get_curvature(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 2), nodes2(3, 2)
    real(c_double) :: s_val, tangent_vec(1, 2), curvature
    integer :: case_id
    character(13) :: name

    case_id = 1
    name = "get_curvature"

    ! CASE 1: Linear curve.
    nodes1(1, :) = 0
    nodes1(2, :) = 1
    s_val = 0.5_dp
    call evaluate_hodograph( &
         s_val, 2, 2, nodes1, tangent_vec)
    call get_curvature( &
         2, 2, nodes1, tangent_vec, s_val, curvature)
    case_success = (curvature == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Linear curve, degree-elevated to quadratic.
    nodes2(1, :) = 0
    nodes2(2, :) = 0.5_dp
    nodes2(3, :) = 1
    s_val = 0.25_dp
    call evaluate_hodograph( &
         s_val, 3, 2, nodes2, tangent_vec)
    call get_curvature( &
         3, 2, nodes2, tangent_vec, s_val, curvature)
    case_success = (curvature == 0.0_dp)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Quadratic curve.
    nodes2(1, :) = 0
    nodes2(2, :) = [0.5_dp, 1.0_dp]
    nodes2(3, :) = [1.0_dp, 0.0_dp]
    s_val = 0.5_dp
    call evaluate_hodograph( &
         s_val, 3, 2, nodes2, tangent_vec)
    call get_curvature( &
         3, 2, nodes2, tangent_vec, s_val, curvature)
    case_success = (curvature == -4.0_dp)
    call print_status(name, case_id, case_success, success)

  end subroutine test_get_curvature

  subroutine test_reduce_pseudo_inverse(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    logical(c_bool) :: not_implemented
    real(c_double) :: nodes1(2, 2), reduced1(1, 2), expected1(1, 2)
    real(c_double) :: nodes2(3, 2), reduced2(2, 2), expected2(2, 2)
    real(c_double) :: nodes3(4, 3), reduced3(3, 3), expected3(3, 3)
    real(c_double) :: nodes4(5, 1), reduced4(4, 1), expected4(4, 1)
    real(c_double) :: reduced5(4, 4), expected5(4, 4)
    real(c_double) :: nodes6(6, 2), reduced6(5, 2)
    integer :: case_id
    character(21) :: name

    case_id = 1
    name = "reduce_pseudo_inverse"

    ! CASE 1: Reduce from line to point/constant.
    nodes1(1, :) = [-2.0_dp, 1.0_dp]
    nodes1(2, :) = nodes1(1, :)
    expected1(1, :) = nodes1(1, :)
    call reduce_pseudo_inverse( &
         2, 2, nodes1, reduced1, not_implemented)
    case_success = (all(reduced1 == expected1) .AND. .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Reduce a degree-elevated line (as quadratic).
    nodes2(1, :) = 0
    nodes2(2, :) = [1.0_dp, 2.0_dp]
    nodes2(3, :) = [2.0_dp, 4.0_dp]
    expected2(1, :) = 0
    expected2(2, :) = [2.0_dp, 4.0_dp]
    call reduce_pseudo_inverse( &
         3, 2, nodes2, reduced2, not_implemented)
    case_success = (all(reduced2 == expected2) .AND. .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Test that degree 2->1 reduction is actually a pseudo-inverse.
    call pseudo_inverse_helper(2, reduced2, expected2, not_implemented)
    case_success = (all(reduced2 == expected2) .AND. .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Reduce a quadratic that is not a deg.-elevated line.
    nodes2(1, :) = 0
    nodes2(2, :) = [1.0_dp, 1.5_dp]
    nodes2(3, :) = [2.0_dp, 0.0_dp]
    expected2(1, :) = [0.0_dp, 0.5_dp]
    expected2(2, :) = [2.0_dp, 0.5_dp]
    call reduce_pseudo_inverse( &
         3, 2, nodes2, reduced2, not_implemented)
    case_success = (all(reduced2 == expected2) .AND. .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Reduce a degree-elevated quadratic (as cubic).
    nodes3(1, :) = [0.0_dp, 0.5_dp, 0.75_dp]
    nodes3(2, :) = [2.0_dp, 0.5_dp, 2.25_dp]
    nodes3(3, :) = [4.0_dp, 0.5_dp, 2.75_dp]
    nodes3(4, :) = [6.0_dp, 0.5_dp, 2.25_dp]
    expected3(1, :) = nodes3(1, :)
    expected3(2, :) = [3.0_dp, 0.5_dp, 3.0_dp]
    expected3(3, :) = nodes3(4, :)
    call reduce_pseudo_inverse( &
         4, 3, nodes3, reduced3, not_implemented)
    case_success = (all(reduced3 == expected3) .AND. .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Test that degree 3->2 reduction is actually a pseudo-inverse.
    call pseudo_inverse_helper(3, reduced3, expected3, not_implemented)
    case_success = ( &
         maxval(abs(reduced3 - expected3)) < MACHINE_EPS .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Reduce a degree-elevated cubic (as quartic).
    nodes4(:, 1) = [0.0_dp, 0.75_dp, 2.0_dp, 2.75_dp, 2.0_dp]
    expected4(:, 1) = [0.0_dp, 1.0_dp, 3.0_dp, 2.0_dp]
    call reduce_pseudo_inverse( &
         5, 1, nodes4, reduced4, not_implemented)
    case_success = (all(reduced4 == expected4) .AND. .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Test that degree 4->3 reduction is actually a pseudo-inverse.
    call pseudo_inverse_helper(4, reduced5, expected5, not_implemented)
    case_success = ( &
         maxval(abs(reduced3 - expected3)) < MACHINE_EPS .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 9: Unsupported degree
    call get_random_nodes(nodes6, 9177093, 16718, num_bits=8)
    call reduce_pseudo_inverse( &
         6, 2, nodes6, reduced6, not_implemented)
    case_success = (not_implemented)
    call print_status(name, case_id, case_success, success)

  end subroutine test_reduce_pseudo_inverse

  subroutine pseudo_inverse_helper(degree, result_, id_mat, not_implemented)
    integer, intent(in) :: degree
    real(c_double), intent(out) :: result_(degree, degree)
    real(c_double), intent(out) :: id_mat(degree, degree)
    logical(c_bool), intent(out) :: not_implemented
    ! Variables outside of signature.
    real(c_double) :: nodes(degree + 1, degree + 1)
    real(c_double) :: reduction_mat(degree, degree + 1)
    real(c_double) :: elevation_mat(degree + 1, degree)

    nodes = get_id_mat(degree + 1)
    call reduce_pseudo_inverse( &
         degree + 1, degree + 1, nodes, reduction_mat, not_implemented)
    ! NOTE: We assume none of the callers will use inputs that force
    !       ``not_implemeted`` to be ``.FALSE.``.

    id_mat = get_id_mat(degree)
    call elevate_nodes( &
         degree, degree, id_mat, elevation_mat)

    result_ = matmul(reduction_mat, elevation_mat)

  end subroutine pseudo_inverse_helper

  subroutine test_full_reduce(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 1), expected1(2, 1), reduced1(2, 1)
    real(c_double) :: nodes2(4, 2), expected2(4, 2), reduced2(4, 2)
    real(c_double) :: nodes3(6, 2), reduced3(6, 2)
    real(c_double) :: nodes4(1, 2), expected4(1, 2), reduced4(1, 2)
    real(c_double) :: nodes5(5, 2), expected5(5, 2), reduced5(5, 2)
    integer(c_int) :: num_reduced_nodes
    logical(c_bool) :: not_implemented
    integer :: case_id
    character(11) :: name

    case_id = 1
    name = "full_reduce"

    ! CASE 1: Linear curve, one reduction.
    nodes1 = 5.5_dp
    expected1(1, :) = 5.5_dp
    call full_reduce( &
         2, 1, nodes1, num_reduced_nodes, &
         reduced1, not_implemented)
    case_success = ( &
         num_reduced_nodes == 1 .AND. &
         all(reduced1(:1, :) == expected1(:1, :)) .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Cubic curve, one reduction.
    nodes2(1, :) = 0
    nodes2(2, :) = [2.0_dp, 4.0_dp]
    nodes2(3, :) = [4.0_dp, 6.0_dp]
    nodes2(4, :) = [6.0_dp, 6.0_dp]
    expected2(1, :) = 0
    expected2(2, :) = [3.0_dp, 6.0_dp]
    expected2(3, :) = [6.0_dp, 6.0_dp]
    call full_reduce( &
         4, 2, nodes2, num_reduced_nodes, &
         reduced2, not_implemented)
    case_success = ( &
         num_reduced_nodes == 3 .AND. &
         all(reduced2(:3, :) == expected2(:3, :)) .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Cubic curve, two reductions.
    nodes2(1, :) = [0.0_dp, 4.0_dp]
    nodes2(2, :) = [1.0_dp, 4.5_dp]
    nodes2(3, :) = [2.0_dp, 5.0_dp]
    nodes2(4, :) = [3.0_dp, 5.5_dp]
    expected2(1, :) = [0.0_dp, 4.0_dp]
    expected2(2, :) = [3.0_dp, 5.5_dp]
    call full_reduce( &
         4, 2, nodes2, num_reduced_nodes, &
         reduced2, not_implemented)
    case_success = ( &
         num_reduced_nodes == 2 .AND. &
         all(reduced2(:2, :) == expected2(:2, :)) .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Cubic curve, no reductions.
    nodes2(1, :) = [0.0_dp, 2.0_dp]
    nodes2(2, :) = [-1.0_dp, 0.0_dp]
    nodes2(3, :) = 1
    nodes2(4, :) = [-0.75_dp, 1.625_dp]
    expected2 = nodes2
    call full_reduce( &
         4, 2, nodes2, num_reduced_nodes, &
         reduced2, not_implemented)
    case_success = ( &
         num_reduced_nodes == 4 .AND. &
         all(reduced2(:4, :) == expected2(:4, :)) .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Unsupported degree (quintic, degree 5).
    nodes3 = 0
    call full_reduce( &
         6, 2, nodes3, num_reduced_nodes, &
         reduced3, not_implemented)
    case_success = (not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Point/constant (can't be any reduction).
    nodes4(1, :) = [0.0_dp, 1.0_dp]
    expected4 = nodes4
    call full_reduce( &
         1, 2, nodes4, num_reduced_nodes, &
         reduced4, not_implemented)
    case_success = ( &
         num_reduced_nodes == 1 .AND. &
         all(reduced4 == expected4) .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Quartic curve, one reduction.
    expected5(1, :) = [0.0_dp, 0.0_dp]
    expected5(2, :) = [2.0_dp, 4.0_dp]
    expected5(3, :) = [4.0_dp, 6.5_dp]
    expected5(4, :) = [18.0_dp, 1.5_dp]
    nodes5(1, :) = expected5(1, :)
    nodes5(2, :) = (expected5(1, :) + 3 * expected5(2, :)) / 4.0_dp
    nodes5(3, :) = (expected5(2, :) + expected5(3, :)) / 2.0_dp
    nodes5(4, :) = (3 * expected5(3, :) + expected5(4, :)) / 4.0_dp
    nodes5(5, :) = expected5(4, :)
    call full_reduce( &
         5, 2, nodes5, num_reduced_nodes, &
         reduced5, not_implemented)
    case_success = ( &
         num_reduced_nodes == 4 .AND. &
         all(reduced5(:4, :) == expected5(:4, :)) .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Quartic curve, no reductions.
    nodes5(1, :) = [0.0_dp, 0.0_dp]
    nodes5(2, :) = [1.0_dp, 2.0_dp]
    nodes5(3, :) = [2.0_dp, 3.5_dp]
    nodes5(4, :) = [5.0_dp, 3.5_dp]
    nodes5(5, :) = [8.0_dp, 18.0_dp]
    call full_reduce( &
         5, 2, nodes5, num_reduced_nodes, &
         reduced5, not_implemented)
    case_success = ( &
         num_reduced_nodes == 5 .AND. &
         all(reduced5 == nodes5) .AND. &
         .NOT. not_implemented)
    call print_status(name, case_id, case_success, success)

  end subroutine test_full_reduce

  subroutine test_compute_length(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    real(c_double) :: nodes1(2, 2), nodes2(3, 2), nodes3(4, 2)
    real(c_double) :: length, expected
    integer(c_int) :: error_val
    integer :: case_id
    character(14) :: name

    case_id = 1
    name = "compute_length"

    ! CASE 1: Linear curve.
    nodes1(1, :) = 0
    nodes1(2, :) = [3.0_dp, 4.0_dp]
    call compute_length( &
         2, 2, nodes1, length, error_val)
    case_success = (length == 5.0_dp .AND. error_val == 0)
    call print_status(name, case_id, case_success, success)

    ! CASE 1: Quadratic curve.
    nodes2(1, :) = 0
    nodes2(2, :) = [1.0_dp, 2.0_dp]
    nodes2(3, :) = [2.0_dp, 0.0_dp]
    call compute_length( &
         3, 2, nodes2, length, error_val)
    ! 2 INT_0^1 SQRT(16 s^2  - 16 s + 5) ds = SQRT(5) + sinh^{-1}(2)/2
    expected = sqrt(5.0_dp) + 0.5_dp * asinh(2.0_dp)
    case_success = ( &
         abs(length - expected) <= spacing(expected) .AND. &
         error_val == 0)
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Cubic curve.
    nodes3(1, :) = 0
    nodes3(2, :) = [1.0_dp, 2.0_dp]
    nodes3(3, :) = [2.0_dp, 0.0_dp]
    nodes3(4, :) = [3.5_dp, 0.0_dp]
    call compute_length( &
         4, 2, nodes3, length, error_val)
    ! x(s) = s (s^2 + 6) / 2
    ! y(s) = 6 s (s - 1)^2
    ! x'(s)^2 + y'(s)^2 = (9/4)(145s^4 - 384s^3 + 356s^2 - 128s + 20)
    ! NOTE: This "literal" has more precision than the type provides.
    expected = 4.0916195514424730685_dp
    case_success = ( &
         abs(length - expected) <= spacing(expected) .AND. &
         error_val == 0)
    call print_status(name, case_id, case_success, success)

  end subroutine test_compute_length

  subroutine test_curves_equal(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    type(CurveData) :: &
         curve1, curve2, curve3, curve4, curve5, curve6, curve7
    character(12) :: name

    case_id = 1
    name = "curves_equal"

    ! curve1 is the basis for all comparisons.
    curve1%start = 0.25_dp
    curve1%end_ = 0.75_dp
    allocate(curve1%nodes(2, 2))
    curve1%nodes(1, :) = [0.5_dp, 2.5_dp]
    curve1%nodes(2, :) = [7.0_dp, -2.0_dp]

    ! CASE 1: Compare curve1 to itself.
    case_success = curves_equal(curve1, curve1)
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Compare curve1 to uninitialized nodes.
    curve2%start = curve1%start
    curve2%end_ = curve1%end_
    case_success = ( &
         .NOT. allocated(curve2%nodes) .AND. &
         .NOT. curves_equal(curve1, curve2) .AND. &
         .NOT. curves_equal(curve2, curve1))
    call print_status(name, case_id, case_success, success)

    ! CASE 3: Compare curve1 to an identical, but different curve.
    curve3%start = curve1%start
    curve3%end_ = curve1%end_
    curve3%nodes = curve1%nodes
    case_success = ( &
         curves_equal(curve1, curve3) .AND. curves_equal(curve3, curve1))
    call print_status(name, case_id, case_success, success)

    ! CASE 4: Compare curve1 to a curve that differs only in ``start``.
    curve4%start = curve1%start - 0.25_dp
    curve4%end_ = curve1%end_
    curve4%nodes = curve1%nodes
    case_success = ( &
         .NOT. curves_equal(curve1, curve4) .AND. &
         .NOT. curves_equal(curve4, curve1))
    call print_status(name, case_id, case_success, success)

    ! CASE 5: Compare curve1 to a curve that differs only in ``end_``.
    curve5%start = curve1%start
    curve5%end_ = curve1%end_ + 0.25_dp
    curve5%nodes = curve1%nodes
    case_success = ( &
         .NOT. curves_equal(curve1, curve5) .AND. &
         .NOT. curves_equal(curve5, curve1))
    call print_status(name, case_id, case_success, success)

    ! CASE 6: Compare curve1 to a curve that differs in node shape
    curve6%start = curve1%start
    curve6%end_ = curve1%end_
    allocate(curve6%nodes(1, 2))
    curve6%nodes(1, :) = curve1%nodes(1, :)
    case_success = ( &
         .NOT. curves_equal(curve1, curve6) .AND. &
         .NOT. curves_equal(curve6, curve1))
    call print_status(name, case_id, case_success, success)

    ! CASE 7: Compare curve1 to a curve that differs in
    !         ``nodes`` value but not shape.
    curve7%start = curve1%start
    curve7%end_ = curve1%end_
    curve7%nodes = curve1%nodes + 2.5_dp
    case_success = ( &
         .NOT. curves_equal(curve1, curve7) .AND. &
         .NOT. curves_equal(curve7, curve1))
    call print_status(name, case_id, case_success, success)

    ! CASE 8: Compare curve with uninitialized nodes to **itself**.
    case_success = ( &
         .NOT. allocated(curve2%nodes) .AND. &
         curves_equal(curve2, curve2))
    call print_status(name, case_id, case_success, success)

  end subroutine test_curves_equal

  subroutine test_subdivide_curve(success)
    logical(c_bool), intent(inout) :: success
    ! Variables outside of signature.
    logical :: case_success
    integer :: case_id
    type(CurveData) :: curve_data, left1, left2, right
    real(c_double) :: left_nodes(4, 2), right_nodes(4, 2)
    character(15) :: name

    case_id = 1
    name = "subdivide_curve"

    ! CASE 1: Curve has not yet been subdivided.
    allocate(curve_data%nodes(4, 2))
    curve_data%nodes(1, :) = [0.0_dp, 1.0_dp]
    curve_data%nodes(2, :) = [1.0_dp, 2.0_dp]
    curve_data%nodes(3, :) = [1.0_dp, 2.0_dp]
    curve_data%nodes(4, :) = [3.0_dp, 0.0_dp]

    left_nodes(1, :) = [0.0_dp, 1.0_dp]
    left_nodes(2, :) = [0.5_dp, 1.5_dp]
    left_nodes(3, :) = [0.75_dp, 1.75_dp]
    left_nodes(4, :) = [1.125_dp, 1.625_dp]
    right_nodes(1, :) = [1.125_dp, 1.625_dp]
    right_nodes(2, :) = [1.5_dp, 1.5_dp]
    right_nodes(3, :) = [2.0_dp, 1.0_dp]
    right_nodes(4, :) = [3.0_dp, 0.0_dp]

    call subdivide_curve(curve_data, left1, right)
    case_success = ( &
         left1%start == 0.0_dp .AND. &
         left1%end_ == 0.5_dp .AND. &
         all(left1%nodes == left_nodes) .AND. &
         right%start == 0.5_dp .AND. &
         right%end_ == 1.0_dp .AND. &
         all(right%nodes == right_nodes))
    call print_status(name, case_id, case_success, success)

    ! CASE 2: Curve has **already** been subdivided.
    left_nodes(1, :) = [0.0_dp, 1.0_dp]
    left_nodes(2, :) = [0.25_dp, 1.25_dp]
    left_nodes(3, :) = [0.4375_dp, 1.4375_dp]
    left_nodes(4, :) = [0.609375_dp, 1.546875_dp]
    right_nodes(1, :) = [0.609375_dp, 1.546875_dp]
    right_nodes(2, :) = [0.78125_dp, 1.65625_dp]
    right_nodes(3, :) = [0.9375_dp, 1.6875_dp]
    right_nodes(4, :) = [1.125_dp, 1.625_dp]

    call subdivide_curve(left1, left2, right)
    case_success = ( &
         left2%start == 0.0_dp .AND. &
         left2%end_ == 0.25_dp .AND. &
         all(left2%nodes == left_nodes) .AND. &
         right%start == 0.25_dp .AND. &
         right%end_ == 0.5_dp .AND. &
         all(right%nodes == right_nodes))
    call print_status(name, case_id, case_success, success)

  end subroutine test_subdivide_curve

end module test_curve
