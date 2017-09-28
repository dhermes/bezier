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

module surface

  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use types, only: dp
  use curve, only: evaluate_curve_barycentric
  implicit none
  private
  public &
       de_casteljau_one_round, evaluate_barycentric, &
       evaluate_barycentric_multi, evaluate_cartesian_multi, jacobian_both, &
       jacobian_det

contains

  subroutine de_casteljau_one_round( &
       num_nodes, dimension_, nodes, degree, &
       lambda1, lambda2, lambda3, new_nodes) &
       bind(c, name='de_casteljau_one_round')

    ! NOTE: This is de Casteljau on a Bezier surface / triangle.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(in) :: lambda1
    real(c_double), intent(in) :: lambda2
    real(c_double), intent(in) :: lambda3
    real(c_double), intent(out) :: &
         new_nodes(num_nodes - degree - 1, dimension_)
    ! Variables outside of signature.
    integer(c_int) :: index_
    integer(c_int) :: parent_i1, parent_i2, parent_i3
    integer(c_int) :: k, j

    index_ = 1
    parent_i1 = 1
    parent_i2 = 2
    parent_i3 = degree + 2
    ! NOTE: Throughout for index_ <--> (i, j, k) in the (degree - 1)
    !       triangle, we have parent_i1 = index_ + k <--> (i + 1, j, k),
    !       parent_i2 = index_ + k + 1 <--> (i, j + 1, k) and
    !       parent_i3 = index_ + degree + 1 <--> (i, j, k + 1).
    do k = 0, degree - 1
       do j = 0, degree - k - 1
          ! NOTE: i = (degree - 1) - j - k
          new_nodes(index_, :) = ( &
               lambda1 * nodes(parent_i1, :) + &
               lambda2 * nodes(parent_i2, :) + &
               lambda3 * nodes(parent_i3, :))
          ! Update all the indices.
          parent_i1 = parent_i1 + 1
          parent_i2 = parent_i2 + 1
          parent_i3 = parent_i3 + 1
          index_ = index_ + 1
       end do
       ! Update the indices that depend on k.
       parent_i1 = parent_i1 + 1
       parent_i2 = parent_i2 + 1
    end do

  end subroutine de_casteljau_one_round

  subroutine evaluate_barycentric( &
       num_nodes, dimension_, nodes, degree, &
       lambda1, lambda2, lambda3, point) &
       bind(c, name='evaluate_barycentric')

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This assumes degree >= 1.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(in) :: lambda1
    real(c_double), intent(in) :: lambda2
    real(c_double), intent(in) :: lambda3
    real(c_double), intent(out) :: point(1, dimension_)
    ! Variables outside of signature.
    real(c_double) :: param_vals(1, 3)

    param_vals(1, 1) = lambda1
    param_vals(1, 2) = lambda2
    param_vals(1, 3) = lambda3
    call evaluate_barycentric_multi( &
         num_nodes, dimension_, nodes, degree, 1, param_vals, point)

  end subroutine evaluate_barycentric

  subroutine evaluate_barycentric_multi( &
       num_nodes, dimension_, nodes, degree, num_vals, param_vals, evaluated) &
       bind(c, name='evaluate_barycentric_multi')

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This assumes degree >= 1.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in) :: degree, num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 3)
    real(c_double), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    integer(c_int) :: k, binom_val, index_, new_index
    real(c_double) :: row_result(num_vals, dimension_)

    index_ = num_nodes
    forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
       evaluated(new_index, :) = nodes(index_, :)
    end forall

    if (degree == 0) then
       return
    end if

    binom_val = 1
    do k = degree - 1, 0, -1
       ! We want to go from (d C (k + 1)) to (d C k).
       binom_val = (binom_val * (k + 1)) / (degree - k)
       index_ = index_ - 1  ! Step to last element in row.
       !     k = d - 1, d - 2, ...
       ! d - k =     1,     2, ...
       ! We know row k has (d - k + 1) elements.
       new_index = index_ - degree + k  ! First element in row.

       ! lambda1 = param_vals(:, 1)
       ! lambda2 = param_vals(:, 2)
       call evaluate_curve_barycentric( &
            degree - k, dimension_, nodes(new_index:index_, :), &
            num_vals, param_vals(:, 1), param_vals(:, 2), row_result)

       ! Update index for next iteration.
       index_ = new_index

       ! lambda3 = param_vals(:, 3)
       forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
          evaluated(new_index, :) = ( &
               param_vals(new_index, 3) * evaluated(new_index, :) + &
               binom_val * row_result(new_index, :))
       end forall
    end do

  end subroutine evaluate_barycentric_multi

  subroutine evaluate_cartesian_multi( &
       num_nodes, dimension_, nodes, degree, num_vals, param_vals, evaluated) &
       bind(c, name='evaluate_cartesian_multi')

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This mostly copies evaluate_barycentric_multi but does not just
    !       call it directly. This is to avoid copying param_vals.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in) :: degree, num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 2)
    real(c_double), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    integer(c_int) :: k, binom_val, index_, new_index
    real(c_double) :: row_result(num_vals, dimension_)
    real(c_double) :: lambda1_vals(num_vals)

    index_ = num_nodes
    forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
       evaluated(new_index, :) = nodes(index_, :)
    end forall

    if (degree == 0) then
       return
    end if

    lambda1_vals = 1.0_dp - param_vals(:, 1) - param_vals(:, 2)

    binom_val = 1
    do k = degree - 1, 0, -1
       ! We want to go from (d C (k + 1)) to (d C k).
       binom_val = (binom_val * (k + 1)) / (degree - k)
       index_ = index_ - 1  ! Step to last element in row.
       !     k = d - 1, d - 2, ...
       ! d - k =     1,     2, ...
       ! We know row k has (d - k + 1) elements.
       new_index = index_ - degree + k  ! First element in row.

       ! lambda1 = param_vals(:, 1)
       ! lambda2 = param_vals(:, 1)
       call evaluate_curve_barycentric( &
            degree - k, dimension_, nodes(new_index:index_, :), &
            num_vals, lambda1_vals, param_vals(:, 1), row_result)

       ! Update index for next iteration.
       index_ = new_index

       ! lambda3 = param_vals(:, 2)
       forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
          evaluated(new_index, :) = ( &
               param_vals(new_index, 2) * evaluated(new_index, :) + &
               binom_val * row_result(new_index, :))
       end forall
    end do

  end subroutine evaluate_cartesian_multi

  subroutine jacobian_both( &
       num_nodes, dimension_, nodes, degree, new_nodes) &
       bind(c, name='jacobian_both')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(out) :: &
         new_nodes(num_nodes - degree - 1, 2 * dimension_)
    ! Variables outside of signature.
    integer(c_int) :: index_, i, j, k, num_vals

    index_ = 1
    i = 1
    j = degree + 2
    do num_vals = degree, 1, -1
       do k = 0, num_vals - 1
          ! jacobian_s
          new_nodes(index_, :dimension_) = nodes(i + 1, :) - nodes(i, :)
          ! jacobian_t
          new_nodes(index_, dimension_ + 1:) = nodes(j, :) - nodes(i, :)
          ! Update the indices
          index_ = index_ + 1
          i = i + 1
          j = j + 1
       end do
       ! In between each row, the index_ gains an extra value.
       i = i + 1
    end do

    new_nodes = degree * new_nodes

  end subroutine jacobian_both

  subroutine jacobian_det( &
       num_nodes, nodes, degree, num_vals, param_vals, evaluated) &
       bind(c, name='jacobian_det')

    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: nodes(num_nodes, 2)
    integer(c_int), intent(in) :: degree, num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 2)
    real(c_double), intent(out) :: evaluated(num_vals)
    ! Variables outside of signature.
    real(c_double) :: jac_nodes(num_nodes - degree - 1, 4)
    real(c_double) :: Bs_Bt_vals(num_vals, 4)
    real(c_double) :: determinant

    call jacobian_both( &
         num_nodes, 2, nodes, degree, jac_nodes)
    if (degree == 1) then
       determinant = ( &
            jac_nodes(1, 1) * jac_nodes(1, 4) - &
            jac_nodes(1, 2) * jac_nodes(1, 3))
       evaluated = determinant
    else
       call evaluate_cartesian_multi( &
            num_nodes - degree - 1, 4, jac_nodes, degree - 1, &
            num_vals, param_vals, Bs_Bt_vals)
       evaluated = ( &
            Bs_Bt_vals(:, 1) * Bs_Bt_vals(:, 4) - &
            Bs_Bt_vals(:, 2) * Bs_Bt_vals(:, 3))
    end if

  end subroutine jacobian_det

end module surface
