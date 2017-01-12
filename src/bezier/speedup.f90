module speedup

  implicit none
  private
  public de_casteljau_one_round, evaluate_multi, linearization_error

  ! NOTE: This still relies on .f2py_f2cmap being present
  !       in the directory that build is called from.
  integer, parameter :: dp=kind(0.d0)

contains

  subroutine de_casteljau_one_round( &
       num_nodes, dimension_, nodes, degree, &
       lambda1, lambda2, lambda3, new_nodes)

    ! NOTE: This is de Casteljau on a Bezier surface / triangle.

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    integer :: num_nodes
    integer :: dimension_
    real(dp), intent(in) :: nodes(num_nodes, dimension_)
    integer :: degree
    real(dp), intent(in) :: lambda1
    real(dp), intent(in) :: lambda2
    real(dp), intent(in) :: lambda3
    real(dp), intent(out) :: new_nodes(num_nodes - degree - 1, dimension_)
    ! Variables outside of signature.
    integer :: index
    integer :: parent_i1
    integer :: parent_i2
    integer :: parent_i3
    integer :: k, j

    index = 1
    parent_i1 = 1
    parent_i2 = 2
    parent_i3 = degree + 2
    ! NOTE: Throughout for index <--> (i, j, k) in the (degree - 1)
    !       triangle, we have parent_i1 = index + k <--> (i + 1, j, k),
    !       parent_i2 = index + k + 1 <--> (i, j + 1, k) and
    !       parent_i3 = index + degree + 1 <--> (i, j, k + 1).
    do k = 0, degree - 1
       do j = 0, degree - k - 1
          ! NOTE: i = (degree - 1) - j - k
          new_nodes(index, :) = ( &
               lambda1 * nodes(parent_i1, :) + &
               lambda2 * nodes(parent_i2, :) + &
               lambda3 * nodes(parent_i3, :))
          ! Update all the indices.
          parent_i1 = parent_i1 + 1
          parent_i2 = parent_i2 + 1
          parent_i3 = parent_i3 + 1
          index = index + 1
       enddo
       ! Update the indices that depend on k.
       parent_i1 = parent_i1 + 1
       parent_i2 = parent_i2 + 1
    enddo

  end subroutine de_casteljau_one_round

  subroutine evaluate_multi( &
       num_nodes, dimension_, nodes, num_vals, s_vals, evaluated)

    ! NOTE: This is evaluate_multi for a Bezier curve.

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    !f2py integer intent(hide), depend(s_vals) :: num_vals = size(s_vals)
    integer :: num_nodes
    integer :: dimension_
    real(dp), intent(in) :: nodes(num_nodes, dimension_)
    integer :: num_vals
    real(dp), intent(in) :: s_vals(num_vals)
    real(dp), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    integer :: curr_deg, index
    real(dp) :: one_less(num_vals)
    real(dp) :: weights_next(num_vals, num_nodes)
    real(dp) :: weights_curr(num_vals, num_nodes)
    real(dp) :: swap(num_vals, num_nodes)

    ! Degree 0 only occupies first column
    weights_curr(:, 1) = 1.0_dp
    weights_curr(:, 2:) = 0.0_dp
    ! Zero out "next" workspace
    weights_next = 0.0_dp

    one_less(:) = 1.0_dp - s_vals(:)

    ! Increase from degree 0 to (degree - 1).
    do curr_deg = 1, num_nodes - 1
       forall (index = 1:curr_deg)
          weights_next(:, index) = weights_curr(:, index) * one_less(:)
          weights_next(:, index + 1) = ( &
               weights_next(:, index + 1) + &
               weights_curr(:, index) * s_vals(:))
       end forall

       swap = weights_curr
       weights_curr = weights_next
       weights_next = swap
    enddo

    evaluated = matmul(weights_curr, nodes)

  end subroutine evaluate_multi

  subroutine linearization_error(nodes, degree, dimension_, error)

    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    real(dp), intent(in) :: nodes(degree + 1, dimension_)
    integer :: dimension_
    integer :: degree
    real(dp), intent(out) :: error
    ! Variables outside of signature.
    real(dp) :: second_deriv(degree - 1, dimension_)
    real(dp) :: worst_case(dimension_)

    if (degree == 1) then
       error = 0.0_dp
       return
    endif

    second_deriv = ( &
         nodes(:degree - 1, :) - &
         2.0_dp * nodes(2:degree, :) + &
         nodes(3:, :))
    worst_case = maxval(abs(second_deriv), 1)
    error = 0.125_dp * degree * (degree - 1) * norm2(worst_case)
 end subroutine linearization_error

end module speedup
