module speedup

  implicit none
  private
  public de_casteljau_one_round, evaluate_multi, linearization_error, &
         evaluate_barycentric, evaluate_barycentric_multi, &
         evaluate_cartesian_multi, cross_product, segment_intersection, &
         bbox

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

  subroutine evaluate_barycentric( &
       num_nodes, dimension_, nodes, degree, &
       lambda1, lambda2, lambda3, point)

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This assumes degree >= 1.

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    integer :: num_nodes
    integer :: dimension_
    real(dp), intent(in) :: nodes(num_nodes, dimension_)
    integer :: degree
    real(dp), intent(in) :: lambda1
    real(dp), intent(in) :: lambda2
    real(dp), intent(in) :: lambda3
    real(dp), intent(out) :: point(1, dimension_)
    ! Variables outside of signature.
    integer :: curr_deg, curr_num_nodes, next_num_nodes
    real(dp) :: workspace_curr(num_nodes - degree - 1, dimension_)
    real(dp) :: workspace_next(num_nodes - 2 * degree - 1, dimension_)
    real(dp) :: swap(num_nodes - 2 * degree - 1, dimension_)

    ! Do the first round of de Casteljau (this assumes degree >= 1).
    call de_casteljau_one_round( &
         num_nodes, dimension_, nodes, degree, &
         lambda1, lambda2, lambda3, workspace_curr)

    curr_num_nodes = num_nodes - degree - 1
    do curr_deg = degree - 1, 1, -1
       next_num_nodes = curr_num_nodes - curr_deg - 1
       call de_casteljau_one_round( &
            curr_num_nodes, dimension_, &
            workspace_curr(:curr_num_nodes, :), curr_deg, &
            lambda1, lambda2, lambda3, &
            workspace_next(:next_num_nodes, :))

       ! Swap the relevant data in the current and next workspaces
       curr_num_nodes = next_num_nodes
       swap(:curr_num_nodes, :) = workspace_curr(:curr_num_nodes, :)
       workspace_curr(:curr_num_nodes, :) = workspace_next(:curr_num_nodes, :)
       workspace_next(:curr_num_nodes, :) = swap(:curr_num_nodes, :)
    enddo

    point = workspace_curr(1:1, :)

  end subroutine evaluate_barycentric

  subroutine evaluate_barycentric_multi( &
       num_nodes, nodes, degree, num_vals, param_vals, dimension_, evaluated)

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This assumes degree >= 1.

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    !f2py integer depend(nodes) :: dimension_ = size(nodes, 2)
    !f2py integer intent(hide), depend(param_vals) :: num_vals &
    !f2py     = size(param_vals, 1)
    integer :: num_nodes
    integer :: dimension_
    real(dp), intent(in) :: nodes(num_nodes, dimension_)
    integer :: degree
    integer :: num_vals
    real(dp), intent(in) :: param_vals(num_vals, 3)
    real(dp), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    integer :: index

    do index = 1, num_vals
       call evaluate_barycentric( &
            num_nodes, dimension_, nodes, degree, &
            param_vals(index, 1), param_vals(index, 2), &
            param_vals(index, 3), evaluated(index, :))
    enddo

  end subroutine evaluate_barycentric_multi

  subroutine evaluate_cartesian_multi( &
       num_nodes, nodes, degree, num_vals, param_vals, dimension_, evaluated)

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This assumes degree >= 1.

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    !f2py integer depend(nodes) :: dimension_ = size(nodes, 2)
    !f2py integer intent(hide), depend(param_vals) :: num_vals &
    !f2py     = size(param_vals, 1)
    integer :: num_nodes
    integer :: dimension_
    real(dp), intent(in) :: nodes(num_nodes, dimension_)
    integer :: degree
    integer :: num_vals
    real(dp), intent(in) :: param_vals(num_vals, 2)
    real(dp), intent(out) :: evaluated(num_vals, dimension_)
    ! Variables outside of signature.
    integer :: index

    do index = 1, num_vals
       call evaluate_barycentric( &
            num_nodes, dimension_, nodes, degree, &
            1.0_dp - param_vals(index, 1) - param_vals(index, 2), &
            param_vals(index, 1), param_vals(index, 2), &
            evaluated(index, :))
    enddo

  end subroutine evaluate_cartesian_multi

  subroutine cross_product(vec0, vec1, result_)

    real(dp), intent(in) :: vec0(1, 2)
    real(dp), intent(in) :: vec1(1, 2)
    real(dp), intent(out) :: result_

    result_ = vec0(1, 1) * vec1(1, 2) - vec0(1, 2) * vec1(1, 1)

  end subroutine cross_product

  subroutine segment_intersection(start0, end0, start1, end1, s, t, success)

    real(dp), intent(in) :: start0(1, 2)
    real(dp), intent(in) :: end0(1, 2)
    real(dp), intent(in) :: start1(1, 2)
    real(dp), intent(in) :: end1(1, 2)
    real(dp), intent(out) :: s, t
    logical(1), intent(out) :: success
    ! Variables outside of signature.
    real(dp) :: delta0(1, 2)
    real(dp) :: delta1(1, 2)
    real(dp) :: start_delta(1, 2)
    real(dp) :: cross_d0_d1
    real(dp) :: other_cross

    delta0 = end0 - start0
    delta1 = end1 - start1
    call cross_product(delta0, delta1, cross_d0_d1)

    if (cross_d0_d1 == 0.0_dp) then
       success = .FALSE.
    else
       start_delta = start1 - start0
       call cross_product(start_delta, delta1, other_cross)
       s = other_cross / cross_d0_d1
       call cross_product(start_delta, delta0, other_cross)
       t = other_cross / cross_d0_d1
       success = .TRUE.
    endif

  end subroutine segment_intersection

  subroutine bbox(num_nodes, nodes, left, right, bottom, top)

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    integer :: num_nodes
    real(dp), intent(in) :: nodes(num_nodes, 2)
    real(dp), intent(out) :: left, right, bottom, top
    ! Variables outside of signature.
    real(dp) :: workspace(2)

    workspace = minval(nodes, 1)
    left = workspace(1)
    bottom = workspace(2)
    workspace = maxval(nodes, 1)
    right = workspace(1)
    top = workspace(2)

  end subroutine bbox

end module speedup
