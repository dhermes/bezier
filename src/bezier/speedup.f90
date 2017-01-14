module speedup

  implicit none
  private
  public de_casteljau_one_round, evaluate_multi, linearization_error, &
         evaluate_barycentric, evaluate_barycentric_multi, &
         evaluate_cartesian_multi, cross_product, segment_intersection, &
         bbox, specialize_curve_generic, specialize_curve_quadratic, &
         specialize_curve, jacobian_both, evaluate_hodograph, &
         newton_refine_intersect

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
    integer :: index_
    integer :: parent_i1
    integer :: parent_i2
    integer :: parent_i3
    integer :: k, j

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
    integer :: curr_deg, index_
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
       forall (index_ = 1:curr_deg)
          weights_next(:, index_) = weights_curr(:, index_) * one_less(:)
          weights_next(:, index_ + 1) = ( &
               weights_next(:, index_ + 1) + &
               weights_curr(:, index_) * s_vals(:))
       end forall

       swap = weights_curr
       weights_curr = weights_next
       weights_next = swap
    end do

    evaluated = matmul(weights_curr, nodes)

  end subroutine evaluate_multi

  subroutine linearization_error(nodes, degree, dimension_, error)

    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    real(dp), intent(in) :: nodes(degree + 1, dimension_)
    integer :: dimension_
    integer, intent(in) :: degree
    real(dp), intent(out) :: error
    ! Variables outside of signature.
    real(dp) :: second_deriv(degree - 1, dimension_)
    real(dp) :: worst_case(dimension_)

    if (degree == 1) then
       error = 0.0_dp
       return
    end if

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

    if (degree == 0) then
       point = nodes
       return
    end if

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
    end do

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
    integer :: index_

    do index_ = 1, num_vals
       call evaluate_barycentric( &
            num_nodes, dimension_, nodes, degree, &
            param_vals(index_, 1), param_vals(index_, 2), &
            param_vals(index_, 3), evaluated(index_, :))
    end do

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
    integer :: index_

    do index_ = 1, num_vals
       call evaluate_barycentric( &
            num_nodes, dimension_, nodes, degree, &
            1.0_dp - param_vals(index_, 1) - param_vals(index_, 2), &
            param_vals(index_, 1), param_vals(index_, 2), &
            evaluated(index_, :))
    end do

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
    end if

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

  subroutine specialize_curve_generic( &
       nodes, degree, dimension_, start, end_, new_nodes)

    ! NOTE: This is a helper for ``specialize_curve`` that works on any degree.

    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    real(dp), intent(in) :: nodes(degree + 1, dimension_)
    integer :: dimension_
    integer, intent(in) :: degree
    real(dp), intent(in) :: start, end_
    real(dp), intent(out) :: new_nodes(degree + 1, dimension_)
    ! Variables outside of signature.
    real(dp) :: workspace(degree, dimension_, degree + 1)
    integer :: index_, curr_size, j
    real(dp) :: minus_start, minus_end

    minus_start = 1.0_dp - start
    minus_end = 1.0_dp - end_
    workspace(:, :, 1) = minus_start * nodes(:degree, :) + start * nodes(2:, :)
    workspace(:, :, 2) = minus_end * nodes(:degree, :) + end_ * nodes(2:, :)

    curr_size = degree
    do index_ = 3, degree + 1
       curr_size = curr_size - 1
       ! First add a new "column" (or whatever the 3rd dimension is called)
       ! at the end using ``end_``.
       workspace(:curr_size, :, index_) = ( &
            minus_end * workspace(:curr_size, :, index_ - 1) + &
            end_ * workspace(2:curr_size + 1, :, index_ - 1))
       ! Update all the values in place by using de Casteljau with the
       ! ``start`` parameter.
       forall (j = 1:index_ - 1)
          workspace(:curr_size, :, j) = ( &
               minus_start * workspace(:curr_size, :, j) + &
               start * workspace(2:curr_size + 1, :, j))
       end forall
    end do

    ! Move the final "column" (or whatever the 3rd dimension is called)
    ! of the workspace into ``new_nodes``.
    forall (index_ = 1:degree + 1)
       new_nodes(index_, :) = workspace(1, :, index_)
    end forall

  end subroutine specialize_curve_generic

  subroutine specialize_curve_quadratic( &
       nodes, dimension_, start, end_, new_nodes)

    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    real(dp), intent(in) :: nodes(3, dimension_)
    integer :: dimension_
    real(dp), intent(in) :: start, end_
    real(dp), intent(out) :: new_nodes(3, dimension_)
    ! Variables outside of signature.
    real(dp) :: minus_start, minus_end, prod_both

    minus_start = 1.0_dp - start
    minus_end = 1.0_dp - end_
    prod_both = start * end_

    new_nodes(1, :) = ( &
         minus_start * minus_start * nodes(1, :) + &
         2.0_dp * start * minus_start * nodes(2, :) + &
         start * start * nodes(3, :))
    new_nodes(2, :) = ( &
         minus_start * minus_end * nodes(1, :) + &
         (end_ + start - 2.0_dp * prod_both) * nodes(2, :) + &
         prod_both * nodes(3, :))
    new_nodes(3, :) = ( &
         minus_end * minus_end * nodes(1, :) + &
         2.0_dp * end_ * minus_end * nodes(2, :) + &
         end_ * end_ * nodes(3, :))

  end subroutine specialize_curve_quadratic

  subroutine specialize_curve( &
       nodes, degree, dimension_, start, end_, curve_start, curve_end, &
       new_nodes, true_start, true_end)

    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    real(dp), intent(in) :: nodes(degree + 1, dimension_)
    integer :: dimension_
    integer, intent(in) :: degree
    real(dp), intent(in) :: start, end_, curve_start, curve_end
    real(dp), intent(out) :: new_nodes(degree + 1, dimension_)
    real(dp), intent(out) :: true_start, true_end
    ! Variables outside of signature.
    real(dp) :: interval_delta

    if (degree == 1) then
       new_nodes(1, :) = (1.0_dp - start) * nodes(1, :) + start * nodes(2, :)
       new_nodes(2, :) = (1.0_dp - end_) * nodes(1, :) + end_ * nodes(2, :)
    else if (degree == 2) then
       call specialize_curve_quadratic( &
            nodes, dimension_, start, end_, new_nodes)
    else
       call specialize_curve_generic( &
            nodes, degree, dimension_, start, end_, new_nodes)
    end if

    ! Now, compute the new interval.
    interval_delta = curve_end - curve_start
    true_start = curve_start + start * interval_delta
    true_end = curve_start + end_ * interval_delta

  end subroutine specialize_curve

  subroutine jacobian_both( &
       num_nodes, dimension_, nodes, degree, new_nodes)

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    !f2py integer, depend(nodes) :: dimension_ = size(nodes, 2)
    integer :: num_nodes
    integer :: dimension_
    real(dp), intent(in) :: nodes(num_nodes, dimension_)
    integer :: degree
    real(dp), intent(out) :: new_nodes(num_nodes - degree - 1, 2 * dimension_)
    ! Variables outside of signature.
    integer :: index_, i, j, k, num_vals

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

  subroutine evaluate_hodograph(s, nodes, dimension_, degree, hodograph)

    !f2py integer intent(hide), depend(nodes) :: dimension_ = size(nodes, 2)
    real(dp), intent(in) :: s
    real(dp), intent(in) :: nodes(degree + 1, dimension_)
    integer :: dimension_
    integer, intent(in) :: degree
    real(dp), intent(out) :: hodograph(1, dimension_)
    ! Variables outside of signature.
    real(dp) :: first_deriv(degree, dimension_)
    real(dp) :: param(1)

    first_deriv = nodes(2:, :) - nodes(:degree, :)
    param = s
    call evaluate_multi( &
         degree, dimension_, first_deriv, 1, param, hodograph)
    hodograph = degree * hodograph

  end subroutine evaluate_hodograph

  subroutine newton_refine_intersect( &
       s, nodes1, degree1, t, nodes2, degree2, new_s, new_t)

    real(dp), intent(in) :: s
    real(dp), intent(in) :: nodes1(degree1 + 1, 2)
    integer, intent(in) :: degree1
    real(dp), intent(in) :: t
    real(dp), intent(in) :: nodes2(degree2 + 1, 2)
    integer, intent(in) :: degree2
    real(dp), intent(out) :: new_s, new_t
    ! Variables outside of signature.
    real(dp) :: param(1)
    real(dp) :: func_val(1, 2)
    real(dp) :: workspace(1, 2)
    real(dp) :: jac_mat(2, 2)
    real(dp) :: determinant, delta_s, delta_t

    param = t
    call evaluate_multi( &
         degree2 + 1, 2, nodes2, 1, param, func_val)
    param = s
    call evaluate_multi( &
         degree1 + 1, 2, nodes1, 1, param, workspace)
    func_val = func_val - workspace

    if (all(func_val == 0.0_dp)) then
       new_s = s
       new_t = t
       return
    end if

    call evaluate_hodograph(s, nodes1, 2, degree1, jac_mat(1:1, :))
    ! NOTE: We actually want the negative, since we want -B2'(t), but
    !       since we manually solve the system, it's just algebra
    !       to figure out how to use the negative values.
    call evaluate_hodograph(t, nodes2, 2, degree2, jac_mat(2:2, :))

    determinant = ( &
         jac_mat(1, 1) * jac_mat(2, 2) - jac_mat(1, 2) * jac_mat(2, 1))

    ! NOTE: We manually invert the 2x2 system ([ds, dt] J)^T = f^T.
    delta_s = ( &
         (jac_mat(2, 2) * func_val(1, 1) - &
         jac_mat(2, 1) * func_val(1, 2)) / determinant)
    new_s = s + delta_s

    delta_t = ( &
         (jac_mat(1, 2) * func_val(1, 1) - &
         jac_mat(1, 1) * func_val(1, 2)) / determinant)
    new_t = t + delta_t

  end subroutine newton_refine_intersect

end module speedup
