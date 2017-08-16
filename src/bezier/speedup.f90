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

module speedup

  use types, only: dp
  use curve, only: &
       evaluate_curve_barycentric, evaluate_multi, evaluate_hodograph
  implicit none
  private in_interval
  public &
       de_casteljau_one_round, linearization_error, evaluate_barycentric, &
       evaluate_barycentric_multi, evaluate_cartesian_multi, cross_product, &
       segment_intersection, bbox, jacobian_both, newton_refine_intersect, &
       jacobian_det, bbox_intersect, wiggle_interval, parallel_different, &
       from_linearized, bbox_line_intersect

  integer, parameter :: BoxIntersectionType_INTERSECTION = 0
  integer, parameter :: BoxIntersectionType_TANGENT = 1
  integer, parameter :: BoxIntersectionType_DISJOINT = 2

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
    real(dp) :: param_vals(1, 3)

    param_vals(1, 1) = lambda1
    param_vals(1, 2) = lambda2
    param_vals(1, 3) = lambda3
    call evaluate_barycentric_multi( &
         num_nodes, nodes, degree, 1, param_vals, dimension_, point)

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
    integer :: k, binom_val, index_, new_index
    real(dp) :: row_result(num_vals, dimension_)

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
             nodes(new_index:index_, :), degree - k, dimension_, &
             param_vals(:, 1), param_vals(:, 2), num_vals, row_result)

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
       num_nodes, nodes, degree, num_vals, param_vals, dimension_, evaluated)

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This mostly copies evaluate_barycentric_multi but does not just
    !       call it directly. This is to avoid copying param_vals.

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
    integer :: k, binom_val, index_, new_index
    real(dp) :: row_result(num_vals, dimension_)
    real(dp) :: lambda1_vals(num_vals)

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
             nodes(new_index:index_, :), degree - k, dimension_, &
             lambda1_vals, param_vals(:, 1), num_vals, row_result)

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

  subroutine newton_refine_intersect( &
       s, nodes1, degree1, t, nodes2, degree2, new_s, new_t)

    !f2py integer intent(hide), depend(nodes1) :: degree1 = size(nodes1, 1) - 1
    !f2py integer intent(hide), depend(nodes2) :: degree2 = size(nodes2, 1) - 1
    real(dp), intent(in) :: s
    real(dp), intent(in) :: nodes1(degree1 + 1, 2)
    integer :: degree1
    real(dp), intent(in) :: t
    real(dp), intent(in) :: nodes2(degree2 + 1, 2)
    integer :: degree2
    real(dp), intent(out) :: new_s, new_t
    ! Variables outside of signature.
    real(dp) :: param(1)
    real(dp) :: func_val(1, 2)
    real(dp) :: workspace(1, 2)
    real(dp) :: jac_mat(2, 2)
    real(dp) :: determinant, delta_s, delta_t

    param = t
    call evaluate_multi( &
         nodes2, degree2, 2, param, 1, func_val)
    param = s
    call evaluate_multi( &
         nodes1, degree1, 2, param, 1, workspace)
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

  subroutine jacobian_det( &
       num_nodes, nodes, degree, num_vals, param_vals, evaluated)

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    !f2py integer intent(hide), depend(param_vals) :: num_vals &
    !f2py     = size(param_vals, 1)
    integer :: num_nodes
    real(dp), intent(in) :: nodes(num_nodes, 2)
    integer :: degree
    integer :: num_vals
    real(dp), intent(in) :: param_vals(num_vals, 2)
    real(dp), intent(out) :: evaluated(num_vals)
    ! Variables outside of signature.
    real(dp) :: jac_nodes(num_nodes - degree - 1, 4)
    real(dp) :: Bs_Bt_vals(num_vals, 4)
    real(dp) :: determinant

    call jacobian_both( &
         num_nodes, 2, nodes, degree, jac_nodes)
    if (degree == 1) then
       determinant = ( &
            jac_nodes(1, 1) * jac_nodes(1, 4) - &
            jac_nodes(1, 2) * jac_nodes(1, 3))
       evaluated = determinant
    else
       call evaluate_cartesian_multi( &
            num_nodes - degree - 1, jac_nodes, degree - 1, &
            num_vals, param_vals, 4, Bs_Bt_vals)
       evaluated = ( &
            Bs_Bt_vals(:, 1) * Bs_Bt_vals(:, 4) - &
            Bs_Bt_vals(:, 2) * Bs_Bt_vals(:, 3))
    end if

  end subroutine jacobian_det

  subroutine bbox_intersect(num_nodes1, nodes1, num_nodes2, nodes2, enum_)

    !f2py integer intent(hide), depend(nodes1) :: num_nodes1 = size(nodes1, 1)
    !f2py integer intent(hide), depend(nodes2) :: num_nodes2 = size(nodes2, 1)
    integer :: num_nodes1, num_nodes2
    real(dp), intent(in) :: nodes1(num_nodes1, 2)
    real(dp), intent(in) :: nodes2(num_nodes2, 2)
    integer, intent(out) :: enum_
    ! Variables outside of signature.
    real(dp) :: left1, right1, bottom1, top1
    real(dp) :: left2, right2, bottom2, top2

    call bbox(num_nodes1, nodes1, left1, right1, bottom1, top1)
    call bbox(num_nodes2, nodes2, left2, right2, bottom2, top2)

    if ( &
         right2 < left1 .OR. right1 < left2 .OR. &
         top2 < bottom1 .OR. top1 < bottom2) then
       enum_ = BoxIntersectionType_DISJOINT
    else if ( &
         right2 == left1 .OR. right1 == left2 .OR. &
         top2 == bottom1 .OR. top1 == bottom2) then
       enum_ = BoxIntersectionType_TANGENT
    else
       enum_ = BoxIntersectionType_INTERSECTION
    end if

  end subroutine bbox_intersect

  subroutine wiggle_interval(value_, result_, success)

    real(dp), intent(in) :: value_
    real(dp), intent(out) :: result_
    logical(1), intent(out) :: success

    success = .TRUE.
    if (-0.5_dp**45 < value_ .AND. value_ < 0.5_dp**45) then
       result_ = 0.0_dp
    else if (0.5_dp**45 <= value_ .AND. value_ <= 1.0_dp - 0.5_dp**45) then
       result_ = value_
    else if ( &
         1.0_dp - 0.5_dp**45 < value_ .AND. value_ < 1.0_dp + 0.5_dp**45) then
       result_ = 1.0_dp
    else
       success = .FALSE.
    end if

  end subroutine wiggle_interval

  subroutine parallel_different(start0, end0, start1, end1, result_)

    real(dp), intent(in) :: start0(1, 2)
    real(dp), intent(in) :: end0(1, 2)
    real(dp), intent(in) :: start1(1, 2)
    real(dp), intent(in) :: end1(1, 2)
    logical(1), intent(out) :: result_
    ! Variables outside of signature.
    real(dp) :: delta0(1, 2)
    real(dp) :: val1, val2, val3

    delta0 = end0 - start0
    call cross_product(start0, delta0, val1)  ! line0_const
    call cross_product(start1, delta0, val2)  ! start1_against

    if (val1 /= val2) then
       result_ = .TRUE.
       return
    end if

    val1 = dot_product(delta0(1, :), delta0(1, :))  ! norm0_sq
    val2 = dot_product(start1(1, :) - start0(1, :), delta0(1, :))  ! start_numer
    !      0 <= start_numer / norm0_sq <= 1
    ! <==> 0 <= start_numer            <= norm0_sq
    if (0.0_dp <= val2 .AND. val2 <= val1) then
       result_ = .FALSE.
       return
    end if

    val2 = dot_product(end1(1, :) - start0(1, :), delta0(1, :))  ! end_numer
    !      0 <= end_numer / norm0_sq <= 1
    ! <==> 0 <= end_numer            <= norm0_sq
    if (0.0_dp <= val2 .AND. val2 <= val1) then
       result_ = .FALSE.
       return
    end if

    ! We know neither the start or end parameters are in [0, 1], but
    ! they may contain [0, 1] between them.
    val3 = min(val1, val2)  ! min_val
    val1 = max(val1, val2)  ! max_val

    ! So we make sure that 0 isn't between them.
    if (val3 <= 0.0_dp .AND. 0.0_dp <= val1) then
       result_ = .FALSE.
    else
       result_ = .TRUE.
    end if

  end subroutine parallel_different

  subroutine from_linearized( &
       error1, start1, end1, start_node1, end_node1, nodes1, degree1, &
       error2, start2, end2, start_node2, end_node2, nodes2, degree2, &
       refined_s, refined_t, does_intersect, py_exc)

    !f2py integer intent(hide), depend(nodes1) :: degree1 = size(nodes1, 1) - 1
    !f2py integer intent(hide), depend(nodes2) :: degree2 = size(nodes2, 1) - 1
    real(dp), intent(in) :: error1, start1, end1
    real(dp), intent(in) :: start_node1(1, 2)
    real(dp), intent(in) :: end_node1(1, 2)
    real(dp), intent(in) :: nodes1(degree1 + 1, 2)
    integer, intent(in) :: degree1
    real(dp), intent(in) :: error2, start2, end2
    real(dp), intent(in) :: start_node2(1, 2)
    real(dp), intent(in) :: end_node2(1, 2)
    real(dp), intent(in) :: nodes2(degree2 + 1, 2)
    integer, intent(in) :: degree2
    real(dp), intent(out) :: refined_s, refined_t
    logical(1), intent(out) :: does_intersect
    integer, intent(out) :: py_exc
    ! Variables outside of signature.
    real(dp) :: s, t
    integer :: enum_
    logical(1) :: success

    py_exc = 0
    call segment_intersection( &
         start_node1, end_node1, start_node2, end_node2, s, t, success)

    if (success) then
       ! Special case for lines, allow no leeway on almost intersections.
       if (error1 == 0.0_dp .AND. (s < 0.0_dp .OR. 1.0_dp < s)) then
          does_intersect = .FALSE.
          return
       end if

       if (error2 == 0.0_dp .AND. (t < 0.0_dp .OR. 1.0_dp < t)) then
          does_intersect = .FALSE.
          return
       end if

       if (s < -(0.5_dp**16) .OR. 1.0_dp + 0.5_dp**16 < s) then
          does_intersect = .FALSE.
          return
       end if

       if (t < -(0.5_dp**16) .OR. 1.0_dp + 0.5_dp**16 < t) then
          does_intersect = .FALSE.
          return
       end if
    else
       ! Handle special case where the curves are actually lines.
       if (error1 == 0.0_dp .AND. error2 == 0.0_dp) then
          call parallel_different( &
               start_node1, end_node1, start_node2, end_node2, success)
          if (success) then
             does_intersect = .FALSE.
             return
          end if
       else
          call bbox_intersect( &
               degree1 + 1, nodes1, degree2 + 1, nodes2, enum_)
          if (enum_ == BoxIntersectionType_DISJOINT) then
             does_intersect = .FALSE.
             return
          end if
       end if

       ! Expect the wrapper code to raise
       ! NotImplementedError('Line segments parallel.')
       py_exc = 1
       return
    end if

    does_intersect = .TRUE.
    ! Now, promote ``s`` and ``t`` onto the original curves.
    s = (1.0_dp - s) * start1 + s * end1  ! orig_s
    t = (1.0_dp - t) * start2 + t * end2  ! orig_t
    ! Perform one step of Newton iteration to refine the computed
    ! values of s and t.
    call newton_refine_intersect( &
         s, nodes1, degree1, t, nodes2, degree2, refined_s, refined_t)

    call wiggle_interval(refined_s, s, success)
    if (.NOT. success) then
       ! py_exc==2 indicates ``wiggle_interval`` failed.
       py_exc = 2
       return
    end if
    refined_s = s

    call wiggle_interval(refined_t, t, success)
    if (.NOT. success) then
       ! py_exc==2 indicates ``wiggle_interval`` failed.
       py_exc = 2
       return
    end if
    refined_t = t

  end subroutine from_linearized

  pure function in_interval(value_, start, end) result(predicate)

    real(dp), intent(in) :: value_, start, end
    logical(1) :: predicate

    predicate = (start <= value_) .AND. (value_ <= end)

  end function in_interval

  subroutine bbox_line_intersect( &
       num_nodes, nodes, line_start, line_end, enum_)

    !f2py integer intent(hide), depend(nodes) :: num_nodes = size(nodes, 1)
    integer :: num_nodes
    real(dp), intent(in) :: nodes(num_nodes, 2)
    real(dp), intent(in) :: line_start(1, 2)
    real(dp), intent(in) :: line_end(1, 2)
    integer, intent(out) :: enum_
    ! Variables outside of signature.
    real(dp) :: left, right, bottom, top
    real(dp) :: segment_start(1, 2)
    real(dp) :: segment_end(1, 2)
    real(dp) :: s_curr, t_curr
    logical(1) :: success

    call bbox(num_nodes, nodes, left, right, bottom, top)

    ! Check if line start is inside the bounding box.
    if ( &
         in_interval(line_start(1, 1), left, right) .AND. &
         in_interval(line_start(1, 2), bottom, top)) then
       enum_ = BoxIntersectionType_INTERSECTION
       return
    end if

    ! Check if line end is inside the bounding box.
    if ( &
         in_interval(line_end(1, 1), left, right) .AND. &
         in_interval(line_end(1, 2), bottom, top)) then
       enum_ = BoxIntersectionType_INTERSECTION
       return
    end if

    ! NOTE: We allow ``segment_intersection`` to fail below (i.e.
    !       ``success=False``). At first, this may appear to "ignore"
    !       some potential intersections of parallel lines. However,
    !       no intersections will be missed. If parallel lines don't
    !       overlap, then there is nothing to miss. If they do overlap,
    !       then either the segment will have endpoints on the box (already
    !       covered by the checks above) or the segment will contain an
    !       entire side of the box, which will force it to intersect the 3
    !       edges that meet at the two ends of those sides. The parallel
    !       edge will be skipped, but the other two will be covered.

    ! Bottom Edge
    segment_start(1, 1) = left
    segment_start(1, 2) = bottom
    segment_end(1, 1) = right
    segment_end(1, 2) = bottom
    call segment_intersection( &
         segment_start, segment_end, line_start, line_end, &
         s_curr, t_curr, success)
    if ( &
         success .AND. &
         in_interval(s_curr, 0.0_dp, 1.0_dp) .AND. &
         in_interval(t_curr, 0.0_dp, 1.0_dp)) then
       enum_ = BoxIntersectionType_INTERSECTION
       return
    end if

    ! Right Edge
    segment_start(1, 1) = right
    segment_start(1, 2) = bottom
    segment_end(1, 1) = right
    segment_end(1, 2) = top
    call segment_intersection( &
         segment_start, segment_end, line_start, line_end, &
         s_curr, t_curr, success)
    if ( &
         success .AND. &
         in_interval(s_curr, 0.0_dp, 1.0_dp) .AND. &
         in_interval(t_curr, 0.0_dp, 1.0_dp)) then
       enum_ = BoxIntersectionType_INTERSECTION
       return
    end if

    ! Top Edge
    segment_start(1, 1) = right
    segment_start(1, 2) = top
    segment_end(1, 1) = left
    segment_end(1, 2) = top
    call segment_intersection( &
         segment_start, segment_end, line_start, line_end, &
         s_curr, t_curr, success)
    if ( &
         success .AND. &
         in_interval(s_curr, 0.0_dp, 1.0_dp) .AND. &
         in_interval(t_curr, 0.0_dp, 1.0_dp)) then
       enum_ = BoxIntersectionType_INTERSECTION
       return
    end if

    ! NOTE: We skip the "last" edge. This is because any curve
    !       that doesn't have an endpoint on a curve must cross
    !       at least two, so we will have already covered such curves
    !       in one of the branches above.

    enum_ = BoxIntersectionType_DISJOINT

  end subroutine bbox_line_intersect

end module speedup
