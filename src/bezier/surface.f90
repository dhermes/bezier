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

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_bool
  use types, only: dp
  use curve, only: evaluate_curve_barycentric
  implicit none
  private specialize_workspace_sizes, specialize_surface_one_round
  public &
       de_casteljau_one_round, evaluate_barycentric, &
       evaluate_barycentric_multi, evaluate_cartesian_multi, jacobian_both, &
       jacobian_det, specialize_surface, subdivide_nodes, compute_edge_nodes

contains

  subroutine de_casteljau_one_round( &
       num_nodes, dimension_, nodes, degree, &
       lambda1, lambda2, lambda3, new_nodes) &
       bind(c, name='de_casteljau_one_round')

    ! NOTE: This is de Casteljau on a Bezier surface / triangle.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(in) :: lambda1
    real(c_double), intent(in) :: lambda2
    real(c_double), intent(in) :: lambda3
    real(c_double), intent(out) :: &
         new_nodes(dimension_, num_nodes - degree - 1)
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
          new_nodes(:, index_) = ( &
               lambda1 * nodes(:, parent_i1) + &
               lambda2 * nodes(:, parent_i2) + &
               lambda3 * nodes(:, parent_i3))
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
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(in) :: lambda1
    real(c_double), intent(in) :: lambda2
    real(c_double), intent(in) :: lambda3
    real(c_double), intent(out) :: point(dimension_, 1)
    ! Variables outside of signature.
    real(c_double) :: param_vals(3, 1)

    param_vals(1, 1) = lambda1
    param_vals(2, 1) = lambda2
    param_vals(3, 1) = lambda3
    call evaluate_barycentric_multi( &
         num_nodes, dimension_, nodes, degree, 1, param_vals, point)

  end subroutine evaluate_barycentric

  subroutine evaluate_barycentric_multi( &
       num_nodes, dimension_, nodes, degree, num_vals, param_vals, evaluated) &
       bind(c, name='evaluate_barycentric_multi')

    ! NOTE: This evaluation is on a Bezier surface / triangle.
    ! NOTE: This assumes degree >= 1.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: degree, num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 3)
    real(c_double), intent(out) :: evaluated(dimension_, num_vals)
    ! Variables outside of signature.
    integer(c_int) :: k, binom_val, index_, new_index
    real(c_double) :: row_result(dimension_, num_vals)

    index_ = num_nodes
    forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
       evaluated(:, new_index) = nodes(:, index_)
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
            degree - k + 1, dimension_, nodes(:, new_index:index_), &
            num_vals, param_vals(:, 1), param_vals(:, 2), row_result)

       ! Update index for next iteration.
       index_ = new_index

       ! lambda3 = param_vals(:, 3)
       forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
          evaluated(:, new_index) = ( &
               param_vals(new_index, 3) * evaluated(:, new_index) + &
               binom_val * row_result(:, new_index))
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
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: degree, num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 2)
    real(c_double), intent(out) :: evaluated(dimension_, num_vals)
    ! Variables outside of signature.
    integer(c_int) :: k, binom_val, index_, new_index
    real(c_double) :: row_result(dimension_, num_vals)
    real(c_double) :: lambda1_vals(num_vals)

    index_ = num_nodes
    forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
       evaluated(:, new_index) = nodes(:, index_)
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
            degree - k + 1, dimension_, nodes(:, new_index:index_), &
            num_vals, lambda1_vals, param_vals(:, 1), row_result)

       ! Update index for next iteration.
       index_ = new_index

       ! lambda3 = param_vals(:, 2)
       forall (new_index = 1:num_vals)  ! Borrow new_index for this loop.
          evaluated(:, new_index) = ( &
               param_vals(new_index, 2) * evaluated(:, new_index) + &
               binom_val * row_result(:, new_index))
       end forall
    end do

  end subroutine evaluate_cartesian_multi

  subroutine jacobian_both( &
       num_nodes, dimension_, nodes, degree, new_nodes) &
       bind(c, name='jacobian_both')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(out) :: &
         new_nodes(2 * dimension_, num_nodes - degree - 1)
    ! Variables outside of signature.
    integer(c_int) :: index_, i, j, k, num_vals

    index_ = 1
    i = 1
    j = degree + 2
    do num_vals = degree, 1, -1
       do k = 0, num_vals - 1
          ! jacobian_s
          new_nodes(:dimension_, index_) = nodes(:, i + 1) - nodes(:, i)
          ! jacobian_t
          new_nodes(dimension_ + 1:, index_) = nodes(:, j) - nodes(:, i)
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
    real(c_double), intent(in) :: nodes(2, num_nodes)
    integer(c_int), intent(in) :: degree, num_vals
    real(c_double), intent(in) :: param_vals(num_vals, 2)
    real(c_double), intent(out) :: evaluated(num_vals)
    ! Variables outside of signature.
    real(c_double) :: jac_nodes(4, num_nodes - degree - 1)
    real(c_double) :: Bs_Bt_vals(4, num_vals)
    real(c_double) :: determinant

    call jacobian_both( &
         num_nodes, 2, nodes, degree, jac_nodes)
    if (degree == 1) then
       determinant = ( &
            jac_nodes(1, 1) * jac_nodes(4, 1) - &
            jac_nodes(2, 1) * jac_nodes(3, 1))
       evaluated = determinant
    else
       call evaluate_cartesian_multi( &
            num_nodes - degree - 1, 4, jac_nodes, degree - 1, &
            num_vals, param_vals, Bs_Bt_vals)
       evaluated = ( &
            Bs_Bt_vals(1, :) * Bs_Bt_vals(4, :) - &
            Bs_Bt_vals(2, :) * Bs_Bt_vals(3, :))
    end if

  end subroutine jacobian_det

  subroutine specialize_workspace_sizes(degree, size_odd, size_even)
    integer(c_int), intent(in) :: degree
    integer(c_int), intent(out) :: size_odd, size_even

    ! NOTE: This is a helper for ``specialize_surface``.
    if (mod(degree, 2) == 1) then
       size_odd = ((degree + 1) * (degree + 3)**2 * (degree + 5)) / 64
       size_even = size_odd
    else if (mod(degree, 4) == 0) then
       size_odd = (degree * (degree + 2) * (degree + 4) * (degree + 6)) / 64
       size_even = (((degree + 2) * (degree + 4)) / 8)**2
    else
       size_odd = (((degree + 2) * (degree + 4)) / 8)**2
       size_even = (degree * (degree + 2) * (degree + 4) * (degree + 6)) / 64
    end if

  end subroutine specialize_workspace_sizes

  subroutine specialize_surface_one_round( &
       dimension_, num_read, read_nodes, num_write, write_nodes, &
       size_read, size_write, step, local_degree, &
       weights_a, weights_b, weights_c)

    ! INVARIANT: `step + local_degree == (total) degree + 1`
    ! INVARIANT: `size_write == size_read - local_degree - 1`

    integer(c_int), intent(in) :: dimension_
    integer(c_int), intent(in) :: num_read
    real(c_double), intent(in) :: read_nodes(dimension_, num_read)
    integer(c_int), intent(in) :: num_write
    real(c_double), intent(inout) :: write_nodes(dimension_, num_write)
    integer(c_int), intent(in) :: size_read, size_write
    integer(c_int), intent(in) :: step, local_degree
    real(c_double), intent(in) :: weights_a(3), weights_b(3), weights_c(3)
    ! Variables outside of signature.
    integer(c_int) :: write_index, new_write, read_index, new_read, i

    ! First: (step, 0, 0)
    call de_casteljau_one_round( &
         size_read, dimension_, read_nodes(:, 1:size_read), &
         local_degree, weights_a(1), weights_a(2), weights_a(3), &
         write_nodes(:, 1:size_write))

    ! Second: (i, j, 0) for j > 0, i + j = step
    write_index = size_write + 1
    read_index = 1
    do i = 1, step
       new_write = write_index + size_write - 1
       new_read = read_index + size_read - 1
       call de_casteljau_one_round( &
            size_read, dimension_, read_nodes(:, read_index:new_read), &
            local_degree, weights_b(1), weights_b(2), weights_b(3), &
            write_nodes(:, write_index:new_write))
       ! Update the indices.
       write_index = new_write + 1
       read_index = new_read + 1
    end do

    ! Third:  (i, j, k) for k > 0, i + j + k = step
    read_index = 1
    do i = 1, ((step + 1) * step) / 2
       new_write = write_index + size_write - 1
       new_read = read_index + size_read - 1
       call de_casteljau_one_round( &
            size_read, dimension_, read_nodes(:, read_index:new_read), &
            local_degree, weights_c(1), weights_c(2), weights_c(3), &
            write_nodes(:, write_index:new_write))
       ! Update the indices.
       write_index = new_write + 1
       read_index = new_read + 1
    end do

  end subroutine specialize_surface_one_round

  subroutine specialize_surface( &
       num_nodes, dimension_, nodes, degree, &
       weights_a, weights_b, weights_c, specialized) &
       bind(c, name='specialize_surface')

    ! This re-parameterizes by "specializing" to a subregion of the unit
    ! triangle. This happens in waves by applying one round of de Casteljau
    ! with the given sets of barycentric weights. For example, in the degree
    ! 1 case, we transform
    !     [n0]          [dC(n0(1:3), a)]
    !     [  ] --> n1 = [dC(n0(1:3), b)]
    !     [  ]          [dC(n0(1:3), c)]
    ! where `n` on the LHS is 3xDIM and each of the transformed sets of
    ! nodes are 1xDIM (since de Casteljau drops the degree by 1). (NOTE:
    ! here I'm using `dC` as a shorthand for "de Casteljau" and `a` as a
    ! shorthand `weights_a`.) In the degree 2 case:
    !     [n0]          [dC(n0(1:6), a)]          [dC(n1(1:3), a)]
    !     [  ]          [              ]          [dC(n1(1:3), b)]
    !     [  ]          [              ] --> n2 = [dC(n1(4:6), b)]
    !     [  ] --> n1 = [dC(n0(1:6), b)]          [dC(n1(1:3), c)]
    !     [  ]          [              ]          [dC(n1(4:6), c)]
    !     [  ]          [              ]          [dC(n1(7:9), c)]
    !                   [dC(n0(1:6), c)]
    !                   [              ]
    !                   [              ]
    ! In the degree 3 case:
    !     [n0]          [dC(n0(1:10), a)]          [dC(n1( 1:6 ), a)]
    !     [  ]          [               ]          [                ]
    !     [  ]          [               ]          [                ]
    !     [  ]          [               ]          [dC(n1( 1:6 ), b)]
    !     [  ]          [               ]          [                ]
    !     [  ] --> n1 = [               ]          [                ]
    !     [  ]          [dC(n0(1:10), b)]          [dC(n1( 7:12), b)]
    !     [  ]          [               ]          [                ]
    !     [  ]          [               ]          [                ]
    !     [  ]          [               ] --> n2 = [dC(n1( 1:6 ), c)]
    !                   [               ]          [                ]
    !                   [               ]          [                ]
    !                   [dC(n0(1:10), c)]          [dC(n1( 7:12), c)]
    !                   [               ]          [                ]
    !                   [               ]          [                ]
    !                   [               ]          [dC(n1(13:18), c)]
    !                   [               ]          [                ]
    !                   [               ]          [                ]
    ! ... CONTINUING TO DEGREE 3 ...
    !                 [dC(n2( 1:3 ), a)]
    !                 [dC(n2( 1:3 ), b)]
    !                 [dC(n2( 4:6 ), b)]
    !                 [dC(n2( 7:9 ), b)]
    !                 [dC(n2( 1:3 ), c)]
    !     n2 --> n3 = [dC(n2( 4:6 ), c)]
    !                 [dC(n2( 7:9 ), c)]
    !                 [dC(n2(10:12), c)]
    !                 [dC(n2(13:15), c)]
    !                 [dC(n2(16:18), c)]
    ! For **all** cases, we are tracking applications of weights at
    ! each step:
    !     STEP 1: {a}, {b}, {c}
    !     STEP 2: {aa}, {ab}, {bb}, {ac}, {bc}, {cc}
    !     STEP 3: {aaa}, {aab}, {abb}, {bbb}, {aac}, ...
    ! and each of these quantities can be computed using the ones computed
    ! at the previous step. So in STEP 2, we compute {aa} by applying
    ! the `a` weights to {a}, we compute {ab} and {bb} by applying the
    ! `b` weights to {a} and {b} and we compute {ac}, {bc}, {cc} by
    ! applying the `c` weights to the whole set.
    !
    ! For a "general" step, we can treat the triple (#a, #b, #c) just like
    ! we treat an index triple (i, j, k), which we know has linear index
    ! (1-based) `1 + j + (k/2) (2(i + j) + k + 3)`. This is very useful because
    ! we can find the **previous** quantity by using this index. For example,
    ! to compute {aabccc} by applying `c` to {aabcc}, we know the index
    ! of the {aabcc} is 13 and the index of {aabccc} is 20. In fact, one can
    ! show that when dropping from (i, j, k) to (i, j, k - 1) the index will
    ! always drop by `(i + j + k) + 1 = step + 1`. Unfortunately, some values
    ! have no `c` weights at all, i.e. `k = 0`. For these special cases we
    ! also use the formula to show `(step, 0, 0) --> (step - 1, 0, 0)` has
    ! no change in index (both are the first index). From `(i, j, 0)` to
    ! `(i, j - 1, 0)` drops the index by `1`.
    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(in) :: weights_a(3), weights_b(3), weights_c(3)
    real(c_double), intent(out) :: specialized(dimension_, num_nodes)
    ! Variables outside of signature.
    integer(c_int) :: num_curves, size, size_odd, size_even
    real(c_double), allocatable :: workspace_odd(:, :), workspace_even(:, :)
    integer(c_int) :: delta_nc, delta_size, size_new, step
    logical(c_bool) :: is_even

    num_curves = 1
    size = ((degree + 1) * (degree + 2)) / 2
    call specialize_workspace_sizes(degree, size_odd, size_even)
    allocate(workspace_odd(dimension_, size_odd))
    allocate(workspace_even(dimension_, size_even))

    ! `num_curves` and `delta_size` are triangular numbers going in
    ! the opposite direction. `num_curves` will go from 1, 3, 6, 10, ...
    ! and `delta_size` will go (d + 2 C 2), (d + 1 C 2), ...
    ! Since the delta between consecutive terms grows linearly, we track
    ! them as well to avoid multiplications.
    delta_nc = 2
    delta_size = -degree - 1
    is_even = .FALSE.  ! At `step == 1`.

    do step = 1, degree
       num_curves = num_curves + delta_nc
       size_new = size + delta_size
       if (step == 1) then
          ! Read from `nodes`, write to `workspace_odd`.
          call specialize_surface_one_round( &
               dimension_, num_nodes, nodes, size_odd, workspace_odd, &
               size, size_new, step, degree + 1 - step, &
               weights_a, weights_b, weights_c)
       else if (is_even) then
          ! Read from `workspace_odd`, write to `workspace_even`.
          call specialize_surface_one_round( &
               dimension_, size_odd, workspace_odd, &
               size_even, workspace_even, &
               size, size_new, step, degree + 1 - step, &
               weights_a, weights_b, weights_c)
       else
          ! Read from `workspace_even`, write to `workspace_odd`.
          call specialize_surface_one_round( &
               dimension_, size_even, workspace_even, &
               size_odd, workspace_odd, &
               size, size_new, step, degree + 1 - step, &
               weights_a, weights_b, weights_c)
       end if

       ! Update our index values.
       size = size_new
       delta_nc = delta_nc + 1
       delta_size = delta_size + 1
       is_even = .NOT. is_even  ! Switch parity.
    end do

    ! Read from the last workspace that was written to.
    if (is_even) then
       specialized = workspace_odd(:, 1:num_nodes)
    else
       specialized = workspace_even(:, 1:num_nodes)
    end if

  end subroutine specialize_surface

  subroutine subdivide_nodes( &
       num_nodes, dimension_, nodes, degree, &
       nodes_a, nodes_b, nodes_c, nodes_d) &
       bind(c, name='subdivide_nodes_surface')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(out) :: nodes_a(dimension_, num_nodes)
    real(c_double), intent(out) :: nodes_b(dimension_, num_nodes)
    real(c_double), intent(out) :: nodes_c(dimension_, num_nodes)
    real(c_double), intent(out) :: nodes_d(dimension_, num_nodes)

    if (degree == 1) then
       nodes_a(:, 1) = nodes(:, 1)
       nodes_a(:, 2) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
       nodes_a(:, 3) = 0.5_dp * (nodes(:, 1) + nodes(:, 3))
       nodes_b(:, 1) = 0.5_dp * (nodes(:, 2) + nodes(:, 3))
       nodes_b(:, 2) = nodes_a(:, 3)
       nodes_b(:, 3) = nodes_a(:, 2)
       nodes_c(:, 1) = nodes_a(:, 2)
       nodes_c(:, 2) = nodes(:, 2)
       nodes_c(:, 3) = nodes_b(:, 1)
       nodes_d(:, 1) = nodes_a(:, 3)
       nodes_d(:, 2) = nodes_b(:, 1)
       nodes_d(:, 3) = nodes(:, 3)
    else if (degree == 2) then
       nodes_a(:, 1) = nodes(:, 1)
       nodes_a(:, 2) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
       nodes_a(:, 3) = 0.25_dp * (nodes(:, 1) + 2 * nodes(:, 2) + nodes(:, 3))
       nodes_a(:, 4) = 0.5_dp * (nodes(:, 1) + nodes(:, 4))
       nodes_a(:, 5) = 0.25_dp * ( &
            nodes(:, 1) + nodes(:, 2) + nodes(:, 4) + nodes(:, 5))
       nodes_a(:, 6) = 0.25_dp * (nodes(:, 1) + 2 * nodes(:, 4) + nodes(:, 6))
       nodes_b(:, 1) = 0.25_dp * (nodes(:, 3) + 2 * nodes(:, 5) + nodes(:, 6))
       nodes_b(:, 2) = 0.25_dp * ( &
            nodes(:, 2) + nodes(:, 4) + nodes(:, 5) + nodes(:, 6))
       nodes_b(:, 3) = nodes_a(:, 6)
       nodes_b(:, 4) = 0.25_dp * ( &
            nodes(:, 2) + nodes(:, 3) + nodes(:, 4) + nodes(:, 5))
       nodes_b(:, 5) = nodes_a(:, 5)
       nodes_b(:, 6) = nodes_a(:, 3)
       nodes_c(:, 1) = nodes_a(:, 3)
       nodes_c(:, 2) = 0.5_dp * (nodes(:, 2) + nodes(:, 3))
       nodes_c(:, 3) = nodes(:, 3)
       nodes_c(:, 4) = nodes_b(:, 4)
       nodes_c(:, 5) = 0.5_dp * (nodes(:, 3) + nodes(:, 5))
       nodes_c(:, 6) = nodes_b(:, 1)
       nodes_d(:, 1) = nodes_a(:, 6)
       nodes_d(:, 2) = nodes_b(:, 2)
       nodes_d(:, 3) = nodes_b(:, 1)
       nodes_d(:, 4) = 0.5_dp * (nodes(:, 4) + nodes(:, 6))
       nodes_d(:, 5) = 0.5_dp * (nodes(:, 5) + nodes(:, 6))
       nodes_d(:, 6) = nodes(:, 6)
    else if (degree == 3) then
       nodes_a(:, 1) = nodes(:, 1)
       nodes_a(:, 2) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
       nodes_a(:, 3) = 0.25_dp * (nodes(:, 1) + 2 * nodes(:, 2) + nodes(:, 3))
       nodes_a(:, 4) = 0.125_dp * ( &
            nodes(:, 1) + 3 * nodes(:, 2) + 3 * nodes(:, 3) + nodes(:, 4))
       nodes_a(:, 5) = 0.5_dp * (nodes(:, 1) + nodes(:, 5))
       nodes_a(:, 6) = 0.25_dp * ( &
            nodes(:, 1) + nodes(:, 2) + nodes(:, 5) + nodes(:, 6))
       nodes_a(:, 7) = 0.125_dp * ( &
            nodes(:, 1) + 2 * nodes(:, 2) + nodes(:, 3) + nodes(:, 5) + &
            2 * nodes(:, 6) + nodes(:, 7))
       nodes_a(:, 8) = 0.25_dp * (nodes(:, 1) + 2 * nodes(:, 5) + nodes(:, 8))
       nodes_a(:, 9) = 0.125_dp * ( &
            nodes(:, 1) + nodes(:, 2) + 2 * nodes(:, 5) + 2 * nodes(:, 6) + &
            nodes(:, 8) + nodes(:, 9))
       nodes_a(:, 10) = 0.125_dp * ( &
            nodes(:, 1) + 3 * nodes(:, 5) + 3 * nodes(:, 8) + nodes(:, 10))
       nodes_b(:, 1) = 0.125_dp * ( &
            nodes(:, 4) + 3 * nodes(:, 7) + 3 * nodes(:, 9) + nodes(:, 10))
       nodes_b(:, 2) = 0.125_dp * ( &
            nodes(:, 3) + 2 * nodes(:, 6) + nodes(:, 7) + nodes(:, 8) + &
            2 * nodes(:, 9) + nodes(:, 10))
       nodes_b(:, 3) = 0.125_dp * ( &
            nodes(:, 2) + nodes(:, 5) + 2 * nodes(:, 6) + 2 * nodes(:, 8) + &
            nodes(:, 9) + nodes(:, 10))
       nodes_b(:, 4) = nodes_a(:, 10)
       nodes_b(:, 5) = 0.125_dp * ( &
            nodes(:, 3) + nodes(:, 4) + 2 * nodes(:, 6) + 2 * nodes(:, 7) + &
            nodes(:, 8) + nodes(:, 9))
       nodes_b(:, 6) = 0.125_dp * ( &
            nodes(:, 2) + nodes(:, 3) + nodes(:, 5) + 2 * nodes(:, 6) + &
            nodes(:, 7) + nodes(:, 8) + nodes(:, 9))
       nodes_b(:, 7) = nodes_a(:, 9)
       nodes_b(:, 8) = 0.125_dp * ( &
            nodes(:, 2) + 2 * nodes(:, 3) + nodes(:, 4) + nodes(:, 5) + &
            2 * nodes(:, 6) + nodes(:, 7))
       nodes_b(:, 9) = nodes_a(:, 7)
       nodes_b(:, 10) = nodes_a(:, 4)
       nodes_c(:, 1) = nodes_a(:, 4)
       nodes_c(:, 2) = 0.25_dp * (nodes(:, 2) + 2 * nodes(:, 3) + nodes(:, 4))
       nodes_c(:, 3) = 0.5_dp * (nodes(:, 3) + nodes(:, 4))
       nodes_c(:, 4) = nodes(:, 4)
       nodes_c(:, 5) = nodes_b(:, 8)
       nodes_c(:, 6) = 0.25_dp * ( &
            nodes(:, 3) + nodes(:, 4) + nodes(:, 6) + nodes(:, 7))
       nodes_c(:, 7) = 0.5_dp * (nodes(:, 4) + nodes(:, 7))
       nodes_c(:, 8) = nodes_b(:, 5)
       nodes_c(:, 9) = 0.25_dp * (nodes(:, 4) + 2 * nodes(:, 7) + nodes(:, 9))
       nodes_c(:, 10) = nodes_b(:, 1)
       nodes_d(:, 1) = nodes_a(:, 10)
       nodes_d(:, 2) = nodes_b(:, 3)
       nodes_d(:, 3) = nodes_b(:, 2)
       nodes_d(:, 4) = nodes_b(:, 1)
       nodes_d(:, 5) = 0.25_dp * (nodes(:, 5) + 2 * nodes(:, 8) + nodes(:, 10))
       nodes_d(:, 6) = 0.25_dp * ( &
            nodes(:, 6) + nodes(:, 8) + nodes(:, 9) + nodes(:, 10))
       nodes_d(:, 7) = 0.25_dp * (nodes(:, 7) + 2 * nodes(:, 9) + nodes(:, 10))
       nodes_d(:, 8) = 0.5_dp * (nodes(:, 8) + nodes(:, 10))
       nodes_d(:, 9) = 0.5_dp * (nodes(:, 9) + nodes(:, 10))
       nodes_d(:, 10) = nodes(:, 10)
    else if (degree == 4) then
       nodes_a(:, 1) = nodes(:, 1)
       nodes_a(:, 2) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
       nodes_a(:, 3) = 0.25_dp * (nodes(:, 1) + 2 * nodes(:, 2) + nodes(:, 3))
       nodes_a(:, 4) = 0.125_dp * ( &
            nodes(:, 1) + 3 * nodes(:, 2) + 3 * nodes(:, 3) + nodes(:, 4))
       nodes_a(:, 5) = 0.0625_dp * ( &
            nodes(:, 1) + 4 * nodes(:, 2) + 6 * nodes(:, 3) + &
            4 * nodes(:, 4) + nodes(:, 5))
       nodes_a(:, 6) = 0.5_dp * (nodes(:, 1) + nodes(:, 6))
       nodes_a(:, 7) = 0.25_dp * ( &
            nodes(:, 1) + nodes(:, 2) + nodes(:, 6) + nodes(:, 7))
       nodes_a(:, 8) = 0.125_dp * ( &
            nodes(:, 1) + 2 * nodes(:, 2) + nodes(:, 3) + nodes(:, 6) + &
            2 * nodes(:, 7) + nodes(:, 8))
       nodes_a(:, 9) = 0.0625_dp * ( &
            nodes(:, 1) + 3 * nodes(:, 2) + 3 * nodes(:, 3) + nodes(:, 4) + &
            nodes(:, 6) + 3 * nodes(:, 7) + 3 * nodes(:, 8) + nodes(:, 9))
       nodes_a(:, 10) = 0.25_dp * ( &
            nodes(:, 1) + 2 * nodes(:, 6) + nodes(:, 10))
       nodes_a(:, 11) = 0.125_dp * ( &
            nodes(:, 1) + nodes(:, 2) + 2 * nodes(:, 6) + 2 * nodes(:, 7) + &
            nodes(:, 10) + nodes(:, 11))
       nodes_a(:, 12) = 0.0625_dp * ( &
            nodes(:, 1) + 2 * nodes(:, 2) + nodes(:, 3) + 2 * nodes(:, 6) + &
            4 * nodes(:, 7) + 2 * nodes(:, 8) + nodes(:, 10) + &
            2 * nodes(:, 11) + nodes(:, 12))
       nodes_a(:, 13) = 0.125_dp * ( &
            nodes(:, 1) + 3 * nodes(:, 6) + 3 * nodes(:, 10) + nodes(:, 13))
       nodes_a(:, 14) = 0.0625_dp * ( &
            nodes(:, 1) + nodes(:, 2) + 3 * nodes(:, 6) + 3 * nodes(:, 7) + &
            3 * nodes(:, 10) + 3 * nodes(:, 11) + nodes(:, 13) + nodes(:, 14))
       nodes_a(:, 15) = 0.0625_dp * ( &
            nodes(:, 1) + 4 * nodes(:, 6) + 6 * nodes(:, 10) + &
            4 * nodes(:, 13) + nodes(:, 15))
       nodes_b(:, 1) = 0.0625_dp * ( &
            nodes(:, 5) + 4 * nodes(:, 9) + 6 * nodes(:, 12) + &
            4 * nodes(:, 14) + nodes(:, 15))
       nodes_b(:, 2) = 0.0625_dp * ( &
            nodes(:, 4) + 3 * nodes(:, 8) + nodes(:, 9) + 3 * nodes(:, 11) + &
            3 * nodes(:, 12) + nodes(:, 13) + 3 * nodes(:, 14) + nodes(:, 15))
       nodes_b(:, 3) = 0.0625_dp * ( &
            nodes(:, 3) + 2 * nodes(:, 7) + 2 * nodes(:, 8) + nodes(:, 10) + &
            4 * nodes(:, 11) + nodes(:, 12) + 2 * nodes(:, 13) + &
            2 * nodes(:, 14) + nodes(:, 15))
       nodes_b(:, 4) = 0.0625_dp * ( &
            nodes(:, 2) + nodes(:, 6) + 3 * nodes(:, 7) + 3 * nodes(:, 10) + &
            3 * nodes(:, 11) + 3 * nodes(:, 13) + nodes(:, 14) + nodes(:, 15))
       nodes_b(:, 5) = nodes_a(:, 15)
       nodes_b(:, 6) = 0.0625_dp * ( &
            nodes(:, 4) + nodes(:, 5) + 3 * nodes(:, 8) + 3 * nodes(:, 9) + &
            3 * nodes(:, 11) + 3 * nodes(:, 12) + nodes(:, 13) + nodes(:, 14))
       nodes_b(:, 7) = 0.0625_dp * ( &
            nodes(:, 3) + nodes(:, 4) + 2 * nodes(:, 7) + 3 * nodes(:, 8) + &
            nodes(:, 9) + nodes(:, 10) + 3 * nodes(:, 11) + &
            2 * nodes(:, 12) + nodes(:, 13) + nodes(:, 14))
       nodes_b(:, 8) = 0.0625_dp * ( &
            nodes(:, 2) + nodes(:, 3) + nodes(:, 6) + 3 * nodes(:, 7) + &
            2 * nodes(:, 8) + 2 * nodes(:, 10) + 3 * nodes(:, 11) + &
            nodes(:, 12) + nodes(:, 13) + nodes(:, 14))
       nodes_b(:, 9) = nodes_a(:, 14)
       nodes_b(:, 10) = 0.0625_dp * ( &
            nodes(:, 3) + 2 * nodes(:, 4) + nodes(:, 5) + 2 * nodes(:, 7) + &
            4 * nodes(:, 8) + 2 * nodes(:, 9) + nodes(:, 10) + &
            2 * nodes(:, 11) + nodes(:, 12))
       nodes_b(:, 11) = 0.0625_dp * ( &
            nodes(:, 2) + 2 * nodes(:, 3) + nodes(:, 4) + nodes(:, 6) + &
            3 * nodes(:, 7) + 3 * nodes(:, 8) + nodes(:, 9) + nodes(:, 10) + &
            2 * nodes(:, 11) + nodes(:, 12))
       nodes_b(:, 12) = nodes_a(:, 12)
       nodes_b(:, 13) = 0.0625_dp * ( &
            nodes(:, 2) + 3 * nodes(:, 3) + 3 * nodes(:, 4) + nodes(:, 5) + &
            nodes(:, 6) + 3 * nodes(:, 7) + 3 * nodes(:, 8) + nodes(:, 9))
       nodes_b(:, 14) = nodes_a(:, 9)
       nodes_b(:, 15) = nodes_a(:, 5)
       nodes_c(:, 1) = nodes_a(:, 5)
       nodes_c(:, 2) = 0.125_dp * ( &
            nodes(:, 2) + 3 * nodes(:, 3) + 3 * nodes(:, 4) + nodes(:, 5))
       nodes_c(:, 3) = 0.25_dp * ( &
            nodes(:, 3) + 2 * nodes(:, 4) + nodes(:, 5))
       nodes_c(:, 4) = 0.5_dp * ( &
            nodes(:, 4) + nodes(:, 5))
       nodes_c(:, 5) = nodes(:, 5)
       nodes_c(:, 6) = nodes_b(:, 13)
       nodes_c(:, 7) = 0.125_dp * ( &
            nodes(:, 3) + 2 * nodes(:, 4) + nodes(:, 5) + nodes(:, 7) + &
            2 * nodes(:, 8) + nodes(:, 9))
       nodes_c(:, 8) = 0.25_dp * ( &
            nodes(:, 4) + nodes(:, 5) + nodes(:, 8) + nodes(:, 9))
       nodes_c(:, 9) = 0.5_dp * (nodes(:, 5) + nodes(:, 9))
       nodes_c(:, 10) = nodes_b(:, 10)
       nodes_c(:, 11) = 0.125_dp * ( &
            nodes(:, 4) + nodes(:, 5) + 2 * nodes(:, 8) + 2 * nodes(:, 9) + &
            nodes(:, 11) + nodes(:, 12))
       nodes_c(:, 12) = 0.25_dp * ( &
            nodes(:, 5) + 2 * nodes(:, 9) + nodes(:, 12))
       nodes_c(:, 13) = nodes_b(:, 6)
       nodes_c(:, 14) = 0.125_dp * ( &
            nodes(:, 5) + 3 * nodes(:, 9) + 3 * nodes(:, 12) + nodes(:, 14))
       nodes_c(:, 15) = nodes_b(:, 1)
       nodes_d(:, 1) = nodes_a(:, 15)
       nodes_d(:, 2) = nodes_b(:, 4)
       nodes_d(:, 3) = nodes_b(:, 3)
       nodes_d(:, 4) = nodes_b(:, 2)
       nodes_d(:, 5) = nodes_b(:, 1)
       nodes_d(:, 6) = 0.125_dp * ( &
            nodes(:, 6) + 3 * nodes(:, 10) + 3 * nodes(:, 13) + nodes(:, 15))
       nodes_d(:, 7) = 0.125_dp * ( &
            nodes(:, 7) + nodes(:, 10) + 2 * nodes(:, 11) + &
            2 * nodes(:, 13) + nodes(:, 14) + nodes(:, 15))
       nodes_d(:, 8) = 0.125_dp * ( &
            nodes(:, 8) + 2 * nodes(:, 11) + nodes(:, 12) + nodes(:, 13) + &
            2 * nodes(:, 14) + nodes(:, 15))
       nodes_d(:, 9) = 0.125_dp * ( &
            nodes(:, 9) + 3 * nodes(:, 12) + 3 * nodes(:, 14) + nodes(:, 15))
       nodes_d(:, 10) = 0.25_dp * ( &
            nodes(:, 10) + 2 * nodes(:, 13) + nodes(:, 15))
       nodes_d(:, 11) = 0.25_dp * ( &
            nodes(:, 11) + nodes(:, 13) + nodes(:, 14) + nodes(:, 15))
       nodes_d(:, 12) = 0.25_dp * ( &
            nodes(:, 12) + 2 * nodes(:, 14) + nodes(:, 15))
       nodes_d(:, 13) = 0.5_dp * (nodes(:, 13) + nodes(:, 15))
       nodes_d(:, 14) = 0.5_dp * (nodes(:, 14) + nodes(:, 15))
       nodes_d(:, 15) = nodes(:, 15)
    else
       call specialize_surface( &
            num_nodes, dimension_, nodes, degree, &
            [1.0_dp, 0.0_dp, 0.0_dp], &
            [0.5_dp, 0.5_dp, 0.0_dp], &
            [0.5_dp, 0.0_dp, 0.5_dp], &
            nodes_a)
       call specialize_surface( &
            num_nodes, dimension_, nodes, degree, &
            [0.0_dp, 0.5_dp, 0.5_dp], &
            [0.5_dp, 0.0_dp, 0.5_dp], &
            [0.5_dp, 0.5_dp, 0.0_dp], &
            nodes_b)
       call specialize_surface( &
            num_nodes, dimension_, nodes, degree, &
            [0.5_dp, 0.5_dp, 0.0_dp], &
            [0.0_dp, 1.0_dp, 0.0_dp], &
            [0.0_dp, 0.5_dp, 0.5_dp], &
            nodes_c)
       call specialize_surface( &
            num_nodes, dimension_, nodes, degree, &
            [0.5_dp, 0.0_dp, 0.5_dp], &
            [0.0_dp, 0.5_dp, 0.5_dp], &
            [0.0_dp, 0.0_dp, 1.0_dp], &
            nodes_d)
    end if

  end subroutine subdivide_nodes

  subroutine compute_edge_nodes( &
       num_nodes, dimension_, nodes, degree, nodes1, nodes2, nodes3) &
       bind(c, name='compute_edge_nodes')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(out) :: nodes1(dimension_, degree + 1)
    real(c_double), intent(out) :: nodes2(dimension_, degree + 1)
    real(c_double), intent(out) :: nodes3(dimension_, degree + 1)
    ! Variables outside of signature.
    integer(c_int) :: index1, index2, index3

    index2 = degree + 1
    index3 = num_nodes
    do index1 = 1, degree + 1
       nodes1(:, index1) = nodes(:, index1)
       nodes2(:, index1) = nodes(:, index2)
       nodes3(:, index1) = nodes(:, index3)
       ! Update the indices.
       index2 = index2 - index1 + degree + 1
       index3 = index3 - index1 - 1
    end do

  end subroutine compute_edge_nodes

end module surface
