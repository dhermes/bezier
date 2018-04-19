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

module curve

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_bool
  use types, only: dp
  use helpers, only: cross_product, contains_nd
  implicit none
  private &
       MAX_LOCATE_SUBDIVISIONS, LOCATE_STD_CAP, &
       SQRT_PREC, REDUCE_THRESHOLD, scalar_func, dqagse, &
       specialize_curve_generic, specialize_curve_quadratic, &
       subdivide_nodes_generic, split_candidate, allocate_candidates, &
       update_candidates, projection_error, can_reduce
  public &
       CurveData, LOCATE_MISS, LOCATE_INVALID, evaluate_curve_barycentric, &
       evaluate_multi, specialize_curve, evaluate_hodograph, subdivide_nodes, &
       newton_refine, locate_point, elevate_nodes, get_curvature, &
       reduce_pseudo_inverse, full_reduce, compute_length, curves_equal, &
       subdivide_curve

  ! NOTE: This (for now) is not meant to be C-interoperable. This is mostly
  !       because the shape is encoded in ``nodes``, so it would be wasteful to
  !       also store ``num_nodes`` and ``dimension``. In addition,
  !       ``allocatable`` (rather than ``type(c_ptr)``) is used for the
  !       ``nodes`` so that the data is scoped with the ``CurveData`` instance
  !       and we don't need to worry about memory management.
  type :: CurveData
     real(c_double) :: start = 0.0_dp
     real(c_double) :: end_ = 1.0_dp
     real(c_double), allocatable :: nodes(:, :)
  end type CurveData

  ! NOTE: These values are also defined in equivalent Python source.
  integer(c_int), parameter :: MAX_LOCATE_SUBDIVISIONS = 20
  real(c_double), parameter :: LOCATE_STD_CAP = 0.5_dp**20
  ! NOTE: Should probably use ``d1mach`` to determine ``SQRT_PREC``.
  real(c_double), parameter :: SQRT_PREC = 0.5_dp**26
  real(c_double), parameter :: REDUCE_THRESHOLD = SQRT_PREC
  real(c_double), parameter :: LOCATE_MISS = -1
  real(c_double), parameter :: LOCATE_INVALID = -2

  ! Interface blocks for QUADPACK:dqagse
  abstract interface
     ! f: real(c_double) --> real(c_double)
     real(c_double) function scalar_func(x)
       use, intrinsic :: iso_c_binding, only: c_double
       implicit none

       real(c_double), intent(in) :: x
     end function scalar_func
  end interface

  interface
     ! D - double precision
     ! Q - quadrature
     ! A - adaptive
     ! G - General integrand (i.e. INT f(x), not weighted INT w(x) f(x))
     ! S - Singularities handled
     ! E - Extended
     ! See: https://en.wikipedia.org/wiki/QUADPACK
     ! QUADPACK is "Public Domain"
     subroutine dqagse( &
          f, a, b, epsabs, epsrel, limit, result_, &
          abserr, neval, ier, alist, blist, rlist, &
          elist, iord, last)
       use, intrinsic :: iso_c_binding, only: c_double, c_int
       implicit none

       procedure(scalar_func) :: f
       real(c_double), intent(in) :: a, b
       real(c_double), intent(in) :: epsabs, epsrel
       integer(c_int), intent(in) :: limit
       real(c_double), intent(out) :: result_, abserr
       integer(c_int), intent(out) :: neval, ier
       real(c_double), intent(out) :: alist(limit)
       real(c_double), intent(out) :: blist(limit)
       real(c_double), intent(out) :: rlist(limit)
       real(c_double), intent(out) :: elist(limit)
       integer(c_int), intent(out) :: iord(limit)
       integer(c_int), intent(out) :: last

     end subroutine dqagse
  end interface

contains

  subroutine evaluate_curve_barycentric( &
       num_nodes, dimension_, nodes, num_vals, lambda1, lambda2, evaluated) &
       bind(c, name='evaluate_curve_barycentric')

    ! NOTE: This is evaluate_multi_barycentric for a Bezier curve.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: num_vals
    real(c_double), intent(in) :: lambda1(num_vals)
    real(c_double), intent(in) :: lambda2(num_vals)
    real(c_double), intent(out) :: evaluated(dimension_, num_vals)
    ! Variables outside of signature.
    integer(c_int) :: i, j
    real(c_double) :: lambda2_pow(num_vals)
    integer(c_int) :: binom_val

    lambda2_pow = 1.0_dp
    binom_val = 1

    forall (i = 1:num_vals)
       evaluated(:, i) = lambda1(i) * nodes(:, 1)
    end forall

    do i = 2, num_nodes - 1
       lambda2_pow = lambda2_pow * lambda2
       binom_val = (binom_val * (num_nodes - i + 1)) / (i - 1)
       forall (j = 1:num_vals)
          evaluated(:, j) = ( &
               evaluated(:, j) + &
               binom_val * lambda2_pow(j) * nodes(:, i)) * lambda1(j)
       end forall
    end do

    forall (i = 1:num_vals)
       evaluated(:, i) = ( &
            evaluated(:, i) + &
            lambda2_pow(i) * lambda2(i) * nodes(:, num_nodes))
    end forall

  end subroutine evaluate_curve_barycentric

  subroutine evaluate_multi( &
       num_nodes, dimension_, nodes, num_vals, s_vals, evaluated) &
       bind(c, name='evaluate_multi')

    ! NOTE: This is evaluate_multi for a Bezier curve.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(in) :: num_vals
    real(c_double), intent(in) :: s_vals(num_vals)
    real(c_double), intent(out) :: evaluated(dimension_, num_vals)
    ! Variables outside of signature.
    real(c_double) :: one_less(num_vals)

    one_less = 1.0_dp - s_vals
    call evaluate_curve_barycentric( &
         num_nodes, dimension_, nodes, num_vals, one_less, s_vals, evaluated)
  end subroutine evaluate_multi

  subroutine specialize_curve_generic( &
       num_nodes, dimension_, nodes, start, end_, new_nodes)

    ! NOTE: This is a helper for ``specialize_curve`` that works on any degree.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(in) :: start, end_
    real(c_double), intent(out) :: new_nodes(dimension_, num_nodes)
    ! Variables outside of signature.
    real(c_double) :: workspace(dimension_, num_nodes - 1, num_nodes)
    integer(c_int) :: index_, curr_size, j
    real(c_double) :: minus_start, minus_end

    minus_start = 1.0_dp - start
    minus_end = 1.0_dp - end_
    workspace(:, :, 1) = ( &
         minus_start * nodes(:, :num_nodes - 1) + start * nodes(:, 2:))
    workspace(:, :, 2) = ( &
         minus_end * nodes(:, :num_nodes - 1) + end_ * nodes(:, 2:))

    curr_size = num_nodes - 1
    do index_ = 3, num_nodes
       curr_size = curr_size - 1
       ! First add a new "column" (or whatever the 3rd dimension is called)
       ! at the end using ``end_``.
       workspace(:, :curr_size, index_) = ( &
            minus_end * workspace(:, :curr_size, index_ - 1) + &
            end_ * workspace(:, 2:curr_size + 1, index_ - 1))
       ! Update all the values in place by using de Casteljau with the
       ! ``start`` parameter.
       forall (j = 1:index_ - 1)
          workspace(:, :curr_size, j) = ( &
               minus_start * workspace(:, :curr_size, j) + &
               start * workspace(:, 2:curr_size + 1, j))
       end forall
    end do

    ! Move the final "column" (or whatever the 3rd dimension is called)
    ! of the workspace into ``new_nodes``.
    forall (index_ = 1:num_nodes)
       new_nodes(:, index_) = workspace(:, 1, index_)
    end forall

  end subroutine specialize_curve_generic

  subroutine specialize_curve_quadratic( &
       dimension_, nodes, start, end_, new_nodes)

    integer(c_int), intent(in) :: dimension_
    real(c_double), intent(in) :: nodes(dimension_, 3)
    real(c_double), intent(in) :: start, end_
    real(c_double), intent(out) :: new_nodes(dimension_, 3)
    ! Variables outside of signature.
    real(c_double) :: minus_start, minus_end, prod_both

    minus_start = 1.0_dp - start
    minus_end = 1.0_dp - end_
    prod_both = start * end_

    new_nodes(:, 1) = ( &
         minus_start * minus_start * nodes(:, 1) + &
         2.0_dp * start * minus_start * nodes(:, 2) + &
         start * start * nodes(:, 3))
    new_nodes(:, 2) = ( &
         minus_start * minus_end * nodes(:, 1) + &
         (end_ + start - 2.0_dp * prod_both) * nodes(:, 2) + &
         prod_both * nodes(:, 3))
    new_nodes(:, 3) = ( &
         minus_end * minus_end * nodes(:, 1) + &
         2.0_dp * end_ * minus_end * nodes(:, 2) + &
         end_ * end_ * nodes(:, 3))

  end subroutine specialize_curve_quadratic

  subroutine specialize_curve( &
       num_nodes, dimension_, nodes, start, end_, new_nodes) &
       bind(c, name='specialize_curve')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(in) :: start, end_
    real(c_double), intent(out) :: new_nodes(dimension_, num_nodes)

    if (num_nodes == 2) then
       new_nodes(:, 1) = (1.0_dp - start) * nodes(:, 1) + start * nodes(:, 2)
       new_nodes(:, 2) = (1.0_dp - end_) * nodes(:, 1) + end_ * nodes(:, 2)
    else if (num_nodes == 3) then
       call specialize_curve_quadratic( &
            dimension_, nodes, start, end_, new_nodes)
    else
       call specialize_curve_generic( &
            num_nodes, dimension_, nodes, start, end_, new_nodes)
    end if

  end subroutine specialize_curve

  subroutine evaluate_hodograph( &
       s, num_nodes, dimension_, nodes, hodograph) &
       bind(c, name='evaluate_hodograph')

    real(c_double), intent(in) :: s
    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: hodograph(dimension_, 1)
    ! Variables outside of signature.
    real(c_double) :: first_deriv(dimension_, num_nodes - 1)

    first_deriv = nodes(:, 2:) - nodes(:, :num_nodes - 1)
    call evaluate_multi( &
         num_nodes - 1, dimension_, first_deriv, 1, [s], hodograph)
    hodograph = (num_nodes - 1) * hodograph

  end subroutine evaluate_hodograph

  subroutine subdivide_nodes_generic( &
       num_nodes, dimension_, nodes, left_nodes, right_nodes)

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: left_nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: right_nodes(dimension_, num_nodes)
    ! Variables outside of signature.
    real(c_double) :: pascals_triangle(num_nodes)
    integer(c_int) :: elt_index, pascal_index

    pascals_triangle = 0  ! Make sure all zero.
    pascals_triangle(1) = 1

    do elt_index = 1, num_nodes
       ! Update Pascal's triangle (intentionally at beginning, not end).
       if (elt_index > 1) then
          pascals_triangle(:elt_index) = 0.5_dp * ( &
               pascals_triangle(:elt_index) + pascals_triangle(elt_index:1:-1))
       end if

       left_nodes(:, elt_index) = 0
       right_nodes(:, num_nodes + 1 - elt_index) = 0
       do pascal_index = 1, elt_index
          left_nodes(:, elt_index) = ( &
               left_nodes(:, elt_index) + &
               pascals_triangle(pascal_index) * nodes(:, pascal_index))
          right_nodes(:, num_nodes + 1 - elt_index) = ( &
               right_nodes(:, num_nodes + 1 - elt_index) + &
               pascals_triangle(pascal_index) * &
               nodes(:, num_nodes + 1 - pascal_index))
       end do
    end do

  end subroutine subdivide_nodes_generic

  subroutine subdivide_nodes( &
       num_nodes, dimension_, nodes, left_nodes, right_nodes) &
       bind(c, name='subdivide_nodes_curve')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: left_nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: right_nodes(dimension_, num_nodes)

    if (num_nodes == 2) then
       left_nodes(:, 1) = nodes(:, 1)
       left_nodes(:, 2) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
       right_nodes(:, 1) = left_nodes(:, 2)
       right_nodes(:, 2) = nodes(:, 2)
    else if (num_nodes == 3) then
       left_nodes(:, 1) = nodes(:, 1)
       left_nodes(:, 2) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
       left_nodes(:, 3) = 0.25_dp * ( &
            nodes(:, 1) + 2 * nodes(:, 2) + nodes(:, 3))
       right_nodes(:, 1) = left_nodes(:, 3)
       right_nodes(:, 2) = 0.5_dp * (nodes(:, 2) + nodes(:, 3))
       right_nodes(:, 3) = nodes(:, 3)
    else if (num_nodes == 4) then
       left_nodes(:, 1) = nodes(:, 1)
       left_nodes(:, 2) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
       left_nodes(:, 3) = 0.25_dp * ( &
            nodes(:, 1) + 2 * nodes(:, 2) + nodes(:, 3))
       left_nodes(:, 4) = 0.125_dp * ( &
            nodes(:, 1) + 3 * nodes(:, 2) + 3 * nodes(:, 3) + nodes(:, 4))
       right_nodes(:, 1) = left_nodes(:, 4)
       right_nodes(:, 2) = 0.25_dp * ( &
            nodes(:, 2) + 2 * nodes(:, 3) + nodes(:, 4))
       right_nodes(:, 3) = 0.5_dp * (nodes(:, 3) + nodes(:, 4))
       right_nodes(:, 4) = nodes(:, 4)
    else
       call subdivide_nodes_generic( &
            num_nodes, dimension_, nodes, left_nodes, right_nodes)
    end if

  end subroutine subdivide_nodes

  subroutine newton_refine( &
       num_nodes, dimension_, nodes, point, s, updated_s) &
       bind(c, name='newton_refine_curve')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(in) :: point(dimension_, 1)
    real(c_double), intent(in) :: s
    real(c_double), intent(out) :: updated_s
    ! Variables outside of signature.
    real(c_double) :: pt_delta(dimension_, 1)
    real(c_double) :: derivative(dimension_, 1)

    call evaluate_multi( &
         num_nodes, dimension_, nodes, 1, [s], pt_delta)
    pt_delta = point - pt_delta
    ! At this point `pt_delta` is `p - B(s)`.
    call evaluate_hodograph( &
         s, num_nodes, dimension_, nodes, derivative)

    updated_s = ( &
         s + &
         (dot_product(pt_delta(:, 1), derivative(:, 1)) / &
         dot_product(derivative(:, 1), derivative(:, 1))))

  end subroutine newton_refine

  subroutine split_candidate( &
       num_nodes, dimension_, candidate, num_next_candidates, next_candidates)

    integer(c_int), intent(in) :: num_nodes, dimension_
    type(CurveData), intent(in) :: candidate
    integer(c_int), intent(in) :: num_next_candidates
    type(CurveData), intent(inout) :: next_candidates(:)

    ! Allocate the new nodes and call sub-divide.
    ! NOTE: This **assumes** but does not check that if the nodes are
    !       allocated, they are also the correct shape.
    if (.NOT. allocated(next_candidates(num_next_candidates - 1)%nodes)) then
       allocate(next_candidates(num_next_candidates - 1)%nodes(dimension_, num_nodes))
    end if
    if (.NOT. allocated(next_candidates(num_next_candidates)%nodes)) then
       allocate(next_candidates(num_next_candidates)%nodes(dimension_, num_nodes))
    end if

    call subdivide_nodes( &
         num_nodes, dimension_, candidate%nodes, &
         next_candidates(num_next_candidates - 1)%nodes, &
         next_candidates(num_next_candidates)%nodes)

    ! Left half.
    next_candidates(num_next_candidates - 1)%start = candidate%start
    next_candidates(num_next_candidates - 1)%end_ = ( &
         0.5_dp * (candidate%start + candidate%end_))
    ! Right half.
    next_candidates(num_next_candidates)%start = ( &
         next_candidates(num_next_candidates - 1)%end_)
    next_candidates(num_next_candidates)%end_ = candidate%end_

  end subroutine split_candidate

  subroutine allocate_candidates(num_candidates, next_candidates)
    integer(c_int), intent(in) :: num_candidates
    type(CurveData), allocatable, intent(inout) :: next_candidates(:)

    if (allocated(next_candidates)) then
       if (size(next_candidates) < 2 * num_candidates) then
          deallocate(next_candidates)
          allocate(next_candidates(2 * num_candidates))
       end if
    else
       allocate(next_candidates(2 * num_candidates))
    end if

  end subroutine allocate_candidates

  subroutine update_candidates( &
       num_nodes, dimension_, point, num_candidates, candidates, &
       num_next_candidates, next_candidates)

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: point(dimension_, 1)
    integer(c_int), intent(in) :: num_candidates
    type(CurveData), intent(in) :: candidates(:)
    integer(c_int), intent(inout) :: num_next_candidates
    type(CurveData), allocatable, intent(inout) :: next_candidates(:)
    ! Variables outside of signature.
    integer(c_int) :: cand_index
    logical(c_bool) :: predicate

    ! Allocate maximum amount of space needed.
    call allocate_candidates(num_candidates, next_candidates)

    do cand_index = 1, num_candidates
       call contains_nd( &
            num_nodes, dimension_, candidates(cand_index)%nodes, &
            point(:, 1), predicate)
       if (predicate) then
          num_next_candidates = num_next_candidates + 2
          call split_candidate( &
               num_nodes, dimension_, candidates(cand_index), &
               num_next_candidates, next_candidates)
       end if
    end do

  end subroutine update_candidates

  subroutine locate_point( &
       num_nodes, dimension_, nodes, point, s_approx) &
       bind(c, name='locate_point_curve')

    ! NOTE: This returns ``-1`` (``LOCATE_MISS``) as a signal for "point is
    !       not on the curve" and ``-2`` (``LOCATE_INVALID``) for "point is
    !       on separate segments" (i.e. the standard deviation of the
    !       parameters is too large).

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(in) :: point(dimension_, 1)
    real(c_double), intent(out) :: s_approx
    ! Variables outside of signature.
    integer(c_int) :: sub_index
    integer(c_int) :: num_candidates, num_next_candidates
    type(CurveData), allocatable :: candidates_odd(:), candidates_even(:)
    real(c_double), allocatable :: s_params(:)
    real(c_double) :: std_dev
    logical(c_bool) :: is_even

    ! Start out with the full curve.
    allocate(candidates_odd(1))  ! First iteration is odd.
    candidates_odd(1) = CurveData(0.0_dp, 1.0_dp, nodes)
    ! NOTE: `num_candidates` will be tracked separately
    !       from `size(candidates)`.
    num_candidates = 1
    s_approx = LOCATE_MISS

    is_even = .TRUE.  ! At zero.
    do sub_index = 1, MAX_LOCATE_SUBDIVISIONS + 1
       is_even = .NOT. is_even  ! Switch parity.
       num_next_candidates = 0

       if (is_even) then
          ! Read from even, write to odd.
          call update_candidates( &
               num_nodes, dimension_, point, num_candidates, candidates_even, &
               num_next_candidates, candidates_odd)
       else
          ! Read from odd, write to even.
          call update_candidates( &
               num_nodes, dimension_, point, num_candidates, candidates_odd, &
               num_next_candidates, candidates_even)
       end if

       ! If there are no more candidates, we are done.
       if (num_next_candidates == 0) then
          return
       end if

       num_candidates = num_next_candidates

    end do

    ! Compute the s-parameter as the mean of the **start** and
    ! **end** parameters.
    allocate(s_params(2 * num_candidates))
    if (is_even) then
       ! NOTE: We "exclude" this block from ``lcov`` because it can never
       !       be invoked while ``MAX_LOCATE_SUBDIVISIONS`` is even.
       ! LCOV_EXCL_START
       s_params(:num_candidates) = ( &
            candidates_odd(:num_candidates)%start)
       s_params(num_candidates + 1:) = ( &
            candidates_odd(:num_candidates)%end_)
       ! LCOV_EXCL_STOP
    else
       s_params(:num_candidates) = candidates_even(:num_candidates)%start
       s_params(num_candidates + 1:) = candidates_even(:num_candidates)%end_
    end if
    s_approx = sum(s_params) / (2 * num_candidates)

    std_dev = sqrt(sum((s_params - s_approx)**2) / (2 * num_candidates))
    if (std_dev > LOCATE_STD_CAP) then
       s_approx = LOCATE_INVALID
       return
    end if

    ! NOTE: Use ``std_dev`` variable as a "placeholder" for the update.
    call newton_refine( &
         num_nodes, dimension_, nodes, point, s_approx, std_dev)
    ! NOTE: Since ``mean(s_params)`` must be in ``[0, 1]`` it's
    !       "safe" to push the Newton-refined value back into the unit
    !       interval.
    if (std_dev < 0.0_dp) then
       s_approx = 0.0_dp
    else if (std_dev > 1.0_dp) then
       s_approx = 1.0_dp
    else
       s_approx = std_dev
    end if

  end subroutine locate_point

  subroutine elevate_nodes( &
       num_nodes, dimension_, nodes, elevated) &
       bind(c, name='elevate_nodes_curve')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: elevated(dimension_, num_nodes + 1)
    ! Variables outside of signature.
    integer(c_int) :: i

    elevated(:, 1) = nodes(:, 1)
    forall (i = 1:num_nodes - 1)
       elevated(:, i + 1) = ( &
            i * nodes(:, i) + (num_nodes - i) * nodes(:, i + 1)) / num_nodes
    end forall
    elevated(:, num_nodes + 1) = nodes(:, num_nodes)

  end subroutine elevate_nodes

  subroutine get_curvature( &
       num_nodes, nodes, tangent_vec, s, curvature) &
       bind(c, name='get_curvature')

    ! NOTE: This **only** computes curvature for plane curves (i.e. curves
    !       in R^2). An equivalent notion of curvature exists for space
    !       curves, but support for that is not implemented here.

    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: nodes(2, num_nodes)
    real(c_double), intent(in) :: tangent_vec(2, 1)
    real(c_double), intent(in) :: s
    real(c_double), intent(out) :: curvature
    ! Variables outside of signature.
    real(c_double) :: work(2, num_nodes - 1)
    real(c_double) :: concavity(2, 1)

    if (num_nodes == 2) then
       curvature = 0
       return
    end if

    ! NOTE: We somewhat replicate code in ``evaluate_hodograph()`` here.

    ! First derivative:
    work = nodes(:, 2:) - nodes(:, :num_nodes - 1)
    ! Second derivative (no need for last element of work array):
    work(:, :num_nodes - 2) = work(:, 2:) - work(:, :num_nodes - 2)

    ! NOTE: The degree being evaluated is ``degree - 2 == num_nodes - 3``.
    call evaluate_multi( &
         num_nodes - 2, 2, work(:, :num_nodes - 2), 1, [s], concavity)
    ! B''(s) = d (d - 1) D(s) where D(s) is defined by the "double hodograph".
    concavity = concavity * (num_nodes - 1) * (num_nodes - 2)

    call cross_product(tangent_vec, concavity, curvature)
    curvature = curvature / norm2(tangent_vec)**3

  end subroutine get_curvature

  subroutine reduce_pseudo_inverse( &
       num_nodes, dimension_, nodes, reduced, not_implemented) &
       bind(c, name='reduce_pseudo_inverse')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: reduced(dimension_, num_nodes - 1)
    logical(c_bool), intent(out) :: not_implemented

    not_implemented = .FALSE.
    if (num_nodes == 2) then
       reduced(:, 1) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
    else if (num_nodes == 3) then
       reduced(:, 1) = (5 * nodes(:, 1) + 2 * nodes(:, 2) - nodes(:, 3)) / 6
       reduced(:, 2) = (-nodes(:, 1) + 2 * nodes(:, 2) + 5 * nodes(:, 3)) / 6
    else if (num_nodes == 4) then
       reduced(:, 1) = ( &
            19 * nodes(:, 1) + 3 * nodes(:, 2) - &
            3 * nodes(:, 3) + nodes(:, 4)) / 20
       reduced(:, 2) = 0.25_dp * ( &
            -nodes(:, 1) + 3 * nodes(:, 2) + &
            3 * nodes(:, 3) - nodes(:, 4))
       reduced(:, 3) = ( &
            nodes(:, 1) - 3 * nodes(:, 2) + &
            3 * nodes(:, 3) + 19 * nodes(:, 4)) / 20
    else if (num_nodes == 5) then
       reduced(:, 1) = ( &
            69 * nodes(:, 1) + 4 * nodes(:, 2) - 6 * nodes(:, 3) + &
            4 * nodes(:, 4) - nodes(:, 5)) / 70
       reduced(:, 2) = ( &
            -53 * nodes(:, 1) + 212 * nodes(:, 2) + 102 * nodes(:, 3) - &
            68 * nodes(:, 4) + 17 * nodes(:, 5)) / 210
       reduced(:, 3) = ( &
            17 * nodes(:, 1) - 68 * nodes(:, 2) + 102 * nodes(:, 3) + &
            212 * nodes(:, 4) - 53 * nodes(:, 5)) / 210
       reduced(:, 4) = ( &
            -nodes(:, 1) + 4 * nodes(:, 2) - 6 * nodes(:, 3) + &
            4 * nodes(:, 4) + 69 * nodes(:, 5)) / 70
    else
       not_implemented = .TRUE.
    end if

  end subroutine reduce_pseudo_inverse

  subroutine projection_error( &
       num_nodes, dimension_, nodes, projected, error)

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(in) :: projected(dimension_, num_nodes)
    real(c_double), intent(out) :: error

    ! If "dim" is not passed to ``norm2``, will be Frobenius norm.
    error = norm2(nodes - projected)
    if (error == 0.0_dp) then
       return
    end if

    ! Make the error relative (in Frobenius norm).
    error = error / norm2(nodes)

  end subroutine projection_error

  subroutine can_reduce( &
       num_nodes, dimension_, nodes, success)

    ! NOTE: This returns ``success = 0`` for "Failure", ``success = 1`` for
    !       "Success" and ``success = -1`` for "Not Implemented".
    ! NOTE: This subroutine is not part of the C ABI for this module,
    !       but it is (for now) public, so that it can be tested.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(out) :: success
    ! Variables outside of signature.
    real(c_double) :: reduced(dimension_, num_nodes)
    real(c_double) :: relative_err

    if (num_nodes < 2 .OR. num_nodes > 5) then
       ! The only caller `full_reduce` will never try to reduce a point
       ! (i.e. `num_nodes == `). So these cases (along with higher degre
       ! than quartics) are "Not Implemented".
       success = -1
       return
    end if

    ! First, put the "projection" in ``reduced``.
    if (num_nodes == 2) then
       reduced(:, 1) = 0.5_dp * (nodes(:, 1) + nodes(:, 2))
       reduced(:, 2) = reduced(:, 1)
    else if (num_nodes == 3) then
       reduced(:, 1) = (5 * nodes(:, 1) + 2 * nodes(:, 2) - nodes(:, 3)) / 6
       reduced(:, 2) = (nodes(:, 1) + nodes(:, 2) + nodes(:, 3)) / 3
       reduced(:, 3) = (-nodes(:, 1) + 2 * nodes(:, 2) + 5 * nodes(:, 3)) / 6
    else if (num_nodes == 4) then
       reduced(:, 1) = ( &
            19 * nodes(:, 1) + 3 * nodes(:, 2) - &
            3 * nodes(:, 3) + nodes(:, 4)) / 20
       reduced(:, 2) = ( &
            3 * nodes(:, 1) + 11 * nodes(:, 2) + &
            9 * nodes(:, 3) - 3 * nodes(:, 4)) / 20
       reduced(:, 3) = ( &
            -3 * nodes(:, 1) + 9 * nodes(:, 2) + &
            11 * nodes(:, 3) + 3 * nodes(:, 4)) / 20
       reduced(:, 4) = ( &
            nodes(:, 1) - 3 * nodes(:, 2) + &
            3 * nodes(:, 3) + 19 * nodes(:, 4)) / 20
    else if (num_nodes == 5) then
       reduced(:, 1) = ( &
            69 * nodes(:, 1) + 4 * nodes(:, 2) - 6 * nodes(:, 3) + &
            4 * nodes(:, 4) - nodes(:, 5)) / 70
       reduced(:, 2) = ( &
            2 * nodes(:, 1) + 27 * nodes(:, 2) + 12 * nodes(:, 3) - &
            8 * nodes(:, 4) + 2 * nodes(:, 5)) / 35
       reduced(:, 3) = ( &
            -3 * nodes(:, 1) + 12 * nodes(:, 2) + 17 * nodes(:, 3) + &
            12 * nodes(:, 4) - 3 * nodes(:, 5)) / 35
       reduced(:, 4) = ( &
            2 * nodes(:, 1) - 8 * nodes(:, 2) + 12 * nodes(:, 3) + &
            27 * nodes(:, 4) + 2 * nodes(:, 5)) / 35
       reduced(:, 5) = ( &
            -nodes(:, 1) + 4 * nodes(:, 2) - 6 * nodes(:, 3) + &
            4 * nodes(:, 4) + 69 * nodes(:, 5)) / 70
    end if

    call projection_error(num_nodes, dimension_, nodes, reduced, relative_err)
    if (relative_err < REDUCE_THRESHOLD) then
       success = 1
    else
       success = 0
    end if

  end subroutine can_reduce

  subroutine full_reduce( &
       num_nodes, dimension_, nodes, num_reduced_nodes, &
       reduced, not_implemented) &
       bind(c, name='full_reduce')

    ! NOTE: The size of ``reduced`` represents the **maximum** possible
    !       size, but ``num_reduced_nodes`` actually reflects the number
    !       of nodes in the fully reduced nodes.

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    integer(c_int), intent(out) :: num_reduced_nodes
    real(c_double), intent(out) :: reduced(dimension_, num_nodes)
    logical(c_bool), intent(out) :: not_implemented
    ! Variables outside of signature.
    integer(c_int) :: i, cr_success
    real(c_double) :: work(dimension_, num_nodes - 1)

    reduced = nodes
    num_reduced_nodes = num_nodes
    not_implemented = .FALSE.
    ! We can make at most ``num_nodes - 1`` reductions since
    ! we can't reduce past one node (i.e. degree zero).
    do i = 1, num_nodes - 1
       call can_reduce( &
            num_reduced_nodes, dimension_, &
            reduced(:, :num_reduced_nodes), cr_success)

       if (cr_success == 1) then
          ! Actually reduce the nodes.
          call reduce_pseudo_inverse( &
               num_reduced_nodes, dimension_, reduced(:, :num_reduced_nodes), &
               work(:, :num_reduced_nodes - 1), not_implemented)
          if (not_implemented) then
             return  ! LCOV_EXCL_LINE
          else
             num_reduced_nodes = num_reduced_nodes - 1
             ! Update `reduced` based on the **new** number of nodes.
             reduced(:, :num_reduced_nodes) = work(:, :num_reduced_nodes)
          end if
       else if (cr_success == 0) then
          return
       else
          ! ``cr_success == -1`` means "Not Implemented"
          not_implemented = .TRUE.
          return
       end if
    end do

  end subroutine full_reduce

  subroutine compute_length( &
       num_nodes, dimension_, nodes, length, error_val) &
       bind(c, name='compute_length')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: length
    integer(c_int), intent(out) :: error_val
    ! Variables outside of signature.
    real(c_double) :: first_deriv(dimension_, num_nodes - 1)
    real(c_double) :: abserr
    integer(c_int) :: neval
    real(c_double) :: alist(50)
    real(c_double) :: blist(50)
    real(c_double) :: rlist(50)
    real(c_double) :: elist(50)
    integer(c_int) :: iord(50)
    integer(c_int) :: last

    ! NOTE: We somewhat replicate code in ``evaluate_hodograph()``
    !       here. This is so we don't re-compute the nodes for the first
    !       derivative every time it is evaluated.
    first_deriv = (num_nodes - 1) * (nodes(:, 2:) - nodes(:, :num_nodes - 1))
    if (num_nodes == 2) then
       length = norm2(first_deriv)
       error_val = 0
       return
    end if

    call dqagse( &
         vec_size, 0.0_dp, 1.0_dp, SQRT_PREC, SQRT_PREC, 50, length, &
         abserr, neval, error_val, alist, blist, rlist, &
         elist, iord, last)

  contains

    ! Define a closure that evaluates ||B'(s)||_2 where ``s``
    ! is the argument and ``B'(s)`` is parameterized by ``first_deriv``.
    real(c_double) function vec_size(s_val) result(norm_)
      real(c_double), intent(in) :: s_val
      ! Variables outside of signature.
      real(c_double) :: evaluated(dimension_, 1)

      ! ``evaluate_multi`` takes degree, which is one less than the number
      ! of nodes, so our derivative is one less than that.
      call evaluate_multi( &
           num_nodes - 1, dimension_, first_deriv, 1, [s_val], evaluated)
      norm_ = norm2(evaluated)

    end function vec_size

  end subroutine compute_length

  logical(c_bool) function curves_equal(curve1, curve2) result(same)

    ! NOTE: This is **explicitly** not intended for C inter-op.
    ! NOTE: This is really intended to be used by diagnostic / test code
    !       and not during the course of computation.

    type(CurveData), intent(in) :: curve1, curve2

    same = .TRUE.

    ! First, check the scalars.
    if (curve1%start /= curve2%start .OR. curve1%end_ /= curve2%end_) then
       same = .FALSE.
       return
    end if

    ! Then, check the nodes.
    if (allocated(curve1%nodes)) then
       if (allocated(curve2%nodes)) then
          if (all(shape(curve1%nodes) == shape(curve2%nodes))) then
             ! Here, since the same shape, we only exit early
             ! if the nodes aren't identical.
             if (.NOT. all(curve1%nodes == curve2%nodes)) then
                same = .FALSE.
                return
             end if
          else
             ! Different shapes can't be equal.
             same = .FALSE.
             return
          end if
       else
          ! nodes1 are allocated, nodes2 are not, can't be equal.
          same = .FALSE.
          return
       end if
    else
       if (allocated(curve2%nodes)) then
          ! nodes2 are allocated, nodes1 are not, can't be equal.
          same = .FALSE.
          return
       end if
       ! In the implicit ``else`` branch here, both sets of nots are
       ! un-allocated, hence they are the "same".
    end if

  end function curves_equal

  subroutine subdivide_curve(curve_data, left, right)

    ! NOTE: This is **explicitly** not intended for C inter-op.
    ! NOTE: This **assumes** but does not check that ``curve_data%nodes``
    !       is allocated.

    type(CurveData), intent(in) :: curve_data
    type(CurveData), intent(out) :: left, right
    ! Variables outside of signature.
    integer(c_int) :: num_nodes, dimension_

    num_nodes = size(curve_data%nodes, 2)
    dimension_ = size(curve_data%nodes, 1)

    left%start = curve_data%start
    left%end_ = 0.5_dp * (curve_data%start + curve_data%end_)
    allocate(left%nodes(dimension_, num_nodes))

    right%start = left%end_
    right%end_ = curve_data%end_
    allocate(right%nodes(dimension_, num_nodes))

    call subdivide_nodes( &
         num_nodes, dimension_, curve_data%nodes, left%nodes, right%nodes)

  end subroutine subdivide_curve

end module curve
