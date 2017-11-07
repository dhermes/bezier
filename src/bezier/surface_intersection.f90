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

module surface_intersection

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_bool
  use curve, only: LOCATE_MISS
  use helpers, only: contains_nd, vector_close
  use types, only: dp
  use surface, only: evaluate_barycentric, jacobian_both, subdivide_nodes
  implicit none
  private &
       LocateCandidate, MAX_LOCATE_SUBDIVISIONS, LOCATE_EPS, &
       newton_refine_solve, split_candidate
  public newton_refine, locate_point

  ! For ``locate_point``.
  type :: LocateCandidate
     real(c_double) :: centroid_x = 1.0_dp  ! Actually triple.
     real(c_double) :: centroid_y = 1.0_dp  ! Actually triple.
     real(c_double) :: width = 1.0_dp
     real(c_double), allocatable :: nodes(:, :)
  end type LocateCandidate

  ! NOTE: These values are also defined in equivalent Python source.
  integer(c_int), parameter :: MAX_LOCATE_SUBDIVISIONS = 20
  real(c_double), parameter :: LOCATE_EPS = 0.5_dp**47

contains

  subroutine newton_refine_solve( &
       jac_both, x_val, surf_x, y_val, surf_y, delta_s, delta_t)

    ! NOTE: This is a helper for ``newton_refine``.

    real(c_double), intent(in) :: jac_both(1, 4)
    real(c_double), intent(in) :: x_val, surf_x
    real(c_double), intent(in) :: y_val, surf_y
    real(c_double), intent(out) :: delta_s, delta_t
    ! Variables outside of signature.
    real(c_double) :: e_val, f_val, denominator

    e_val = x_val - surf_x
    f_val = y_val - surf_y
    denominator = ( &
         jac_both(1, 1) * jac_both(1, 4) - jac_both(1, 2) * jac_both(1, 3))
    delta_s = (jac_both(1, 4) * e_val - jac_both(1, 3) * f_val) / denominator
    delta_t = (jac_both(1, 1) * f_val - jac_both(1, 2) * e_val) / denominator

  end subroutine newton_refine_solve

  subroutine newton_refine( &
       num_nodes, nodes, degree, x_val, y_val, &
       s, t, updated_s, updated_t) &
       bind(c, name='newton_refine_surface')

    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: nodes(num_nodes, 2)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(in) :: x_val, y_val
    real(c_double), intent(in) :: s, t
    real(c_double), intent(out) :: updated_s, updated_t
    ! Variables outside of signature.
    real(c_double) :: point(1, 2), jac_both(1, 4)
    real(c_double) :: jac_nodes(num_nodes - degree - 1, 4)
    real(c_double) :: delta_s, delta_t

    call evaluate_barycentric( &
         num_nodes, 2, nodes, degree, &
         1.0_dp - s - t, s, t, point)

    if (point(1, 1) == x_val .AND. point(1, 2) == y_val) then
       ! No refinement is needed.
       updated_s = s
       updated_t = t
       return
    end if

    call jacobian_both( &
         num_nodes, 2, nodes, degree, jac_nodes)
    call evaluate_barycentric( &
         num_nodes - degree - 1, 4, jac_nodes, degree - 1, &
         1.0_dp - s - t, s, t, jac_both)
    call newton_refine_solve( &
         jac_both, x_val, point(1, 1), &
         y_val, point(1, 2), delta_s, delta_t)
    updated_s = s + delta_s
    updated_t = t + delta_t

  end subroutine newton_refine

  subroutine split_candidate( &
       num_nodes, degree, candidate, num_next_candidates, next_candidates)

    ! NOTE: This assumes nodes are 2-dimensional.
    ! NOTE: This assumes that the nodes in each sub-candidate are
    !       not yet allocated.

    integer(c_int), intent(in) :: num_nodes
    integer(c_int), intent(in) :: degree
    type(LocateCandidate), intent(in) :: candidate
    integer(c_int), intent(in) :: num_next_candidates
    type(LocateCandidate), intent(inout) :: next_candidates(:)
    ! Variables outside of signature.
    real(c_double) :: half_width

    ! Allocate the new nodes and call sub-divide.
    allocate(next_candidates(num_next_candidates - 3)%nodes(num_nodes, 2))
    allocate(next_candidates(num_next_candidates - 2)%nodes(num_nodes, 2))
    allocate(next_candidates(num_next_candidates - 1)%nodes(num_nodes, 2))
    allocate(next_candidates(num_next_candidates)%nodes(num_nodes, 2))

    call subdivide_nodes( &
         num_nodes, 2, candidate%nodes, degree, &
         next_candidates(num_next_candidates - 3)%nodes, &
         next_candidates(num_next_candidates - 2)%nodes, &
         next_candidates(num_next_candidates - 1)%nodes, &
         next_candidates(num_next_candidates)%nodes)

    half_width = 0.5_dp * candidate%width

    ! Subdivision A.
    next_candidates(num_next_candidates - 3)%centroid_x = ( &
         candidate%centroid_x - half_width)
    next_candidates(num_next_candidates - 3)%centroid_y = ( &
         candidate%centroid_y - half_width)
    next_candidates(num_next_candidates - 3)%width = half_width

    ! Subdivision B.
    next_candidates(num_next_candidates - 2)%centroid_x = candidate%centroid_x
    next_candidates(num_next_candidates - 2)%centroid_y = candidate%centroid_y
    next_candidates(num_next_candidates - 2)%width = -half_width

    ! Subdivision C.
    next_candidates(num_next_candidates - 1)%centroid_x = ( &
         candidate%centroid_x + candidate%width)
    next_candidates(num_next_candidates - 1)%centroid_y = ( &
         next_candidates(num_next_candidates - 3)%centroid_y)
    next_candidates(num_next_candidates - 1)%width = half_width

    ! Subdivision D.
    next_candidates(num_next_candidates)%centroid_x = ( &
         next_candidates(num_next_candidates - 3)%centroid_x)
    next_candidates(num_next_candidates)%centroid_y = ( &
         candidate%centroid_y + candidate%width)
    next_candidates(num_next_candidates)%width = half_width

  end subroutine split_candidate

  subroutine locate_point( &
       num_nodes, nodes, degree, x_val, y_val, s_val, t_val) &
       bind(c, name='locate_point_surface')

    ! NOTE: This solves the inverse problem B(s, t) = (x, y) (if it can be
    !       solved). Does so by subdividing the surface until the sub-surfaces
    !       are sufficiently small, then using Newton's method to narrow
    !       in on the pre-image of the point.
    ! NOTE: This returns ``-1`` (``LOCATE_MISS``) for ``s_val`` as a signal
    !       for "point is not on the surface".
    ! NOTE: This assumes, but does not check, that the surface is "valid",
    !       i.e. all pre-image are unique.

    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: nodes(num_nodes, 2)
    integer(c_int), intent(in) :: degree
    real(c_double), intent(in) :: x_val, y_val
    real(c_double), intent(out) :: s_val, t_val
    ! Variables outside of signature.
    real(c_double) :: point(1, 2)
    integer(c_int) :: sub_index, cand_index
    integer(c_int) :: num_candidates, num_next_candidates
    type(LocateCandidate), allocatable :: candidates(:), next_candidates(:)
    logical(c_bool) :: predicate
    real(c_double) :: s_approx, t_approx
    real(c_double) :: actual(1, 2)

    point(1, :) = [x_val, y_val]
    ! Start out with the full curve.
    allocate(candidates(1))
    candidates(1) = LocateCandidate(1.0_dp, 1.0_dp, 1.0_dp, nodes)
    num_candidates = 1
    s_val = LOCATE_MISS

    do sub_index = 1, MAX_LOCATE_SUBDIVISIONS + 1
       num_next_candidates = 0
       ! Allocate maximum amount of space needed.
       allocate(next_candidates(4 * num_candidates))
       do cand_index = 1, num_candidates
          call contains_nd( &
               num_nodes, 2, candidates(cand_index)%nodes, &
               point(1, :), predicate)
          if (predicate) then
             num_next_candidates = num_next_candidates + 4
             call split_candidate( &
                  num_nodes, degree, candidates(cand_index), &
                  num_next_candidates, next_candidates)
          end if
       end do

       ! If there are no more candidates, we are done.
       if (num_next_candidates == 0) then
          return
       end if

       ! NOTE: This may copy empty slots, but this is OK since we track
       !       `num_candidates` separately.
       call move_alloc(next_candidates, candidates)
       num_candidates = num_next_candidates

    end do

    ! Compute the s- and t-parameters as the mean of the centroid positions.
    ! We divide by `3n` rather than `n` since we have tracked thrice the
    ! centroid rather than the centroid itself (to avoid round-off until
    ! right now).
    s_approx = ( &
         sum(candidates(:num_candidates)%centroid_x) / &
         (3.0_dp * num_candidates))
    t_approx = ( &
         sum(candidates(:num_candidates)%centroid_y) / &
         (3.0_dp * num_candidates))

    call newton_refine( &
         num_nodes, nodes, degree, x_val, y_val, &
         s_approx, t_approx, s_val, t_val)

    ! Check if the solution is "close enough" ...
    call evaluate_barycentric( &
         num_nodes, 2, nodes, degree, &
         1.0_dp - s_val - t_val, s_val, t_val, actual)

    if (vector_close(2, point, actual, LOCATE_EPS)) then
       return
    end if

    ! ... and if **not** close enough, do one more round of
    ! Newton's method if not.
    s_approx = s_val
    t_approx = t_val
    call newton_refine( &
         num_nodes, nodes, degree, x_val, y_val, &
         s_approx, t_approx, s_val, t_val)

  end subroutine locate_point

end module surface_intersection
