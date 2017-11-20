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
  use status, only: &
       Status_SUCCESS, Status_SAME_CURVATURE, Status_BAD_TANGENT, &
       Status_EDGE_END, Status_UNKNOWN
  use curve, only: CurveData, LOCATE_MISS, evaluate_hodograph, get_curvature
  use curve_intersection, only: &
       BoxIntersectionType_INTERSECTION, bbox_intersect, all_intersections
  use helpers, only: cross_product, contains_nd, vector_close
  use types, only: dp
  use surface, only: &
       evaluate_barycentric, jacobian_both, subdivide_nodes, compute_edge_nodes
  implicit none
  private &
       LocateCandidate, MAX_LOCATE_SUBDIVISIONS, LOCATE_EPS, &
       newton_refine_solve, split_candidate, allocate_candidates, &
       update_candidates, ignored_edge_corner, ignored_double_corner, &
       ignored_corner, classify_tangent_intersection, no_intersections, &
       remove_node
  public &
       Intersection, SegmentNode, IntersectionClassification_FIRST, &
       IntersectionClassification_SECOND, IntersectionClassification_OPPOSED, &
       IntersectionClassification_TANGENT_FIRST, &
       IntersectionClassification_TANGENT_SECOND, &
       IntersectionClassification_IGNORED_CORNER, SurfaceContained_NEITHER, &
       SurfaceContained_FIRST, SurfaceContained_SECOND, newton_refine, &
       locate_point, classify_intersection, add_st_vals, &
       surfaces_intersection_points, get_next, surfaces_intersect

  ! NOTE: This (for now) is not meant to be C-interoperable.
  type :: Intersection
     real(c_double) :: s = -1.0_dp
     real(c_double) :: t = -1.0_dp
     integer(c_int) :: index_first = -1
     integer(c_int) :: index_second = -1
     integer(c_int) :: interior_curve = -99  ! Hopefully an unused enum value
  end type Intersection

  ! For ``locate_point``; this is not meant to be C-interoperable.
  type :: LocateCandidate
     real(c_double) :: centroid_x = 1.0_dp  ! Actually triple.
     real(c_double) :: centroid_y = 1.0_dp  ! Actually triple.
     real(c_double) :: width = 1.0_dp
     real(c_double), allocatable :: nodes(:, :)
  end type LocateCandidate

  ! For ``surfaces_intersect`` and its helpers; this is not meant to be
  ! C-interoperable.
  ! NOTE: This type is not meant to be part of the interface for this
  !       module, but it is (for now) public, so that ``get_next()`` and
  !       ``to_front()`` can be tested.
  type :: SegmentNode
     real(c_double) :: edge_param = -1.0_dp
     integer(c_int) :: edge_index = -1
     integer(c_int) :: interior_curve = -99  ! Hopefully an unused enum value
  end type SegmentNode

  ! NOTE: These values are also defined in equivalent Python source.
  integer(c_int), parameter :: MAX_LOCATE_SUBDIVISIONS = 20
  real(c_double), parameter :: LOCATE_EPS = 0.5_dp**47
  ! Values of IntersectionClassification enum:
  integer(c_int), parameter :: IntersectionClassification_FIRST = 0
  integer(c_int), parameter :: IntersectionClassification_SECOND = 1
  integer(c_int), parameter :: IntersectionClassification_OPPOSED = 2
  integer(c_int), parameter :: IntersectionClassification_TANGENT_FIRST = 3
  integer(c_int), parameter :: IntersectionClassification_TANGENT_SECOND = 4
  integer(c_int), parameter :: IntersectionClassification_IGNORED_CORNER = 5
  ! Values of SurfaceContained enum:
  integer(c_int), parameter :: SurfaceContained_NEITHER = 0
  integer(c_int), parameter :: SurfaceContained_FIRST = 1
  integer(c_int), parameter :: SurfaceContained_SECOND = 2

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
    ! NOTE: This **assumes** but does not check that if the nodes are
    !       allocated, they are also the correct shape.
    if (.NOT. allocated(next_candidates(num_next_candidates - 3)%nodes)) then
       allocate(next_candidates(num_next_candidates - 3)%nodes(num_nodes, 2))
    end if
    if (.NOT. allocated(next_candidates(num_next_candidates - 2)%nodes)) then
       allocate(next_candidates(num_next_candidates - 2)%nodes(num_nodes, 2))
    end if
    if (.NOT. allocated(next_candidates(num_next_candidates - 1)%nodes)) then
       allocate(next_candidates(num_next_candidates - 1)%nodes(num_nodes, 2))
    end if
    if (.NOT. allocated(next_candidates(num_next_candidates)%nodes)) then
       allocate(next_candidates(num_next_candidates)%nodes(num_nodes, 2))
    end if

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

  subroutine allocate_candidates(num_candidates, next_candidates)
    integer(c_int), intent(in) :: num_candidates
    type(LocateCandidate), allocatable, intent(inout) :: next_candidates(:)

    if (allocated(next_candidates)) then
       if (size(next_candidates) < 4 * num_candidates) then
          deallocate(next_candidates)
          allocate(next_candidates(4 * num_candidates))
       end if
    else
       allocate(next_candidates(4 * num_candidates))
    end if

  end subroutine allocate_candidates

  subroutine update_candidates( &
       num_nodes, degree, point, num_candidates, candidates, &
       num_next_candidates, next_candidates)

    integer(c_int), intent(in) :: num_nodes, degree
    real(c_double), intent(in) :: point(1, 2)
    integer(c_int), intent(in) :: num_candidates
    type(LocateCandidate), intent(in) :: candidates(:)
    integer(c_int), intent(inout) :: num_next_candidates
    type(LocateCandidate), allocatable, intent(inout) :: next_candidates(:)
    ! Variables outside of signature.
    integer(c_int) :: cand_index
    logical(c_bool) :: predicate

    call allocate_candidates(num_candidates, next_candidates)

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

  end subroutine update_candidates

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
    integer(c_int) :: sub_index
    integer(c_int) :: num_candidates, num_next_candidates
    type(LocateCandidate), allocatable :: candidates_odd(:), candidates_even(:)
    real(c_double) :: s_approx, t_approx
    real(c_double) :: actual(1, 2)
    logical(c_bool) :: is_even

    point(1, :) = [x_val, y_val]
    ! Start out with the full curve.
    allocate(candidates_odd(1))  ! First iteration is odd.
    candidates_odd(1) = LocateCandidate(1.0_dp, 1.0_dp, 1.0_dp, nodes)
    num_candidates = 1
    s_val = LOCATE_MISS

    is_even = .TRUE.  ! At zero.
    do sub_index = 1, MAX_LOCATE_SUBDIVISIONS + 1
       is_even = .NOT. is_even  ! Switch parity.
       num_next_candidates = 0

       if (is_even) then
          ! Read from even, write to odd.
          call update_candidates( &
               num_nodes, degree, point, num_candidates, candidates_even, &
               num_next_candidates, candidates_odd)
       else
          ! Read from odd, write to even.
          call update_candidates( &
               num_nodes, degree, point, num_candidates, candidates_odd, &
               num_next_candidates, candidates_even)
       end if

       ! If there are no more candidates, we are done.
       if (num_next_candidates == 0) then
          return
       end if

       num_candidates = num_next_candidates

    end do

    ! Compute the s- and t-parameters as the mean of the centroid positions.
    ! We divide by `3n` rather than `n` since we have tracked thrice the
    ! centroid rather than the centroid itself (to avoid round-off until
    ! right now).
    if (is_even) then
       ! NOTE: We "exclude" this block from ``lcov`` because it can never
       !       be invoked while ``MAX_LOCATE_SUBDIVISIONS`` is even.
       ! LCOV_EXCL_START
       s_approx = ( &
            sum(candidates_odd(:num_candidates)%centroid_x) / &
            (3.0_dp * num_candidates))
       t_approx = ( &
            sum(candidates_odd(:num_candidates)%centroid_y) / &
            (3.0_dp * num_candidates))
       ! LCOV_EXCL_STOP
    else
       s_approx = ( &
            sum(candidates_even(:num_candidates)%centroid_x) / &
            (3.0_dp * num_candidates))
       t_approx = ( &
            sum(candidates_even(:num_candidates)%centroid_y) / &
            (3.0_dp * num_candidates))
    end if

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

  logical(c_bool) function ignored_edge_corner( &
       edge_tangent, corner_tangent, num_nodes, &
       previous_edge_nodes) result(predicate)

    real(c_double), intent(in) :: edge_tangent(1, 2)
    real(c_double), intent(in) :: corner_tangent(1, 2)
    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: previous_edge_nodes(num_nodes, 2)
    ! Variables outside of signature.
    real(c_double) :: cross_prod
    real(c_double) :: alt_corner_tangent(1, 2)

    call cross_product( &
         edge_tangent, corner_tangent, cross_prod)
    ! A negative cross product indicates that ``edge_tangent`` is
    ! "inside" / "to the left" of ``corner_tangent`` (due to right-hand rule).
    if (cross_prod > 0.0_dp) then
       predicate = .FALSE.
       return
    end if

    ! Do the same for the **other** tangent at the corner.
    call evaluate_hodograph( &
         1.0_dp, num_nodes, 2, &
         previous_edge_nodes, alt_corner_tangent)
    ! Change the direction of the "in" tangent so that it points "out".
    alt_corner_tangent = -alt_corner_tangent
    call cross_product( &
         edge_tangent, alt_corner_tangent, cross_prod)
    predicate = (cross_prod <= 0.0_dp)

  end function ignored_edge_corner

  logical(c_bool) function ignored_double_corner( &
       edges_first, edges_second, intersection_, &
       tangent_s, tangent_t) result(predicate)

    ! NOTE: This assumes that ``intersection_%index_(first|second)``
    !       are in [1, 2, 3] (marked as the edges of a surface).

    type(CurveData), intent(in) :: edges_first(3), edges_second(3)
    type(Intersection), intent(in) :: intersection_
    real(c_double), intent(in) :: tangent_s(1, 2), tangent_t(1, 2)
    ! Variables outside of signature.
    integer(c_int) :: index, num_nodes
    real(c_double) :: alt_tangent_s(1, 2), alt_tangent_t(1, 2)
    real(c_double) :: cross_prod1, cross_prod2, cross_prod3

    ! Compute the other edge for the ``s`` surface.
    index = 1 + modulo(intersection_%index_first - 2, 3)
    num_nodes = size(edges_first(index)%nodes, 1)
    call evaluate_hodograph( &
         1.0_dp, num_nodes, 2, &
         edges_first(index)%nodes, alt_tangent_s)

    ! First check if ``tangent_t`` is interior to the ``s`` surface.
    call cross_product( &
         tangent_s, tangent_t, cross_prod1)
    ! A positive cross product indicates that ``tangent_t`` is
    ! interior to ``tangent_s``. Similar for ``alt_tangent_s``.
    ! If ``tangent_t`` is interior to both, then the surfaces
    ! do more than just "kiss" at the corner, so the corner should
    ! not be ignored.
    if (cross_prod1 >= 0.0_dp) then
       ! Only compute ``cross_prod2`` if we need to.
       call cross_product( &
            alt_tangent_s, tangent_t, cross_prod2)
       if (cross_prod2 >= 0.0_dp) then
          predicate = .FALSE.
          return
       end if
    end if

    ! If ``tangent_t`` is not interior, we check the other ``t``
    ! edge that ends at the corner.
    index = 1 + modulo(intersection_%index_second - 2, 3)
    num_nodes = size(edges_second(index)%nodes, 1)
    call evaluate_hodograph( &
         1.0_dp, num_nodes, 2, &
         edges_second(index)%nodes, alt_tangent_t)
    ! Change the direction of the "in" tangent so that it points "out".
    alt_tangent_t = -alt_tangent_t

    call cross_product( &
         tangent_s, alt_tangent_t, cross_prod2)
    if (cross_prod2 >= 0.0_dp) then
       ! Only compute ``cross_prod3`` if we need to.
       call cross_product( &
            alt_tangent_s, alt_tangent_t, cross_prod3)
       if (cross_prod3 >= 0.0_dp) then
          predicate = .FALSE.
          return
       end if
    end if

    ! If neither of ``tangent_t`` or ``alt_tangent_t`` are interior
    ! to the ``s`` surface, one of two things is true. Either
    ! the two surfaces have no interior intersection (1) or the
    ! ``s`` surface is bounded by both edges of the ``t`` surface
    ! at the corner intersection (2). To detect (2), we only need
    ! check if ``tangent_s`` is interior to both ``tangent_t``
    ! and ``alt_tangent_t``. ``cross_prod1`` contains
    ! (tangent_s) x (tangent_t), so it's negative will tell if
    ! ``tangent_s`` is interior. Similarly, ``cross_prod3``
    ! contains (tangent_s) x (alt_tangent_t), but we also reversed
    ! the sign on ``alt_tangent_t`` so switching the sign back
    ! and reversing the arguments in the cross product cancel out.
    predicate = (cross_prod1 > 0.0_dp .OR. cross_prod2 < 0.0_dp)

  end function ignored_double_corner

  logical(c_bool) function ignored_corner( &
       edges_first, edges_second, intersection_, &
       tangent_s, tangent_t) result(predicate)

    ! NOTE: This assumes that ``intersection_%index_(first|second)``
    !       are in [1, 2, 3] (marked as the edges of a surface).

    type(CurveData), intent(in) :: edges_first(3), edges_second(3)
    type(Intersection), intent(in) :: intersection_
    real(c_double), intent(in) :: tangent_s(1, 2), tangent_t(1, 2)
    ! Variables outside of signature.
    integer(c_int) :: index, num_nodes

    if (intersection_%s == 0.0_dp) then
       if (intersection_%t == 0.0_dp) then
          ! Double corner.
          predicate = ignored_double_corner( &
               edges_first, edges_second, intersection_, &
               tangent_s, tangent_t)
       else
          ! s-only corner.
          index = 1 + modulo(intersection_%index_first - 2, 3)
          ! NOTE: This is a "trick" which requires `modulo` (not `mod`)
          !       since it will return values in {0, 1, 2}. The goal is
          !       to send [1, 2, 3] --> [3, 1, 2] (i.e. move indices to
          !       the left and wrap around).
          num_nodes = size(edges_first(index)%nodes, 1)
          predicate = ignored_edge_corner( &
               tangent_t, tangent_s, num_nodes, edges_first(index)%nodes)
       end if
    else if (intersection_%t == 0.0_dp) then
       ! t-only corner.
       index = 1 + modulo(intersection_%index_second - 2, 3)
       ! NOTE: This is a "trick" which requires `modulo` (not `mod`)
       !       since it will return values in {0, 1, 2}. The goal is
       !       to send [1, 2, 3] --> [3, 1, 2] (i.e. move indices to
       !       the left and wrap around).
       num_nodes = size(edges_second(index)%nodes, 1)
       predicate = ignored_edge_corner( &
            tangent_s, tangent_t, num_nodes, edges_second(index)%nodes)
    else
       predicate = .FALSE.
    end if

  end function ignored_corner

  subroutine classify_tangent_intersection( &
       edges_first, edges_second, intersection_, &
       tangent_s, tangent_t, enum_, status)

    ! NOTE: This **assumes**, but does not check that
    !       ``intersection_%index_(first|second)`` are valid indices within
    !       ``edges_(first|second)`` and that each of those ``CurveData``
    !       instances have already allocated ``%nodes``.

    ! Possible error states:
    ! * Status_SUCCESS       : On success.
    ! * Status_BAD_TANGENT   : If the curves move in an opposite direction
    !                          while defining overlapping arcs.
    ! * Status_SAME_CURVATURE: If the curves have identical curvature at the
    !                          point of tangency.

    type(CurveData), intent(in) :: edges_first(3), edges_second(3)
    type(Intersection), intent(in) :: intersection_
    real(c_double), intent(in) :: tangent_s(1, 2), tangent_t(1, 2)
    integer(c_int), intent(out) :: enum_
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    real(c_double) :: dot_prod
    integer(c_int) :: num_nodes
    real(c_double) :: curvature1, curvature2
    real(c_double) :: sign1, sign2
    real(c_double) :: delta_c

    status = Status_SUCCESS
    dot_prod = dot_product(tangent_s(1, :), tangent_t(1, :))
    ! NOTE: When computing curvatures we assume that we don't have lines
    !       here, because lines that are tangent at an intersection are
    !       parallel and we don't handle that case.
    num_nodes = size(edges_first(intersection_%index_first)%nodes, 1)
    call get_curvature( &
         num_nodes, 2, edges_first(intersection_%index_first)%nodes, &
         tangent_s, intersection_%s, curvature1)

    num_nodes = size(edges_second(intersection_%index_second)%nodes, 1)
    call get_curvature( &
         num_nodes, 2, edges_second(intersection_%index_second)%nodes, &
         tangent_t, intersection_%t, curvature2)

    if (dot_prod < 0.0_dp) then
       ! If the tangent vectors are pointing in the opposite direction,
       ! then the curves are facing opposite directions.
       sign1 = sign(1.0_dp, curvature1)
       sign2 = sign(1.0_dp, curvature2)
       if (sign1 == sign2) then
          ! If both curvatures are positive, since the curves are
          ! moving in opposite directions, the tangency isn't part of
          ! the surface intersection.
          if (sign1 == 1.0_dp) then
             enum_ = IntersectionClassification_OPPOSED
          else
             status = Status_BAD_TANGENT
          end if
       else
          delta_c = abs(curvature1) - abs(curvature2)
          if (delta_c == 0.0_dp) then
             status = Status_SAME_CURVATURE
          else
             sign2 = sign(1.0_dp, delta_c)
             if (sign1 == sign2) then
                enum_ = IntersectionClassification_OPPOSED
             else
                status = Status_BAD_TANGENT
             end if
          end if
       end if
    else
       if (curvature1 > curvature2) then
          enum_ = IntersectionClassification_TANGENT_FIRST
       else if (curvature1 < curvature2) then
          enum_ = IntersectionClassification_TANGENT_SECOND
       else
          status = Status_SAME_CURVATURE
       end if
    end if

  end subroutine classify_tangent_intersection

  subroutine classify_intersection( &
       edges_first, edges_second, intersection_, enum_, status)

    ! NOTE: This is **explicitly** not intended for C inter-op.
    ! NOTE: This subroutine is not part of the C ABI for this module,
    !       but it is (for now) public, so that it can be tested.
    ! NOTE: This **assumes**, but does not check that
    !       ``intersection_%index_(first|second)`` are valid indices within
    !       ``edges_(first|second)`` and that each of those ``CurveData``
    !       instances have already allocated ``%nodes``.

    ! Possible error states:
    ! * Status_SUCCESS       : On success.
    ! * Status_EDGE_END      : If either the s- or t-parameter of the
    !                          intersection is equal to ``1.0``, i.e. if the
    !                          intersection is at the end of an edge. Since
    !                          another edge of the surface will also contain
    !                          that point (at it's beginning), that edge
    !                          should be used.
    ! * Status_BAD_TANGENT   : Via ``classify_tangent_intersection()``.
    ! * Status_SAME_CURVATURE: Via ``classify_tangent_intersection()``.

    type(CurveData), intent(in) :: edges_first(3), edges_second(3)
    type(Intersection), intent(in) :: intersection_
    integer(c_int), intent(out) :: enum_
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    integer(c_int) :: num_nodes
    real(c_double) :: tangent_s(1, 2), tangent_t(1, 2)
    real(c_double) :: cross_prod

    status = Status_SUCCESS
    if (intersection_%s == 1.0_dp .OR. intersection_%t == 1.0_dp) then
       status = Status_EDGE_END
       return
    end if

    ! NOTE: We assume, but don't check that ``%nodes`` is allocated
    !       and that it has 2 columns.
    num_nodes = size(edges_first(intersection_%index_first)%nodes, 1)
    call evaluate_hodograph( &
         intersection_%s, num_nodes, 2, &
         edges_first(intersection_%index_first)%nodes, tangent_s)

    ! NOTE: We assume, but don't check that ``%nodes`` is allocated
    !       and that it has 2 columns.
    num_nodes = size(edges_second(intersection_%index_second)%nodes, 1)
    call evaluate_hodograph( &
         intersection_%t, num_nodes, 2, &
         edges_second(intersection_%index_second)%nodes, tangent_t)

    if (ignored_corner( &
         edges_first, edges_second, intersection_, &
         tangent_s, tangent_t)) then
       enum_ = IntersectionClassification_IGNORED_CORNER
       return
    end if

    ! Take the cross product of tangent vectors to determine which one
    ! is more "inside" / "to the left".
    call cross_product( &
         tangent_s, tangent_t, cross_prod)
    if (cross_prod < 0.0_dp) then
       enum_ = IntersectionClassification_FIRST
    else if (cross_prod > 0.0_dp) then
       enum_ = IntersectionClassification_SECOND
    else
       call classify_tangent_intersection( &
            edges_first, edges_second, intersection_, &
            tangent_s, tangent_t, enum_, status)
    end if

  end subroutine classify_intersection

  subroutine add_st_vals( &
       edges_first, edges_second, num_st_vals, st_vals, &
       index_first, index_second, intersections, num_intersections, &
       all_types, status)

    ! NOTE: This is **explicitly** not intended for C inter-op.
    ! NOTE: This subroutine is not part of the C ABI for this module,
    !       but it is (for now) public, so that it can be tested.
    ! NOTE: This assumes but does not check that ``num_st_vals > 0``.
    ! NOTE: This assumes but does not check that each of ``index_first``
    !       and ``index_second`` is in {1, 2, 3}.

    ! Possible error states:
    ! * Status_SUCCESS       : On success.
    ! * Status_EDGE_END      : Via ``classify_intersection()``.
    ! * Status_BAD_TANGENT   : Via ``classify_intersection()``.
    ! * Status_SAME_CURVATURE: Via ``classify_intersection()``.

    type(CurveData), intent(in) :: edges_first(3), edges_second(3)
    integer(c_int), intent(in) :: num_st_vals
    real(c_double), intent(in) :: st_vals(2, num_st_vals)
    integer(c_int), intent(in) :: index_first, index_second
    type(Intersection), allocatable, intent(inout) :: intersections(:)
    integer(c_int), intent(inout) :: num_intersections
    integer(c_int), intent(inout) :: all_types
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    integer(c_int) :: curr_size, i, intersection_index
    type(Intersection), allocatable :: intersections_swap(:)
    integer(c_int) :: enum_

    status = Status_SUCCESS

    intersection_index = num_intersections
    ! NOTE: Intersections at the end of an edge will be skipped, so this
    !       may over-estimate the number of intersections.
    num_intersections = num_intersections + num_st_vals
    if (allocated(intersections)) then
       curr_size = size(intersections)
       if (curr_size < num_intersections) then
          allocate(intersections_swap(num_intersections))
          intersections_swap(:curr_size) = intersections(:curr_size)
          call move_alloc(intersections_swap, intersections)
       end if
    else
       allocate(intersections(num_intersections))
    end if

    do i = 1, num_st_vals
       if (st_vals(1, i) == 1.0_dp .OR. st_vals(2, i) == 1.0_dp) then
          ! If the intersection is at the end of an edge, ignore it.
          ! This intersection **could** be saved to verify that it matches an
          ! equivalent intersection from the beginning of the previous edge.
          cycle
       end if
       intersection_index = intersection_index + 1
       intersections(intersection_index)%s = st_vals(1, i)
       intersections(intersection_index)%t = st_vals(2, i)
       intersections(intersection_index)%index_first = index_first
       intersections(intersection_index)%index_second = index_second
       call classify_intersection( &
            edges_first, edges_second, &
            intersections(intersection_index), enum_, status)
       if (status /= Status_SUCCESS) then
          return
       end if
       ! NOTE: This assumes, but does not check, that ``enum_`` is in
       !       {0, 1, 2, 3, 4, 5} (should be an enum value from
       !       ``IntersectionClassification``) which limits the value of
       !       ``all_types`` to [0, 63] (inclusive).
       all_types = ior(all_types, 2**enum_)
       if ( &
            enum_ == IntersectionClassification_FIRST .OR. &
            enum_ == IntersectionClassification_SECOND) then
          ! "Keep" the intersection if it is ``FIRST`` or ``SECOND``.
          intersections(intersection_index)%interior_curve = enum_
       else
          ! "Discard" the intersection (rather, allow it to be over-written).
          intersection_index = intersection_index - 1
       end if
    end do

    ! Actually set ``num_intersections`` based on the number of **accepted**
    ! ``s-t`` pairs.
    num_intersections = intersection_index

  end subroutine add_st_vals

  subroutine surfaces_intersection_points( &
       num_nodes1, nodes1, degree1, &
       num_nodes2, nodes2, degree2, &
       intersections, num_intersections, all_types, status)

    ! NOTE: This is **explicitly** not intended for C inter-op.
    ! NOTE: This (and ``add_st_vals``) ignores duplicate nodes (caused when
    !       an intersection happens at a corner of the surface, which will
    !       be at the end of one edge and the start of another). If desired,
    !       a "verify" option could be added (as is done in Python) to make
    !       sure that any "duplicate" intersections (i.e. ones that occur at
    !       the end of an edge) do match with the corresponding intersection
    !       at the beginning of an edge (in the case of a double corner, this
    !       means **three** duplicate intersections).
    ! NOTE: The returned ``all_types`` uses the first 6 bits to identify
    !       which classified states are among the intersections. (The
    !       possible classifications must be in {0, 1, 2, 3, 4, 5}).

    ! Possible error states:
    ! * Status_SUCCESS       : On success.
    ! * Status_PARALLEL      : Via ``curve_intersection.all_intersections()``.
    ! * Status_WIGGLE_FAIL   : Via ``curve_intersection.all_intersections()``.
    ! * Status_NO_CONVERGE   : Via ``curve_intersection.all_intersections()``.
    ! * (N >= 64)            : Via ``curve_intersection.all_intersections()``.
    ! * Status_EDGE_END      : Via ``add_st_vals()``.
    ! * Status_BAD_TANGENT   : Via ``add_st_vals()``.
    ! * Status_SAME_CURVATURE: Via ``add_st_vals()``.

    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(num_nodes1, 2)
    integer(c_int), intent(in) :: degree1
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(num_nodes2, 2)
    integer(c_int), intent(in) :: degree2
    type(Intersection), allocatable, intent(inout) :: intersections(:)
    integer(c_int), intent(out) :: num_intersections
    integer(c_int), intent(out) :: all_types
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    type(CurveData) :: edges_first(3), edges_second(3)
    integer(c_int) :: index1, index2
    real(c_double), allocatable :: st_vals(:, :)
    integer(c_int) :: num_st_vals

    ! Compute the edge nodes for the first surface.
    allocate(edges_first(1)%nodes(degree1 + 1, 2))
    allocate(edges_first(2)%nodes(degree1 + 1, 2))
    allocate(edges_first(3)%nodes(degree1 + 1, 2))
    call compute_edge_nodes( &
         num_nodes1, 2, nodes1, degree1, &
         edges_first(1)%nodes, edges_first(2)%nodes, edges_first(3)%nodes)

    ! Compute the edge nodes for the second surface.
    allocate(edges_second(1)%nodes(degree2 + 1, 2))
    allocate(edges_second(2)%nodes(degree2 + 1, 2))
    allocate(edges_second(3)%nodes(degree2 + 1, 2))
    call compute_edge_nodes( &
         num_nodes2, 2, nodes2, degree2, &
         edges_second(1)%nodes, edges_second(2)%nodes, edges_second(3)%nodes)

    num_intersections = 0
    all_types = 0
    do index1 = 1, 3
       do index2 = 1, 3
          call all_intersections( &
               degree1 + 1, edges_first(index1)%nodes, &
               degree2 + 1, edges_second(index2)%nodes, &
               st_vals, num_st_vals, status)
          if (status == Status_SUCCESS) then
             if (num_st_vals > 0) then
                call add_st_vals( &
                     edges_first, edges_second, num_st_vals, st_vals, &
                     index1, index2, intersections, num_intersections, &
                     all_types, status)
                if (status /= Status_SUCCESS) then
                   return
                end if
             end if
          else
             return
          end if
       end do
    end do

  end subroutine surfaces_intersection_points

  subroutine no_intersections( &
       num_nodes1, nodes1, degree1, &
       num_nodes2, nodes2, degree2, contained)

    ! NOTE: This is a helper for ``surfaces_intersect()``.

    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(num_nodes1, 2)
    integer(c_int), intent(in) :: degree1
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(num_nodes2, 2)
    integer(c_int), intent(in) :: degree2
    integer(c_int), intent(out) :: contained
    ! Variables outside of signature.
    real(c_double) :: s_val, t_val

    ! If the first corner of ``surface1`` is contained in ``surface2``,
    ! then the whole surface must be since there are no intersections.
    call locate_point( &
         num_nodes2, nodes2, degree2, &
         nodes1(1, 1), nodes1(1, 2), s_val, t_val)
    if (s_val /= LOCATE_MISS) then
       contained = SurfaceContained_FIRST
       return
    end if

    ! If the first corner of ``surface2`` is contained in ``surface1``,
    ! then the whole surface must be since there are no intersections.
    call locate_point( &
         num_nodes1, nodes1, degree1, &
         nodes2(1, 1), nodes2(1, 2), s_val, t_val)
    if (s_val /= LOCATE_MISS) then
       contained = SurfaceContained_SECOND
       return
    end if

    ! If neither corner is contained, then there is no intersection.
    contained = SurfaceContained_NEITHER

  end subroutine no_intersections

  subroutine remove_node(node, values, remaining)

    ! Removes a value from a list, if it is contained there.
    ! For example, if ``node == 4``, ``values == [1, 2, 4, 7, ...]`` and
    ! ``remainining == 4``, then after this subroutine will update to
    ! ``values == [1, 2, 7, ...]`` and ``remainining == 3``.
    ! NOTE: This assumes the values in ``values`` are in ascending order
    !       and unique.
    ! NOTE: This assumes that ``remaining`` does not exceed ``size(values)``

    integer(c_int), intent(in) :: node
    integer(c_int), intent(inout) :: values(:)
    integer(c_int), intent(inout) :: remaining
    ! Variables outside of signature.
    integer(c_int) :: i

    ! NOTE: Since we assume the values are in ascending order and unique, a
    !       binary search could speed this up. In practice, ``values``
    !       will likely be too small for this to matter.
    do i = 1, remaining
       if (values(i) == node) then
          ! Remove ``values(i)`` by shifting down all remaining values.
          values(i:remaining - 1) = values(i + 1:remaining)
          remaining = remaining - 1
          return
       end if
    end do

  end subroutine remove_node

  subroutine get_next( &
       num_intersections, intersections, unused, remaining, &
       start, curr_node, next_node, at_start)

    ! NOTE: This subroutine is not meant to be part of the interface for this
    !       module, but it is (for now) public, so that it can be tested.

    integer(c_int), intent(in) :: num_intersections
    type(Intersection), intent(in) :: intersections(num_intersections)
    ! ``unused`` contains the indices of intersections that have not yet been
    ! used as a ``SegmentNode``.
    integer(c_int), intent(inout) :: unused(:)
    integer(c_int), intent(inout) :: remaining
    integer(c_int), intent(in) :: start
    type(SegmentNode), intent(in) :: curr_node
    type(SegmentNode), intent(out) :: next_node
    logical(c_bool), intent(out) :: at_start
    ! Variables outside of signature.
    integer(c_int) :: i, local_index, intersection_index

    intersection_index = -1
    at_start = .FALSE.
    if (curr_node%edge_index <= 3) then
       ! NOTE: This assumes but does not check that ``curr_node%edge_index``
       !       is in {1, 2, 3}.
       local_index = curr_node%edge_index
       do i = 1, num_intersections
          if ( &
               intersections(i)%index_first == local_index .AND. &
               intersections(i)%s > curr_node%edge_param) then
             if (intersection_index == -1) then
                intersection_index = i
                next_node%edge_index = curr_node%edge_index
                next_node%edge_param = intersections(i)%s
                next_node%interior_curve = intersections(i)%interior_curve
             else
                if (intersections(i)%s < next_node%edge_param) then
                   ! Assumes ``edge_index`` is already set.
                   intersection_index = i
                   next_node%edge_param = intersections(i)%s
                   next_node%interior_curve = intersections(i)%interior_curve
                end if
             end if
          end if
       end do

       ! If there is no other intersection on the edge, just return
       ! the segment end.
       if (intersection_index == -1) then
          next_node%edge_index = curr_node%edge_index
          next_node%edge_param = 1.0_dp
          next_node%interior_curve = IntersectionClassification_FIRST
       else
          at_start = (intersection_index == start)
          ! Remove the index from the set of ``unused`` intersections, if
          ! it is contained there.
          call remove_node(intersection_index, unused, remaining)
       end if
    else
       ! NOTE: This assumes but does not check that ``curr_node%edge_index``
       !       is in {4, 5, 6}.
       local_index = curr_node%edge_index - 3
       do i = 1, num_intersections
          if ( &
               intersections(i)%index_second == local_index .AND. &
               intersections(i)%t > curr_node%edge_param) then
             if (intersection_index == -1) then
                intersection_index = i
                next_node%edge_index = curr_node%edge_index
                next_node%edge_param = intersections(i)%t
                next_node%interior_curve = intersections(i)%interior_curve
             else
                if (intersections(i)%t < next_node%edge_param) then
                   ! Assumes ``edge_index`` is already set.
                   intersection_index = i
                   next_node%edge_param = intersections(i)%t
                   next_node%interior_curve = intersections(i)%interior_curve
                end if
             end if
          end if
       end do

       ! If there is no other intersection on the edge, just return
       ! the segment end.
       if (intersection_index == -1) then
          next_node%edge_index = curr_node%edge_index
          next_node%edge_param = 1.0_dp
          next_node%interior_curve = IntersectionClassification_SECOND
       else
          at_start = (intersection_index == start)
          ! Remove the index from the set of ``unused`` intersections, if
          ! it is contained there.
          call remove_node(intersection_index, unused, remaining)
       end if
    end if

  end subroutine get_next

  subroutine to_front( &
       num_intersections, intersections, start, curr_node, next_node, at_start)

    ! NOTE: This subroutine is not meant to be part of the interface for this
    !       module, but it is (for now) public, so that it can be tested.

    integer(c_int), intent(in) :: num_intersections
    type(Intersection), intent(in) :: intersections(num_intersections)
    integer(c_int), intent(in) :: start
    type(SegmentNode), intent(in) :: curr_node
    type(SegmentNode), intent(out) :: next_node
    logical(c_bool), intent(out) :: at_start
    ! Variables outside of signature.
    integer(c_int) :: i, local_index

    at_start = .FALSE.
    if (curr_node%edge_param == 1.0_dp) then
       next_node%edge_param = 0.0_dp
       if (curr_node%edge_index == 3 .OR. curr_node%edge_index == 6) then
          next_node%edge_index = curr_node%edge_index - 2
       else
          next_node%edge_index = curr_node%edge_index + 1
       end if
       next_node%interior_curve = curr_node%interior_curve
    else
       ! The node doesn't need to be moved to the front of the edge.
       next_node = curr_node
       return
    end if

    ! NOTE: If we've reached this point, we know that ``edge_param == 0.0``.
    if (next_node%edge_index <= 3) then
       ! NOTE: This assumes but does not check that ``next_node%edge_index``
       !       is in {1, 2, 3}.
       local_index = next_node%edge_index
       do i = 1, num_intersections
          if ( &
               intersections(i)%s == 0.0_dp .AND. &
               intersections(i)%index_first == local_index) then
             ! Assumes ``edge_param`` and ``edge_index`` are already set.
             next_node%interior_curve = intersections(i)%interior_curve
             at_start = (i == start)
             return
          end if
       end do
    else
       ! NOTE: This assumes but does not check that ``next_node%edge_index``
       !       is in {4, 5, 6}.
       local_index = next_node%edge_index - 3
       do i = 1, num_intersections
          if ( &
               intersections(i)%t == 0.0_dp .AND. &
               intersections(i)%index_second == local_index) then
             ! Assumes ``edge_param`` and ``edge_index`` are already set.
             next_node%interior_curve = intersections(i)%interior_curve
             at_start = (i == start)
             return
          end if
       end do
    end if

  end subroutine to_front

  subroutine surfaces_intersect( &
       num_nodes1, nodes1, degree1, &
       num_nodes2, nodes2, degree2, &
       contained, status)

    ! NOTE: ``contained`` will be a ``SurfaceContained`` enum.

    ! Possible error states:
    ! * Status_SUCCESS       : On success.
    ! * Status_PARALLEL      : Via ``surfaces_intersection_points()``.
    ! * Status_WIGGLE_FAIL   : Via ``surfaces_intersection_points()``.
    ! * Status_NO_CONVERGE   : Via ``surfaces_intersection_points()``.
    ! * (N >= 64)            : Via ``surfaces_intersection_points()``.
    ! * Status_EDGE_END      : Via ``surfaces_intersection_points()``.
    ! * Status_BAD_TANGENT   : Via ``surfaces_intersection_points()``.
    ! * Status_SAME_CURVATURE: Via ``surfaces_intersection_points()``.
    ! * Status_UNKNOWN       : If all of the intersections are classified
    !                          as ``OPPOSED / IGNORED_CORNER / TANGENT_*``
    !                          but not uniquely one type. (This should
    !                          never occur).

    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(num_nodes1, 2)
    integer(c_int), intent(in) :: degree1
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(num_nodes2, 2)
    integer(c_int), intent(in) :: degree2
    integer(c_int), intent(out) :: contained
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    integer(c_int) :: bbox_int
    type(Intersection), allocatable :: intersections(:)
    integer(c_int) :: num_intersections, all_types

    ! If the bounded boxes do not intersect, the surfaces cannot.
    call bbox_intersect( &
         num_nodes1, nodes1, num_nodes2, nodes2, bbox_int)
    if (bbox_int /= BoxIntersectionType_INTERSECTION) then
       contained = SurfaceContained_NEITHER
       status = Status_SUCCESS
       return
    end if

    call surfaces_intersection_points( &
       num_nodes1, nodes1, degree1, &
       num_nodes2, nodes2, degree2, &
       intersections, num_intersections, all_types, status)
    if (status /= Status_SUCCESS) then
       return
    end if

    if (num_intersections == 0) then
       if (all_types == 0) then
          call no_intersections( &
               num_nodes1, nodes1, degree1, &
               num_nodes2, nodes2, degree2, contained)
          return
       else if (all_types == 2**IntersectionClassification_OPPOSED) then
          contained = SurfaceContained_NEITHER
          return
       else if (all_types == 2**IntersectionClassification_IGNORED_CORNER) then
          contained = SurfaceContained_NEITHER
          return
       else if (all_types == 2**IntersectionClassification_TANGENT_FIRST) then
          contained = SurfaceContained_FIRST
          return
       else if (all_types == 2**IntersectionClassification_TANGENT_SECOND) then
          contained = SurfaceContained_SECOND
          return
       else
          ! NOTE: We "exclude" this block from ``lcov`` because it **should**
          !       never occur. In the case that **none** of the intersections
          !       are ``FIRST`` or ``SECOND``, then they all **should** fall
          !       into the exact same category (which would be one of the
          !       four states above).
          ! LCOV_EXCL_START
          status = Status_UNKNOWN
          return
          ! LCOV_EXCL_STOP
       end if
    end if

    ! The remaining parts are not yet implemented.
    status = -1234

  end subroutine surfaces_intersect

end module surface_intersection
