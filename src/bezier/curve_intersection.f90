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

module curve_intersection

  use, intrinsic :: iso_c_binding, only: &
       c_double, c_int, c_bool, c_f_pointer, c_loc
  use types, only: dp
  use status, only: &
       Status_SUCCESS, Status_PARALLEL, Status_NO_CONVERGE, &
       Status_INSUFFICIENT_SPACE
  use helpers, only: &
       VECTOR_CLOSE_EPS, cross_product, bbox, wiggle_interval, &
       vector_close, in_interval, ulps_away, convex_hull, polygon_collide
  use curve, only: &
       CurveData, LOCATE_MISS, LOCATE_INVALID, evaluate_multi, &
       specialize_curve, evaluate_hodograph, locate_point, elevate_nodes, &
       subdivide_curve
  implicit none
  private &
       MAX_INTERSECT_SUBDIVISIONS, MIN_INTERVAL_WIDTH, MAX_CANDIDATES, &
       SIMILAR_ULPS, CANDIDATES_ODD, CANDIDATES_EVEN, make_candidates, &
       prune_candidates, elevate_helper
  public &
       BoxIntersectionType_INTERSECTION, BoxIntersectionType_TANGENT, &
       BoxIntersectionType_DISJOINT, Subdivide_FIRST, Subdivide_SECOND, &
       Subdivide_BOTH, Subdivide_NEITHER, LINEARIZATION_THRESHOLD, &
       INTERSECTIONS_WORKSPACE, linearization_error, segment_intersection, &
       newton_refine_intersect, bbox_intersect, parallel_different, &
       from_linearized, bbox_line_intersect, add_intersection, &
       add_from_linearized, endpoint_check, tangent_bbox_intersection, &
       add_candidates, intersect_one_round, make_same_degree, &
       add_coincident_parameters, all_intersections, all_intersections_abi, &
       set_max_candidates, get_max_candidates, set_similar_ulps, &
       get_similar_ulps, free_curve_intersections_workspace

  ! Values of BoxIntersectionType enum:
  integer(c_int), parameter :: BoxIntersectionType_INTERSECTION = 0
  integer(c_int), parameter :: BoxIntersectionType_TANGENT = 1
  integer(c_int), parameter :: BoxIntersectionType_DISJOINT = 2
  ! Values of Subdivide enum:
  integer(c_int), parameter :: Subdivide_FIRST = 0
  integer(c_int), parameter :: Subdivide_SECOND = 1
  integer(c_int), parameter :: Subdivide_BOTH = 2
  integer(c_int), parameter :: Subdivide_NEITHER = -1
  ! Set the threshold for linearization error at half the bits available.
  real(c_double), parameter :: LINEARIZATION_THRESHOLD = 0.5_dp**26
  integer(c_int), parameter :: MAX_INTERSECT_SUBDIVISIONS = 20
  real(c_double), parameter :: MIN_INTERVAL_WIDTH = 0.5_dp**40
  ! Run-time parameters that can be modified. If multiple threads are used,
  ! these **should** be thread-local (though it's expected that callers will
  ! update these values **before** beginning computation).
  integer(c_int) :: MAX_CANDIDATES = 64
  integer(c_int) :: SIMILAR_ULPS = 1
  ! Long-lived workspaces for ``all_intersections()`` and
  ! ``all_intersections_abi()``. If multiple threads are used, each of these
  ! **should** be thread-local.
  type(CurveData), allocatable :: CANDIDATES_ODD(:, :)
  type(CurveData), allocatable :: CANDIDATES_EVEN(:, :)
  real(c_double), allocatable :: INTERSECTIONS_WORKSPACE(:, :)

contains

  subroutine linearization_error( &
       num_nodes, dimension_, nodes, error)

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(dimension_, num_nodes)
    real(c_double), intent(out) :: error
    ! Variables outside of signature.
    real(c_double) :: second_deriv(dimension_, num_nodes - 2)
    real(c_double) :: worst_case(dimension_)

    ! A line has no linearization error.
    if (num_nodes == 2) then
       error = 0.0_dp
       return
    end if

    second_deriv = ( &
         nodes(:, :num_nodes - 2) - &
         2.0_dp * nodes(:, 2:num_nodes - 1) + &
         nodes(:, 3:))
    worst_case = maxval(abs(second_deriv), 2)
    error = 0.125_dp * (num_nodes - 1) * (num_nodes - 2) * norm2(worst_case)

  end subroutine linearization_error

  subroutine segment_intersection( &
       start0, end0, start1, end1, s, t, success)

    real(c_double), intent(in) :: start0(2)
    real(c_double), intent(in) :: end0(2)
    real(c_double), intent(in) :: start1(2)
    real(c_double), intent(in) :: end1(2)
    real(c_double), intent(out) :: s, t
    logical(c_bool), intent(out) :: success
    ! Variables outside of signature.
    real(c_double) :: delta0(2)
    real(c_double) :: delta1(2)
    real(c_double) :: start_delta(2)
    real(c_double) :: cross_d0_d1
    real(c_double) :: other_cross

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

  subroutine newton_refine_intersect( &
       s, num_nodes1, nodes1, t, num_nodes2, nodes2, new_s, new_t) &
       bind(c, name='newton_refine_curve_intersect')

    real(c_double), intent(in) :: s
    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(2, num_nodes1)
    real(c_double), intent(in) :: t
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(2, num_nodes2)
    real(c_double), intent(out) :: new_s, new_t
    ! Variables outside of signature.
    real(c_double) :: param(1)
    real(c_double) :: func_val(2, 1)
    real(c_double) :: workspace(2, 1)
    real(c_double) :: jac_mat(2, 2)
    real(c_double) :: determinant, delta_s, delta_t

    param = t
    call evaluate_multi( &
         num_nodes2, 2, nodes2, 1, param, func_val)
    param = s
    call evaluate_multi( &
         num_nodes1, 2, nodes1, 1, param, workspace)
    func_val = func_val - workspace

    if (all(func_val == 0.0_dp)) then
       new_s = s
       new_t = t
       return
    end if

    call evaluate_hodograph(s, num_nodes1, 2, nodes1, jac_mat(1:1, :))
    ! NOTE: We actually want the negative, since we want -B2'(t), but
    !       since we manually solve the system, it's just algebra
    !       to figure out how to use the negative values.
    call evaluate_hodograph(t, num_nodes2, 2, nodes2, jac_mat(2:2, :))

    determinant = ( &
         jac_mat(1, 1) * jac_mat(2, 2) - jac_mat(1, 2) * jac_mat(2, 1))

    ! NOTE: We manually invert the 2x2 system ([ds, dt] J)^T = f^T.
    delta_s = ( &
         (jac_mat(2, 2) * func_val(1, 1) - &
         jac_mat(2, 1) * func_val(2, 1)) / determinant)
    new_s = s + delta_s

    delta_t = ( &
         (jac_mat(1, 2) * func_val(1, 1) - &
         jac_mat(1, 1) * func_val(2, 1)) / determinant)
    new_t = t + delta_t

  end subroutine newton_refine_intersect

  subroutine bbox_intersect( &
       num_nodes1, nodes1, num_nodes2, nodes2, enum_) &
       bind(c, name='bbox_intersect')

    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(2, num_nodes1)
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(2, num_nodes2)
    integer(c_int), intent(out) :: enum_
    ! Variables outside of signature.
    real(c_double) :: left1, right1, bottom1, top1
    real(c_double) :: left2, right2, bottom2, top2

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

  subroutine parallel_different( &
       start0, end0, start1, end1, result_)

    real(c_double), intent(in) :: start0(2)
    real(c_double), intent(in) :: end0(2)
    real(c_double), intent(in) :: start1(2)
    real(c_double), intent(in) :: end1(2)
    logical(c_bool), intent(out) :: result_
    ! Variables outside of signature.
    real(c_double) :: delta0(2)
    real(c_double) :: val1, val2  ! Workspace

    delta0 = end0 - start0
    call cross_product(start0, delta0, val1)  ! line0_const
    call cross_product(start1, delta0, val2)  ! start1_against

    if (val1 /= val2) then
       result_ = .TRUE.
       return
    end if

    val1 = dot_product(delta0, delta0)  ! norm0_sq
    ! val2 == start_numer
    val2 = dot_product(start1 - start0, delta0)
    !      0 <= start_numer / norm0_sq <= 1
    ! <==> 0 <= start_numer            <= norm0_sq
    if (0.0_dp <= val2 .AND. val2 <= val1) then
       result_ = .FALSE.
       return
    end if

    val2 = dot_product(end1 - start0, delta0)  ! end_numer
    !      0 <= end_numer / norm0_sq <= 1
    ! <==> 0 <= end_numer            <= norm0_sq
    if (0.0_dp <= val2 .AND. val2 <= val1) then
       result_ = .FALSE.
       return
    end if

    ! We know neither the start or end parameters are in [0, 1], but
    ! they may contain [0, 1] between them, so we make sure that 0
    ! isn't between them.
    result_ = (0.0_dp < min(val1, val2) .OR. max(val1, val2) < 0.0_dp)

  end subroutine parallel_different

  subroutine from_linearized( &
       error1, curve1, num_nodes1, root_nodes1, &
       error2, curve2, num_nodes2, root_nodes2, &
       refined_s, refined_t, does_intersect, status)

    ! Possible error states:
    ! * Status_SUCCESS : On success.
    ! * Status_PARALLEL: If ``segment_intersection()`` fails (which means
    !                    the linearized segments are parallel). This
    !                    can still be avoided if the "root" curves are
    !                    also (parallel) lines that don't overlap or if
    !                    the "root" curves have disjoint bounding boxes.

    real(c_double), intent(in) :: error1
    type(CurveData), intent(in) :: curve1
    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: root_nodes1(2, num_nodes1)
    real(c_double), intent(in) :: error2
    type(CurveData), intent(in) :: curve2
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: root_nodes2(2, num_nodes2)
    real(c_double), intent(out) :: refined_s, refined_t
    logical(c_bool), intent(out) :: does_intersect
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    real(c_double) :: s, t
    integer(c_int) :: enum_
    logical(c_bool) :: success

    status = Status_SUCCESS
    does_intersect = .FALSE.  ! Default value.
    call segment_intersection( &
         curve1%nodes(:, 1), curve1%nodes(:, num_nodes1), &
         curve2%nodes(:, 1), curve2%nodes(:, num_nodes2), &
         s, t, success)

    if (success) then
       ! Special case for lines, allow no leeway on almost intersections.
       if (error1 == 0.0_dp .AND. (s < 0.0_dp .OR. 1.0_dp < s)) then
          return
       end if

       if (error2 == 0.0_dp .AND. (t < 0.0_dp .OR. 1.0_dp < t)) then
          return
       end if

       if (s < -(0.5_dp**16) .OR. 1.0_dp + 0.5_dp**16 < s) then
          return
       end if

       if (t < -(0.5_dp**16) .OR. 1.0_dp + 0.5_dp**16 < t) then
          return
       end if
    else
       ! Handle special case where the curves are actually lines.
       if (error1 == 0.0_dp .AND. error2 == 0.0_dp) then
          call parallel_different( &
               curve1%nodes(:, 1), curve1%nodes(:, num_nodes1), &
               curve2%nodes(:, 1), curve2%nodes(:, num_nodes2), &
               success)
          if (success) then
             return
          end if
       else
          call bbox_intersect( &
               num_nodes1, root_nodes1, num_nodes2, root_nodes2, enum_)
          if (enum_ == BoxIntersectionType_DISJOINT) then
             return
          end if
       end if

       status = Status_PARALLEL
       return
    end if

    ! Now, promote ``s`` and ``t`` onto the original curves.
    s = (1.0_dp - s) * curve1%start + s * curve1%end_  ! orig_s
    t = (1.0_dp - t) * curve2%start + t * curve2%end_  ! orig_t
    ! Perform one step of Newton iteration to refine the computed
    ! values of s and t.
    call newton_refine_intersect( &
         s, num_nodes1, root_nodes1, t, &
         num_nodes2, root_nodes2, refined_s, refined_t)

    call wiggle_interval(refined_s, s, success)
    if (.NOT. success) then
       return
    end if

    call wiggle_interval(refined_t, t, success)
    if (.NOT. success) then
       return
    end if

    does_intersect = .TRUE.
    refined_s = s
    refined_t = t

  end subroutine from_linearized

  subroutine bbox_line_intersect( &
       num_nodes, nodes, line_start, line_end, enum_)

    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: nodes(2, num_nodes)
    real(c_double), intent(in) :: line_start(2)
    real(c_double), intent(in) :: line_end(2)
    integer(c_int), intent(out) :: enum_
    ! Variables outside of signature.
    real(c_double) :: left, right, bottom, top
    real(c_double) :: segment_start(2)
    real(c_double) :: segment_end(2)
    real(c_double) :: s_curr, t_curr
    logical(c_bool) :: success

    call bbox(num_nodes, nodes, left, right, bottom, top)

    ! Check if line start is inside the bounding box.
    if ( &
         in_interval(line_start(1), left, right) .AND. &
         in_interval(line_start(2), bottom, top)) then
       enum_ = BoxIntersectionType_INTERSECTION
       return
    end if

    ! Check if line end is inside the bounding box.
    if ( &
         in_interval(line_end(1), left, right) .AND. &
         in_interval(line_end(2), bottom, top)) then
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
    segment_start(1) = left
    segment_start(2) = bottom
    segment_end(1) = right
    segment_end(2) = bottom
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
    segment_start(1) = right
    segment_start(2) = bottom
    segment_end(1) = right
    segment_end(2) = top
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
    segment_start(1) = right
    segment_start(2) = top
    segment_end(1) = left
    segment_end(2) = top
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

  subroutine add_intersection( &
       s, t, num_intersections, intersections)

    ! Adds an intersection to list of ``intersections``.

    real(c_double), intent(in) :: s
    real(c_double), intent(in) :: t
    integer(c_int), intent(inout) :: num_intersections
    real(c_double), allocatable, intent(inout) :: intersections(:, :)
    ! Variables outside of signature.
    integer(c_int) :: curr_size, index_
    real(c_double), allocatable :: intersections_swap(:, :)

    ! First, check if the intersection is a duplicate (up to precision,
    ! determined by ``ulps_away``).
    do index_ = 1, num_intersections
       if ( &
            ulps_away(intersections(1, index_), s, &
            SIMILAR_ULPS, VECTOR_CLOSE_EPS) .AND. &
            ulps_away(intersections(2, index_), t, &
            SIMILAR_ULPS, VECTOR_CLOSE_EPS)) then
          return
       end if
    end do

    ! Update the number of intersections.
    num_intersections = num_intersections + 1

    if (allocated(intersections)) then
       curr_size = size(intersections, 2)
       if (curr_size < num_intersections) then
          allocate(intersections_swap(2, num_intersections))
          intersections_swap(:, :curr_size) = intersections(:, :curr_size)
          call move_alloc(intersections_swap, intersections)
       end if
    else
       allocate(intersections(2, num_intersections))
    end if

    intersections(1, num_intersections) = s
    intersections(2, num_intersections) = t

  end subroutine add_intersection

  subroutine add_from_linearized( &
       first, root_nodes1, linearization_error1, &
       second, root_nodes2, linearization_error2, &
       num_intersections, intersections, status)

    ! Adds an intersection from two linearizations.
    !
    ! NOTE: This is **explicitly** not intended for C inter-op.

    ! Possible error states:
    ! * Status_SUCCESS : On success.
    ! * Status_PARALLEL: Via ``from_linearized()``.

    type(CurveData), intent(in) :: first
    real(c_double), intent(in) :: root_nodes1(:, :)
    real(c_double), intent(in) :: linearization_error1
    type(CurveData), intent(in) :: second
    real(c_double), intent(in) :: root_nodes2(:, :)
    real(c_double), intent(in) :: linearization_error2
    integer(c_int), intent(inout) :: num_intersections
    real(c_double), allocatable, intent(inout) :: intersections(:, :)
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    integer(c_int) :: num_nodes1, num_nodes2
    real(c_double) :: refined_s, refined_t
    logical(c_bool) :: does_intersect

    num_nodes1 = size(first%nodes, 2)
    num_nodes2 = size(second%nodes, 2)

    call from_linearized( &
         linearization_error1, first, num_nodes1, root_nodes1, &
         linearization_error2, second, num_nodes2, root_nodes2, &
         refined_s, refined_t, does_intersect, status)

    if (status /= Status_SUCCESS) then
       return
    end if

    if (.NOT. does_intersect) then
       return
    end if

    call add_intersection( &
         refined_s, refined_t, num_intersections, intersections)

  end subroutine add_from_linearized

  subroutine endpoint_check( &
       first, node_first, s, second, node_second, t, &
       num_intersections, intersections)

    ! NOTE: This is **explicitly** not intended for C inter-op.

    type(CurveData), intent(in) :: first
    real(c_double), intent(in) :: node_first(2)
    real(c_double), intent(in) :: s
    type(CurveData), intent(in) :: second
    real(c_double), intent(in) :: node_second(2)
    real(c_double), intent(in) :: t
    integer(c_int), intent(inout) :: num_intersections
    real(c_double), allocatable, intent(inout) :: intersections(:, :)
    ! Variables outside of signature.
    real(c_double) :: orig_s, orig_t

    if (.NOT. vector_close(2, node_first, node_second, VECTOR_CLOSE_EPS)) then
       return
    end if

    orig_s = (1 - s) * first%start + s * first%end_
    orig_t = (1 - t) * second%start + t * second%end_
    call add_intersection(orig_s, orig_t, num_intersections, intersections)

  end subroutine endpoint_check

  subroutine tangent_bbox_intersection( &
       first, second, num_intersections, intersections)

    ! NOTE: This is **explicitly** not intended for C inter-op.

    type(CurveData), intent(in) :: first, second
    integer(c_int), intent(inout) :: num_intersections
    real(c_double), allocatable, intent(inout) :: intersections(:, :)
    ! Variables outside of signature.
    integer(c_int) :: num_nodes1, num_nodes2
    real(c_double) :: nodes1_start(2), nodes1_end(2)
    real(c_double) :: nodes2_start(2), nodes2_end(2)

    num_nodes1 = size(first%nodes, 2)
    nodes1_start = first%nodes(:2, 1)
    nodes2_start = first%nodes(:2, num_nodes1)

    num_nodes2 = size(second%nodes, 2)
    nodes1_end = second%nodes(:2, 1)
    nodes2_end = second%nodes(:2, num_nodes2)

    call endpoint_check( &
         first, nodes1_start, 0.0_dp, &
         second, nodes1_end, 0.0_dp, &
         num_intersections, intersections)
    call endpoint_check( &
         first, nodes1_start, 0.0_dp, &
         second, nodes2_end, 1.0_dp, &
         num_intersections, intersections)
    call endpoint_check( &
         first, nodes2_start, 1.0_dp, &
         second, nodes1_end, 0.0_dp, &
         num_intersections, intersections)
    call endpoint_check( &
         first, nodes2_start, 1.0_dp, &
         second, nodes2_end, 1.0_dp, &
         num_intersections, intersections)

  end subroutine tangent_bbox_intersection

  subroutine add_candidates( &
       candidates, num_candidates, first, second, enum_)

    ! Helper for ``intersect_one_round``.

    type(CurveData), allocatable, intent(inout) :: candidates(:, :)
    integer(c_int), intent(inout) :: num_candidates
    type(CurveData), intent(in) :: first, second
    integer(c_int), intent(in) :: enum_
    ! Variables outside of signature.
    integer(c_int) :: curr_size
    type(CurveData), allocatable :: candidates_swap(:, :)

    ! First, update the number of candidates.
    if (enum_ == Subdivide_FIRST .OR. enum_ == Subdivide_SECOND) then
       num_candidates = num_candidates + 2
    else if (enum_ == Subdivide_BOTH) then
       num_candidates = num_candidates + 4
    else
       return
    end if

    if (allocated(candidates)) then
       curr_size = size(candidates, 2)
       if (curr_size < num_candidates) then
          allocate(candidates_swap(2, num_candidates))
          ! NOTE: This assumes, but does not check ``candidates`` has two rows.
          candidates_swap(:, :curr_size) = candidates(:, :curr_size)
          call move_alloc(candidates_swap, candidates)
       end if
    else
       allocate(candidates(2, num_candidates))
    end if

    ! NOTE: This assumes, but does not check, that ``candidates`` is MxN with
    !       M == 2 and ``num_candidates`` <= N.

    if (enum_ == Subdivide_FIRST) then
       candidates(2, num_candidates - 1) = second
       candidates(2, num_candidates) = second
       call subdivide_curve( &
            first, &
            candidates(1, num_candidates - 1), &
            candidates(1, num_candidates))
    else if (enum_ == Subdivide_SECOND) then
       candidates(1, num_candidates - 1) = first
       candidates(1, num_candidates) = first
       call subdivide_curve( &
            second, &
            candidates(2, num_candidates - 1), &
            candidates(2, num_candidates))
    else if (enum_ == Subdivide_BOTH) then
       call subdivide_curve( &
            first, &
            candidates(1, num_candidates - 3), &
            candidates(1, num_candidates - 1))
       call subdivide_curve( &
            second, &
            candidates(2, num_candidates - 3), &
            candidates(2, num_candidates - 2))

       candidates(1, num_candidates - 2) = candidates(1, num_candidates - 3)
       candidates(2, num_candidates - 1) = candidates(2, num_candidates - 3)
       candidates(1, num_candidates) = candidates(1, num_candidates - 1)
       candidates(2, num_candidates) = candidates(2, num_candidates - 2)
    end if

  end subroutine add_candidates

  subroutine intersect_one_round( &
       root_nodes_first, root_nodes_second, num_candidates, candidates, &
       num_intersections, intersections, &
       next_candidates, num_next_candidates, status)

    ! NOTE: This is **explicitly** not intended for C inter-op.
    ! NOTE: This assumes, but does not check, that ``candidates`` has
    !       two rows and has at **least** ``num_candidates`` columns.

    ! Possible error states:
    ! * Status_SUCCESS    : On success.
    ! * Status_PARALLEL   : Via ``add_from_linearized()``.

    real(c_double), intent(in) :: root_nodes_first(:, :)
    real(c_double), intent(in) :: root_nodes_second(:, :)
    integer(c_int), intent(in) :: num_candidates
    type(CurveData), intent(in) :: candidates(:, :)
    integer(c_int), intent(inout) :: num_intersections
    real(c_double), allocatable, intent(inout) :: intersections(:, :)
    type(CurveData), allocatable, intent(inout) :: next_candidates(:, :)
    integer(c_int), intent(out) :: num_next_candidates
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    type(CurveData) :: first, second
    real(c_double) :: linearization_error1, linearization_error2
    integer(c_int) :: bbox_int, index_, num_nodes1, num_nodes2
    integer(c_int) :: subdivide_enum

    num_next_candidates = 0
    status = Status_SUCCESS
    subdivide_loop: do index_ = 1, num_candidates
       ! NOTE: We **hope** that the compiler avoids turning this alias (for
       !       the sake of typing fewer characters) into a copy.
       first = candidates(1, index_)
       second = candidates(2, index_)
       ! Compute the linearization error for each curve.
       num_nodes1 = size(first%nodes, 2)
       num_nodes2 = size(second%nodes, 2)
       call linearization_error( &
            num_nodes1, 2, first%nodes, linearization_error1)
       call linearization_error( &
            num_nodes2, 2, second%nodes, linearization_error2)

       if (linearization_error1 < LINEARIZATION_THRESHOLD) then
          if (linearization_error2 < LINEARIZATION_THRESHOLD) then
             ! If both ``first`` and ``second`` are linearizations, then
             ! we can (attempt to) intersect them immediately.
             subdivide_enum = Subdivide_NEITHER
             call add_from_linearized( &
                  first, root_nodes_first, linearization_error1, &
                  second, root_nodes_second, linearization_error2, &
                  num_intersections, intersections, status)

             ! If there was a failure, exit this subroutine.
             if (status /= Status_SUCCESS) then
                return
             end if

             ! If there was no failure, move to the next iteration.
             cycle subdivide_loop
          else
             subdivide_enum = Subdivide_SECOND
             call bbox_line_intersect( &
                  num_nodes2, second%nodes, &
                  first%nodes(:, 1), first%nodes(:, num_nodes1), bbox_int)
          end if
       else
          if (linearization_error2 < LINEARIZATION_THRESHOLD) then
             subdivide_enum = Subdivide_FIRST
             call bbox_line_intersect( &
                  num_nodes1, first%nodes, &
                  second%nodes(:, 1), second%nodes(:, num_nodes2), bbox_int)
          else
             subdivide_enum = Subdivide_BOTH
             ! If neither curve is close to a line, we can still reject the
             ! pair if the bounding boxes are disjoint.
             call bbox_intersect( &
                  num_nodes1, first%nodes, num_nodes2, second%nodes, bbox_int)
          end if
       end if

       ! Reject if the bounding boxes do not intersect.
       if (bbox_int == BoxIntersectionType_DISJOINT) then
          cycle subdivide_loop
       else if (bbox_int == BoxIntersectionType_TANGENT) then
          call tangent_bbox_intersection( &
               first, second, num_intersections, intersections)
          cycle subdivide_loop
       end if

       ! If we haven't ``cycle``-d this iteration, add the
       ! ``next_candidates`` pair.
       call add_candidates( &
            next_candidates, num_next_candidates, &
            first, second, subdivide_enum)
    end do subdivide_loop

  end subroutine intersect_one_round

  subroutine make_candidates( &
       nodes_first, nodes_second, candidates)

    ! NOTE: This is a (private) helper for ``all_intersections``.

    real(c_double), intent(in) :: nodes_first(:, :)
    real(c_double), intent(in) :: nodes_second(:, :)
    type(CurveData), allocatable, intent(inout) :: candidates(:, :)

    ! NOTE: This assumes (but does not verify) that if ``candidates`` has
    !       been allocated, then it is ``2 x N`` with ``N > 0``.
    if (.NOT. allocated(candidates)) then
       allocate(candidates(2, 1))
    end if

    ! NOTE: Since we **might** be re-using ``candidates``, we can't rely on
    !       the default values of all fields being intact.
    candidates(1, 1)%start = 0.0_dp
    candidates(1, 1)%end_ = 1.0_dp
    candidates(1, 1)%nodes = nodes_first
    candidates(2, 1)%start = 0.0_dp
    candidates(2, 1)%end_ = 1.0_dp
    candidates(2, 1)%nodes = nodes_second

  end subroutine make_candidates

  subroutine prune_candidates(candidates, num_candidates)

    ! Uses more strict bounding box intersection predicate by forming the
    ! actual convex hull of each candidate curve segment and then checking
    ! if those convex hulls collide.

    type(CurveData), allocatable, intent(inout) :: candidates(:, :)
    integer(c_int), intent(inout) :: num_candidates
    ! Variables outside of signature.
    integer(c_int) :: accepted, i
    integer(c_int) :: num_nodes1, num_nodes2
    integer(c_int) :: polygon_size1, polygon_size2
    real(c_double), allocatable :: polygon1(:, :), polygon2(:, :)
    logical(c_bool) :: collision

    accepted = 0
    ! NOTE: This **assumes** but does not check that ``candidates`` is
    !       allocated and size ``2 x NC`` where ``NC >= num_candidates``.
    do i = 1, num_candidates
       ! NOTE: This **assumes** that ``%nodes`` is allocated and size
       !       ``N x 2``.
       num_nodes1 = size(candidates(1, i)%nodes, 2)
       if (allocated(polygon1)) then
          if (size(polygon1, 2) < num_nodes1) then
             ! NOTE: We "exclude" this block from ``lcov`` because it
             !       **should** never occur. (All candidates in row 1 should
             !       have the same ``num_nodes``).
             ! LCOV_EXCL_START
             deallocate(polygon1)
             allocate(polygon1(2, num_nodes1))
             ! LCOV_EXCL_STOP
          end if
       else
          allocate(polygon1(2, num_nodes1))
       end if
       call convex_hull( &
            num_nodes1, candidates(1, i)%nodes, &
            polygon_size1, polygon1(:, :num_nodes1))

       ! NOTE: This **assumes** that ``%nodes`` is allocated and size
       !       ``N x 2``.
       num_nodes2 = size(candidates(2, i)%nodes, 2)
       if (allocated(polygon2)) then
          if (size(polygon2, 2) < num_nodes2) then
             ! NOTE: We "exclude" this block from ``lcov`` because it
             !       **should** never occur. (All candidates in row 2 should
             !       have the same ``num_nodes``).
             ! LCOV_EXCL_START
             deallocate(polygon2)
             allocate(polygon2(2, num_nodes2))
             ! LCOV_EXCL_STOP
          end if
       else
          allocate(polygon2(2, num_nodes2))
       end if
       call convex_hull( &
            num_nodes2, candidates(2, i)%nodes, &
            polygon_size2, polygon2(:, :num_nodes2))

       ! Now check if the convex hulls actually collide.
       call polygon_collide( &
            polygon_size1, polygon1(:, :polygon_size1), &
            polygon_size2, polygon2(:, :polygon_size2), collision)
       if (collision) then
          accepted = accepted + 1
          ! NOTE: This relies on the invariant ``accepted <= i``. In the
          !       ``accepted == i`` case, there is nothing to do.
          if (accepted < i) then
             candidates(1, accepted) = candidates(1, i)
             candidates(2, accepted) = candidates(2, i)
          end if
       end if
    end do
    num_candidates = accepted

  end subroutine prune_candidates

  subroutine elevate_helper( &
       curr_size, nodes, final_size, workspace, elevated)

    ! NOTE: This is a helper for ``make_same_degree``.
    ! NOTE: This assumes, but does not check, that ``final_size > curr_size``.

    integer(c_int), intent(in) :: curr_size
    real(c_double), intent(in) :: nodes(2, curr_size)
    integer(c_int), intent(in) :: final_size
    real(c_double), intent(inout) :: workspace(2, final_size)
    real(c_double), intent(inout) :: elevated(2, final_size)
    ! Variables outside of signature.
    logical(c_bool) :: to_workspace
    integer(c_int) :: i

    ! Elevate ``nodes`` by using both ``workspace`` and ``elevated`` as
    ! workspaces.
    if (mod(final_size - curr_size, 2) == 1) then
       ! If we have an odd number of steps, then we'll write to e1
       ! first, for example:
       ! 1 step : ws -> el
       ! 3 steps: ws -> el -> ws -> el
       workspace(:, :curr_size) = nodes
       to_workspace = .FALSE.
    else
       ! If we have an odd number of steps, then we'll write to ws
       ! first, for example:
       ! 2 steps: el -> ws -> el
       ! 4 steps: el -> ws -> el -> ws -> el
       elevated(:, :curr_size) = nodes
       to_workspace = .TRUE.
    end if

    do i = curr_size, final_size - 1
       if (to_workspace) then
          call elevate_nodes( &
               i, 2, elevated(:, :i), workspace(:, :i + 1))
       else
          call elevate_nodes( &
               i, 2, workspace(:, :i), elevated(:, :i + 1))
       end if

       to_workspace = .NOT. to_workspace  ! Switch parity.
    end do

  end subroutine elevate_helper

  subroutine make_same_degree( &
       num_nodes1, nodes1, num_nodes2, nodes2, &
       num_nodes, elevated1, elevated2)

    ! NOTE: This subroutine is not part of the C ABI for this module,
    !       but it is (for now) public, so that it can be tested.

    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(2, num_nodes1)
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(2, num_nodes2)
    integer(c_int), intent(out) :: num_nodes
    real(c_double), allocatable, intent(out) :: elevated1(:, :)
    real(c_double), allocatable, intent(out) :: elevated2(:, :)

    if (num_nodes1 > num_nodes2) then
       num_nodes = num_nodes1
       allocate(elevated1(2, num_nodes1))
       allocate(elevated2(2, num_nodes1))
       ! Populate ``elevated2`` by elevating (using ``elevated1`` as a
       ! helper workspace).
       call elevate_helper( &
            num_nodes2, nodes2, num_nodes1, elevated1, elevated2)
       ! Populate ``elevated1`` without re-allocating it.
       elevated1(:, :) = nodes1
    else if (num_nodes2 > num_nodes1) then
       num_nodes = num_nodes2
       allocate(elevated1(2, num_nodes2))
       allocate(elevated2(2, num_nodes2))
       ! Populate ``elevated1`` by elevating (using ``elevated2`` as a
       ! helper workspace).
       call elevate_helper( &
            num_nodes1, nodes1, num_nodes2, elevated2, elevated1)
       ! Populate ``elevated2`` without re-allocating it.
       elevated2(:, :) = nodes2
    else
       num_nodes = num_nodes2
       elevated1 = nodes1
       elevated2 = nodes2
    end if

  end subroutine make_same_degree

  subroutine add_coincident_parameters( &
       num_nodes1, nodes1, num_nodes2, nodes2, &
       num_intersections, intersections, coincident)

    ! NOTE: This subroutine is not part of the C ABI for this module,
    !       but it is (for now) public, so that it can be tested.

    ! NOTE: If ``locate_point()`` fails with ``LOCATE_INVALID`` when trying
    !       to place **any** of the curve endpoints on the other segment,
    !       this method will just return ``coincident == FALSE``. This
    !       is "intentional" in the sense that callers (likely) won't benefit
    !       from a **different** error status than the one they are already
    !       trying to avoid by showing the curves are coincident.

    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(2, num_nodes1)
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(2, num_nodes2)
    integer(c_int), intent(inout) :: num_intersections
    real(c_double), allocatable, intent(inout) :: intersections(:, :)
    logical(c_bool), intent(out) :: coincident
    ! Variables outside of signature.
    real(c_double), target, allocatable :: elevated1(:, :)
    real(c_double), target, allocatable :: elevated2(:, :)
    real(c_double), target, allocatable :: specialized(:, :)
    integer(c_int) :: num_nodes
    real(c_double) :: point(2)
    real(c_double) :: s_initial, s_final, t_initial, t_final
    real(c_double), pointer :: as_vec1(:), as_vec2(:)

    coincident = .FALSE.
    ! First, make sure the nodes are the same degree.
    call make_same_degree( &
         num_nodes1, nodes1, num_nodes2, nodes2, &
         num_nodes, elevated1, elevated2)

    point = nodes2(:, 1)
    call locate_point( &
         num_nodes, 2, elevated1, point, s_initial)
    point = nodes2(:, num_nodes2)
    call locate_point( &
         num_nodes, 2, elevated1, point, s_final)
    ! Bail out if the "locate" failed.
    if (s_initial == LOCATE_INVALID .OR. s_final == LOCATE_INVALID) then
       return
    end if

    if (s_initial /= LOCATE_MISS .AND. s_final /= LOCATE_MISS) then
       ! In this case, if the curves were coincident, then ``curve2``
       ! would be "fully" contained in ``curve1``, so we specialize
       ! ``curve1`` down to that interval to check.
       allocate(specialized(2, num_nodes))
       call specialize_curve( &
            num_nodes, 2, elevated1, s_initial, s_final, specialized)
       call c_f_pointer(c_loc(specialized), as_vec1, [2 * num_nodes])
       call c_f_pointer(c_loc(elevated2), as_vec2, [2 * num_nodes])

       if (vector_close( &
            2 * num_nodes, as_vec1, as_vec2, VECTOR_CLOSE_EPS)) then
          coincident = .TRUE.
          ! Empty out candidates since anything in there until now was
          ! "accidental".
          num_intersections = 0
          call add_intersection( &
               s_initial, 0.0_dp, num_intersections, intersections)
          call add_intersection( &
               s_final, 1.0_dp, num_intersections, intersections)
       end if

       ! In either case (``vector_close()`` or not), we are done.
       return
    end if

    point = nodes1(:, 1)
    call locate_point( &
         num_nodes, 2, elevated2, point, t_initial)
    point = nodes1(:, num_nodes1)
    call locate_point( &
         num_nodes, 2, elevated2, point, t_final)
    ! Bail out if the "locate" failed.
    if (t_initial == LOCATE_INVALID .OR. t_final == LOCATE_INVALID) then
       return
    end if

    if (t_initial == LOCATE_MISS .AND. t_final == LOCATE_MISS) then
       ! An overlap must have two endpoints and since at most one of the
       ! endpoints of ``curve2`` lies on ``curve1`` (as indicated by at
       ! least one of the ``s``-parameters being ``None``), we need (at least)
       ! one endpoint of ``curve1`` on ``curve2``.
       return
    end if

    if (t_initial /= LOCATE_MISS .AND. t_final /= LOCATE_MISS) then
       ! In this case, if the curves were coincident, then ``curve1``
       ! would be "fully" contained in ``curve2``, so we specialize
       ! ``curve2`` down to that interval to check.
       allocate(specialized(2, num_nodes))
       call specialize_curve( &
            num_nodes, 2, elevated2, t_initial, t_final, specialized)
       call c_f_pointer(c_loc(elevated1), as_vec1, [2 * num_nodes])
       call c_f_pointer(c_loc(specialized), as_vec2, [2 * num_nodes])

       if (vector_close( &
            2 * num_nodes, as_vec1, as_vec2, VECTOR_CLOSE_EPS)) then
          coincident = .TRUE.
          ! Empty out candidates since anything in there until now was
          ! "accidental".
          num_intersections = 0
          call add_intersection( &
               0.0_dp, t_initial, num_intersections, intersections)
          call add_intersection( &
               1.0_dp, t_final, num_intersections, intersections)
       end if

       ! In either case (``vector_close()`` or not), we are done.
       return
    end if

    if (s_initial == LOCATE_MISS .AND. s_final == LOCATE_MISS) then
       ! An overlap must have two endpoints and since exactly one of the
       ! endpoints of ``curve1`` lies on ``curve2`` (as indicated by exactly
       ! one of the ``t``-parameters being ``None``), we need (at least)
       ! one endpoint of ``curve1`` on ``curve2``.
       return
    end if

    ! At this point, we know exactly one of the ``s``-parameters and exactly
    ! one of the ``t``-parameters is not ``None``. So we fill in all four
    ! of ``(s|t)_(initial|final)`` with the known values.
    if (s_initial == LOCATE_MISS) then
       if (t_initial == LOCATE_MISS) then
          ! B1(s_final) = B2(1) AND B1(1) = B2(t_final)
          s_initial = s_final
          s_final = 1.0_dp
          t_initial = 1.0_dp
          ! t_final is fine as-is
       else
          ! B1(0) = B2(t_initial) AND B1(s_final) = B2(1)
          s_initial = 0.0_dp
          ! s_final is fine as-is
          ! t_initial is fine as-is
          t_final = 1.0_dp
       end if
    else
       if (t_initial == LOCATE_MISS) then
          ! B1(s_initial) = B2(0) AND B1(1 ) = B2(t_final)
          ! s_initial is fine as-is
          s_final = 1.0_dp
          t_initial = 0.0_dp
          ! t_final is fine as-is
       else
          ! B1(0) = B2(t_initial) AND B1(s_initial) = B2(0)
          s_final = s_initial
          s_initial = 0.0_dp
          ! t_initial is fine as-is
          t_final = 0.0_dp
       end if
    end if

    if ( &
         abs(s_initial - s_final) < MIN_INTERVAL_WIDTH .AND. &
         abs(t_initial - t_final) < MIN_INTERVAL_WIDTH) then
       return
    end if

    allocate(specialized(2, num_nodes))
    ! First specialize ``elevated1`` onto ``specialized``.
    call specialize_curve( &
         num_nodes, 2, elevated1, s_initial, s_final, specialized)
    ! Then specialize ``elevated2`` onto ``elevated1`` since we are
    ! done with it.
    call specialize_curve( &
         num_nodes, 2, elevated2, t_initial, t_final, elevated1)

    call c_f_pointer(c_loc(specialized), as_vec1, [2 * num_nodes])
    call c_f_pointer(c_loc(elevated1), as_vec2, [2 * num_nodes])

    if (vector_close( &
         2 * num_nodes, as_vec1, as_vec2, VECTOR_CLOSE_EPS)) then
       coincident = .TRUE.
       ! Empty out candidates since anything in there until now was
       ! "accidental".
       num_intersections = 0
       call add_intersection( &
            s_initial, t_initial, num_intersections, intersections)
       call add_intersection( &
            s_final, t_final, num_intersections, intersections)
    end if

  end subroutine add_coincident_parameters

  subroutine all_intersections( &
       num_nodes_first, nodes_first, num_nodes_second, nodes_second, &
       intersections, num_intersections, status)

    ! NOTE: This is **explicitly** not intended for C inter-op, but
    !       a C compatible interface is exposed as ``all_intersections_abi``.

    ! Possible error states:
    ! * Status_SUCCESS       : On success.
    ! * Status_PARALLEL      : Via ``intersect_one_round()``.
    ! * Status_NO_CONVERGE   : If the curves don't converge to linear after
    !                          ``MAX_INTERSECT_SUBDIVISIONS``.
    ! * (N >= MAX_CANDIDATES): The number of candidates if it exceeds the limit
    !                          ``MAX_CANDIDATES`` (64 is the default).

    integer(c_int), intent(in) :: num_nodes_first
    real(c_double), intent(in) :: nodes_first(2, num_nodes_first)
    integer(c_int), intent(in) :: num_nodes_second
    real(c_double), intent(in) :: nodes_second(2, num_nodes_second)
    real(c_double), allocatable, intent(inout) :: intersections(:, :)
    integer(c_int), intent(out) :: num_intersections
    integer(c_int), intent(out) :: status
    ! Variables outside of signature.
    integer(c_int) :: num_candidates, num_next_candidates
    integer(c_int) :: index_, intersect_status
    logical(c_bool) :: is_even, coincident

    num_intersections = 0
    ! First iteration is odd (i.e. ``index_ == 1``).
    num_candidates = 1
    call make_candidates( &
         nodes_first, nodes_second, CANDIDATES_ODD)
    status = Status_SUCCESS  ! Default.

    is_even = .TRUE.  ! At zero.
    do index_ = 1, MAX_INTERSECT_SUBDIVISIONS
       is_even = .NOT. is_even  ! Switch parity.

       if (is_even) then
          ! Since ``index_`` is even, we READ from ``CANDIDATES_EVEN``
          ! and WRITE to ``CANDIDATES_ODD``.
          call intersect_one_round( &
               nodes_first, nodes_second, &
               num_candidates, CANDIDATES_EVEN, &
               num_intersections, intersections, &
               CANDIDATES_ODD, num_next_candidates, intersect_status)
       else
          ! Since ``index_`` is odd, we READ from ``CANDIDATES_ODD``
          ! and WRITE to ``CANDIDATES_EVEN``.
          call intersect_one_round( &
               nodes_first, nodes_second, &
               num_candidates, CANDIDATES_ODD, &
               num_intersections, intersections, &
               CANDIDATES_EVEN, num_next_candidates, intersect_status)
       end if

       ! NOTE: This only checks for two error statuses from
       !       ``intersect_one_round()``, so it is inherently brittle
       !       to changes there.
       if (intersect_status /= Status_SUCCESS) then
          status = intersect_status
          return
       end if

       ! Update the number of candidates.
       num_candidates = num_next_candidates

       ! Bail out of there are too many candidates.
       if (num_candidates > MAX_CANDIDATES) then
          if (is_even) then
             call prune_candidates(CANDIDATES_ODD, num_candidates)
          else
             call prune_candidates(CANDIDATES_EVEN, num_candidates)
          end if
          ! If pruning didn't fix anything, we check if the curves are
          ! coincident and "fail" if they aren't.
          if (num_candidates > MAX_CANDIDATES) then
             call add_coincident_parameters( &
                  num_nodes_first, nodes_first, &
                  num_nodes_second, nodes_second, &
                  num_intersections, intersections, coincident)
             if (.NOT. coincident) then
                ! NOTE: This assumes that all of the status enum values are
                !       less than ``MAX_CANDIDATES + 1``.
                ! NOTE: This line is very difficult to trigger due to the
                !       ``prune_candidates()`` and coincident check
                !       mitigations. As a result, there is no unit test
                !       to trigger this line (no case has been discovered
                !       yet).
                status = num_candidates  ! LCOV_EXCL_LINE
             end if
             ! Return either way since pruning didn't "fix" anything.
             return
          end if
       end if

       ! If none of the candidate pairs have been accepted, then there are
       ! no more intersections to find.
       if (num_candidates == 0) then
          return
       end if
    end do

    ! If we've reached this point, then the curve intersection failed to
    ! converge to "approximately linear" subdivided curves after
    ! ``MAX_INTERSECT_SUBDIVISIONS``.
    status = Status_NO_CONVERGE

  end subroutine all_intersections

  subroutine all_intersections_abi( &
       num_nodes_first, nodes_first, num_nodes_second, nodes_second, &
       intersections_size, intersections, num_intersections, status) &
       bind(c, name='curve_intersections')

    ! NOTE: The number of intersections cannot be known beforehand (though it
    !       will be **at most** the product of the degrees of the two curves,
    !       by Bezout's theorem). If ``intersections`` is not large enough
    !       (i.e. if ``intersections_size`` < num_intersections``), then
    !       ``intersections`` will not be populated and the ``status`` will be
    !       set to ``Status_INSUFFICIENT_SPACE``. However, the value of
    !       ``num_intersections`` will be accurate and ``intersections``
    !       should be re-sized by the caller to accommodate.

    ! Possible error states:
    ! * Status_SUCCESS           : On success.
    ! * Status_PARALLEL          : Via ``all_intersections()``.
    ! * Status_NO_CONVERGE       : Via ``all_intersections()``.
    ! * Status_INSUFFICIENT_SPACE: If ``intersections_size`` is smaller than
    !                              the number of intersections.
    ! * (N >= MAX_CANDIDATES)    : Via ``all_intersections()``.

    integer(c_int), intent(in) :: num_nodes_first
    real(c_double), intent(in) :: nodes_first(2, num_nodes_first)
    integer(c_int), intent(in) :: num_nodes_second
    real(c_double), intent(in) :: nodes_second(2, num_nodes_second)
    integer(c_int), intent(in) :: intersections_size
    real(c_double), intent(out) :: intersections(2, intersections_size)
    integer(c_int), intent(out) :: num_intersections
    integer(c_int), intent(out) :: status

    call all_intersections( &
         num_nodes_first, nodes_first, num_nodes_second, nodes_second, &
         INTERSECTIONS_WORKSPACE, num_intersections, status)

    if (status /= Status_SUCCESS) then
       return
    end if

    if (num_intersections > intersections_size) then
       status = Status_INSUFFICIENT_SPACE
    else if (num_intersections > 0) then
       ! NOTE: This assumes, but doesn't check that
       !       ``INTERSECTIONS_WORKSPACE`` has been allocated
       !       and that is has 2 rows, just like ``intersections``.
       intersections(:, :num_intersections) = ( &
            INTERSECTIONS_WORKSPACE(:, :num_intersections))
    end if

  end subroutine all_intersections_abi

  subroutine set_max_candidates(num_candidates) &
       bind(c, name='set_max_candidates')

    integer(c_int), intent(in) :: num_candidates

    MAX_CANDIDATES = num_candidates

  end subroutine set_max_candidates

  subroutine get_max_candidates(num_candidates) &
       bind(c, name='get_max_candidates')

    integer(c_int), intent(out) :: num_candidates

    num_candidates = MAX_CANDIDATES

  end subroutine get_max_candidates

  subroutine set_similar_ulps(num_bits) &
       bind(c, name='set_similar_ulps')

    integer(c_int), intent(in) :: num_bits

    SIMILAR_ULPS = num_bits

  end subroutine set_similar_ulps

  subroutine get_similar_ulps(num_bits) &
       bind(c, name='get_similar_ulps')

    integer(c_int), intent(out) :: num_bits

    num_bits = SIMILAR_ULPS

  end subroutine get_similar_ulps

  subroutine free_curve_intersections_workspace() &
       bind(c, name='free_curve_intersections_workspace')

    ! NOTE: This **should** be run during clean-up for any code which
    !       invokes ``all_intersections()``.

    if (allocated(CANDIDATES_ODD)) then
       deallocate(CANDIDATES_ODD)
    end if

    if (allocated(CANDIDATES_EVEN)) then
       deallocate(CANDIDATES_EVEN)
    end if

    if (allocated(INTERSECTIONS_WORKSPACE)) then
       deallocate(INTERSECTIONS_WORKSPACE)
    end if

  end subroutine free_curve_intersections_workspace

end module curve_intersection
