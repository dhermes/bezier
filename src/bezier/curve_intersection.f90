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

  use, intrinsic :: iso_c_binding, only: c_double, c_int, c_bool
  use types, only: dp
  use helpers, only: cross_product, bbox, wiggle_interval, in_interval
  use curve, only: CurveData, evaluate_multi, evaluate_hodograph
  implicit none
  private
  public &
       Intersection, BoxIntersectionType_INTERSECTION, &
       BoxIntersectionType_TANGENT, BoxIntersectionType_DISJOINT, &
       FROM_LINEARIZED_SUCCESS, FROM_LINEARIZED_PARALLEL, &
       FROM_LINEARIZED_WIGGLE_FAIL, linearization_error, &
       segment_intersection, newton_refine_intersect, bbox_intersect, &
       parallel_different, from_linearized, bbox_line_intersect, &
       add_intersection

  ! NOTE: This (for now) is not meant to be C-interoperable.
  type :: Intersection
     real(c_double) :: s = -1.0_dp
     real(c_double) :: t = -1.0_dp
     type(CurveData), pointer :: first => null()
     type(CurveData), pointer :: second => null()
  end type Intersection

  integer(c_int), parameter :: BoxIntersectionType_INTERSECTION = 0
  integer(c_int), parameter :: BoxIntersectionType_TANGENT = 1
  integer(c_int), parameter :: BoxIntersectionType_DISJOINT = 2
  integer(c_int), parameter :: FROM_LINEARIZED_SUCCESS = 0
  ! In ``from_linearized``, ``py_exc == 1`` corresponds to
  ! ``NotImplementedError('Line segments parallel.')``
  integer(c_int), parameter :: FROM_LINEARIZED_PARALLEL = 1
  ! In ``from_linearized``, ``py_exc == 2`` indicates that
  ! ``wiggle_interval`` failed.
  integer(c_int), parameter :: FROM_LINEARIZED_WIGGLE_FAIL = 2

contains

  subroutine linearization_error( &
       num_nodes, dimension_, nodes, error) &
       bind(c, name='linearization_error')

    integer(c_int), intent(in) :: num_nodes, dimension_
    real(c_double), intent(in) :: nodes(num_nodes, dimension_)
    real(c_double), intent(out) :: error
    ! Variables outside of signature.
    real(c_double) :: second_deriv(num_nodes - 2, dimension_)
    real(c_double) :: worst_case(dimension_)

    ! A line has no linearization error.
    if (num_nodes == 2) then
       error = 0.0_dp
       return
    end if

    second_deriv = ( &
         nodes(:num_nodes - 2, :) - &
         2.0_dp * nodes(2:num_nodes - 1, :) + &
         nodes(3:, :))
    worst_case = maxval(abs(second_deriv), 1)
    error = 0.125_dp * (num_nodes - 1) * (num_nodes - 2) * norm2(worst_case)

  end subroutine linearization_error

  subroutine segment_intersection( &
       start0, end0, start1, end1, s, t, success) &
       bind(c, name='segment_intersection')

    real(c_double), intent(in) :: start0(1, 2)
    real(c_double), intent(in) :: end0(1, 2)
    real(c_double), intent(in) :: start1(1, 2)
    real(c_double), intent(in) :: end1(1, 2)
    real(c_double), intent(out) :: s, t
    logical(c_bool), intent(out) :: success
    ! Variables outside of signature.
    real(c_double) :: delta0(1, 2)
    real(c_double) :: delta1(1, 2)
    real(c_double) :: start_delta(1, 2)
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
       bind(c, name='newton_refine_intersect')

    real(c_double), intent(in) :: s
    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(num_nodes1, 2)
    real(c_double), intent(in) :: t
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(num_nodes2, 2)
    real(c_double), intent(out) :: new_s, new_t
    ! Variables outside of signature.
    real(c_double) :: param(1)
    real(c_double) :: func_val(1, 2)
    real(c_double) :: workspace(1, 2)
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
         jac_mat(2, 1) * func_val(1, 2)) / determinant)
    new_s = s + delta_s

    delta_t = ( &
         (jac_mat(1, 2) * func_val(1, 1) - &
         jac_mat(1, 1) * func_val(1, 2)) / determinant)
    new_t = t + delta_t

  end subroutine newton_refine_intersect

  subroutine bbox_intersect( &
       num_nodes1, nodes1, num_nodes2, nodes2, enum_) &
       bind(c, name='bbox_intersect')

    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(num_nodes1, 2)
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(num_nodes2, 2)
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
       start0, end0, start1, end1, result_) &
       bind(c, name='parallel_different')

    real(c_double), intent(in) :: start0(1, 2)
    real(c_double), intent(in) :: end0(1, 2)
    real(c_double), intent(in) :: start1(1, 2)
    real(c_double), intent(in) :: end1(1, 2)
    logical(c_bool), intent(out) :: result_
    ! Variables outside of signature.
    real(c_double) :: delta0(1, 2)
    real(c_double) :: val1, val2  ! Workspace

    delta0 = end0 - start0
    call cross_product(start0, delta0, val1)  ! line0_const
    call cross_product(start1, delta0, val2)  ! start1_against

    if (val1 /= val2) then
       result_ = .TRUE.
       return
    end if

    val1 = dot_product(delta0(1, :), delta0(1, :))  ! norm0_sq
    ! val2 == start_numer
    val2 = dot_product(start1(1, :) - start0(1, :), delta0(1, :))
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
    ! they may contain [0, 1] between them, so we make sure that 0
    ! isn't between them.
    result_ = (0.0_dp < min(val1, val2) .OR. max(val1, val2) < 0.0_dp)

  end subroutine parallel_different

  subroutine from_linearized( &
       error1, start1, end1, start_node1, end_node1, num_nodes1, nodes1, &
       error2, start2, end2, start_node2, end_node2, num_nodes2, nodes2, &
       refined_s, refined_t, does_intersect, py_exc) &
       bind(c, name='from_linearized')

    real(c_double), intent(in) :: error1, start1, end1
    real(c_double), intent(in) :: start_node1(1, 2)
    real(c_double), intent(in) :: end_node1(1, 2)
    integer(c_int), intent(in) :: num_nodes1
    real(c_double), intent(in) :: nodes1(num_nodes1, 2)
    real(c_double), intent(in) :: error2, start2, end2
    real(c_double), intent(in) :: start_node2(1, 2)
    real(c_double), intent(in) :: end_node2(1, 2)
    integer(c_int), intent(in) :: num_nodes2
    real(c_double), intent(in) :: nodes2(num_nodes2, 2)
    real(c_double), intent(out) :: refined_s, refined_t
    logical(c_bool), intent(out) :: does_intersect
    integer(c_int), intent(out) :: py_exc
    ! Variables outside of signature.
    real(c_double) :: s, t
    integer(c_int) :: enum_
    logical(c_bool) :: success

    py_exc = FROM_LINEARIZED_SUCCESS
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
               num_nodes1, nodes1, num_nodes2, nodes2, enum_)
          if (enum_ == BoxIntersectionType_DISJOINT) then
             does_intersect = .FALSE.
             return
          end if
       end if

       ! Expect the wrapper code to raise.
       py_exc = FROM_LINEARIZED_PARALLEL
       return
    end if

    does_intersect = .TRUE.
    ! Now, promote ``s`` and ``t`` onto the original curves.
    s = (1.0_dp - s) * start1 + s * end1  ! orig_s
    t = (1.0_dp - t) * start2 + t * end2  ! orig_t
    ! Perform one step of Newton iteration to refine the computed
    ! values of s and t.
    call newton_refine_intersect( &
         s, num_nodes1, nodes1, t, num_nodes2, nodes2, refined_s, refined_t)

    call wiggle_interval(refined_s, s, success)
    if (.NOT. success) then
       py_exc = FROM_LINEARIZED_WIGGLE_FAIL  ! LCOV_EXCL_LINE
       return  ! LCOV_EXCL_LINE
    end if
    refined_s = s

    call wiggle_interval(refined_t, t, success)
    if (.NOT. success) then
       py_exc = FROM_LINEARIZED_WIGGLE_FAIL  ! LCOV_EXCL_LINE
       return  ! LCOV_EXCL_LINE
    end if
    refined_t = t

  end subroutine from_linearized

  subroutine bbox_line_intersect( &
       num_nodes, nodes, line_start, line_end, enum_) &
       bind(c, name='bbox_line_intersect')

    integer(c_int), intent(in) :: num_nodes
    real(c_double), intent(in) :: nodes(num_nodes, 2)
    real(c_double), intent(in) :: line_start(1, 2)
    real(c_double), intent(in) :: line_end(1, 2)
    integer(c_int), intent(out) :: enum_
    ! Variables outside of signature.
    real(c_double) :: left, right, bottom, top
    real(c_double) :: segment_start(1, 2)
    real(c_double) :: segment_end(1, 2)
    real(c_double) :: s_curr, t_curr
    logical(c_bool) :: success

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

  subroutine add_intersection(first, s, second, t, intersections)

    ! Adds an intersection to list of ``intersections``.
    !
    ! **DOES NOT** accounts for repeated intersections, but should.

    type(CurveData), target, intent(in) :: first
    real(c_double), intent(in) :: s
    type(CurveData), target, intent(in) :: second
    real(c_double), intent(in) :: t
    type(Intersection), allocatable, intent(inout) :: intersections(:)
    ! Variables outside of signature.
    integer(c_int) :: num_intersections, index_
    type(Intersection), allocatable :: intersections_swap(:)

    if (allocated(intersections)) then
       ! NOTE: It's expensive to re-size just to add a single intersection,
       !       but we expect the number of intersections to be low. The
       !       trade-off is that we don't have to track the "length" of
       !       ``intersections`` separate from the capacity.
       num_intersections = size(intersections)
       allocate(intersections_swap(num_intersections + 1))
       intersections_swap(:num_intersections) = ( &
            intersections(:num_intersections))
       call move_alloc(intersections_swap, intersections)
       index_ = num_intersections + 1
    else
       allocate(intersections(1))
       index_ = 1
    end if

    intersections(index_)%s = s
    intersections(index_)%t = t
    intersections(index_)%first => first
    intersections(index_)%second => second

  end subroutine add_intersection

end module curve_intersection
