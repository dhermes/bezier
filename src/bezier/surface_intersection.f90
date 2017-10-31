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

  use, intrinsic :: iso_c_binding, only: c_double, c_int
  use types, only: dp
  use surface, only: evaluate_barycentric, jacobian_both
  implicit none
  private newton_refine_solve
  public newton_refine

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

end module surface_intersection
