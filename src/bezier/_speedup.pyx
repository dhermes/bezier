#!python
#cython: boundscheck=False, wraparound=False, language_level=3
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Cython "wrapped" interface for Fortran ABI.

In particular, this module is split into sections for

- ``curve.f90``
- ``curve_intersection.f90``
- ``helpers.f90``
- ``surface.f90``
- ``surface_intersection.f90``

These are grouped in a **single** ``.pyx`` file to avoid unintended
programming errors caused by sharing the **same** object file across
multiple Cython-generated modules. (This is problematic if the object
files **expect** to share some global state.)
"""

import warnings

from libc.stdlib cimport free
from libc.stdlib cimport malloc
from libcpp cimport bool as bool_t
import numpy as np
from numpy cimport dtype as dtype_t
from numpy cimport ndarray as ndarray_t

cimport bezier._curve
cimport bezier._curve_intersection
from bezier._curve_intersection cimport BoxIntersectionType
cimport bezier._helpers
cimport bezier._status
cimport bezier._surface
cimport bezier._surface_intersection
from bezier._surface_intersection cimport CurvedPolygonSegment
from bezier._surface_intersection cimport SurfaceContained


cdef double EPS = 0.5**40
cdef double LOCATE_MISS = -1.0
cdef double LOCATE_INVALID = -2.0
cdef double[::1, :] CURVES_WORKSPACE = np.empty((2, 2), order="F")
cdef int[::1] SEGMENT_ENDS_WORKSPACE = np.empty(3, dtype=np.intc)
cdef dtype_t SEGMENT_DTYPE = np.dtype(
    [
        ("start", np.double),
        ("end", np.double),
        ("edge_index", np.intc),
    ],
    align=True,
)
cdef CurvedPolygonSegment[::1] SEGMENTS_WORKSPACE = np.empty(
    6, dtype=SEGMENT_DTYPE)

DQAGSE_ERR_MSGS = (
    "Maximum number of subdivisions allowed has been achieved.",
    "Roundoff error detected, which prevents convergence to tolerance.",
    'Integrand behaves "extremely" at some point(s) in the interval.',
    "Assumed: the requested tolerance cannot be achieved",
    "Integral is probably divergent or converges too slowly.",
    "Invalid input.",
)
TOO_MANY_TEMPLATE = (
    "The number of candidate intersections is too high.\n"
    "{:d} candidate pairs.")
TOO_SMALL_TEMPLATE = (
    "Did not have enough space for intersections. Needed space "
    "for {:d} intersections but only had space for {:d}.")
SEGMENT_ENDS_TOO_SMALL = (
    "Did not have enough space for segment ends. Needed space "
    "for {:d} integers but only had space for {:d}.")
SEGMENTS_TOO_SMALL = (
    "Did not have enough space for segments. Needed space "
    "for {:d} `CurvedPolygonSegment`-s but only had space for {:d}.")
# NOTE: The ``SUBDIVISION_NO_CONVERGE`` error message is copied from
#       ``_geometric_intersection.py::_NO_CONVERGE_TEMPLATE``. This
#       assumes, but does not verify, that the Fortran subroutine
#       ``curve_intersection.f90::all_intersections_abi()`` uses
#       ``MAX_INTERSECT_SUBDIVISIONS = 20``.
SUBDIVISION_NO_CONVERGE = (
    "Curve intersection failed to converge to approximately linear "
    "subdivisions after 20 iterations.")
# NOTE: The ``NEWTON_NO_CONVERGE`` error message is copied from
#       ``_intersection_helpers.py::NEWTON_NO_CONVERGE``.
NEWTON_NO_CONVERGE = """\
Unsupported multiplicity.

Newton's method failed to converge to a solution under the
following assumptions:

- The starting ``s-t`` values were already near a solution
- The root / solution has multiplicity 1 or 2
  - 1: The root is "simple", i.e. the curves are not tangent
       and have no self-intersections at the point of intersection.
  - 2: The root is a double root, i.e. the curves are tangent
       but have different curvatures at the point of intersection.

The failure to converge may have been caused by one of:

- The root was of multiplicity greater than 2
- The curves don't actually intersect, though they come very close
- Numerical issues caused the iteration to leave the region
  of convergence
"""


##########################
# Section: ``curve.f90`` #
##########################

def evaluate_multi_barycentric(
        double[::1, :] nodes, double[::1] lambda1, double[::1] lambda2):
    cdef int num_nodes, dimension, num_vals
    cdef ndarray_t[double, ndim=2, mode="fortran"] evaluated

    dimension, num_nodes = np.shape(nodes)
    num_vals, = np.shape(lambda1)
    evaluated = np.empty((dimension, num_vals), order="F")
    bezier._curve.evaluate_curve_barycentric(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &num_vals,
        &lambda1[0],
        &lambda2[0],
        &evaluated[0, 0],
    )
    return evaluated


def evaluate_multi(
        double[::1, :] nodes, double[::1] s_vals):
    cdef int num_nodes, dimension, num_vals
    cdef ndarray_t[double, ndim=2, mode="fortran"] evaluated

    dimension, num_nodes = np.shape(nodes)
    num_vals, = np.shape(s_vals)
    evaluated = np.empty((dimension, num_vals), order="F")
    bezier._curve.evaluate_multi(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &num_vals,
        &s_vals[0],
        &evaluated[0, 0],
    )
    return evaluated


def specialize_curve(double[::1, :] nodes, double start, double end):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] new_nodes

    dimension, num_nodes = np.shape(nodes)
    new_nodes = np.empty((dimension, num_nodes), order="F")

    bezier._curve.specialize_curve(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &start,
        &end,
        &new_nodes[0, 0],
    )

    return new_nodes


def evaluate_hodograph(double s, double[::1, :] nodes):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] hodograph

    dimension, num_nodes = np.shape(nodes)
    hodograph = np.empty((dimension, 1), order="F")

    bezier._curve.evaluate_hodograph(
        &s,
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &hodograph[0, 0],
    )

    return hodograph


def subdivide_nodes_curve(double[::1, :] nodes):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] left_nodes, right_nodes

    dimension, num_nodes = np.shape(nodes)

    left_nodes = np.empty((dimension, num_nodes), order="F")
    right_nodes = np.empty((dimension, num_nodes), order="F")

    bezier._curve.subdivide_nodes_curve(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &left_nodes[0, 0],
        &right_nodes[0, 0],
    )

    return left_nodes, right_nodes


def newton_refine_curve(
        double[::1, :] nodes, double[::1, :] point, double s):
    cdef int num_nodes, dimension
    cdef double updated_s

    dimension, num_nodes = np.shape(nodes)
    # NOTE: We don't check that ``np.shape(point) == (dimension, 1)``.

    bezier._curve.newton_refine_curve(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &point[0, 0],
        &s,
        &updated_s,
    )

    return updated_s


def locate_point_curve(double[::1, :] nodes, double[::1, :] point):
    cdef int num_nodes, dimension
    cdef double s_approx

    dimension, num_nodes = np.shape(nodes)
    # NOTE: We don't check that ``np.shape(point) == (dimension, 1)``.

    bezier._curve.locate_point_curve(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &point[0, 0],
        &s_approx,
    )

    if s_approx == LOCATE_MISS:
        return None
    elif s_approx == LOCATE_INVALID:
        raise ValueError(
            "Parameters not close enough to one another")
    else:
        return s_approx


def elevate_nodes(double[::1, :] nodes):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] elevated

    dimension, num_nodes = np.shape(nodes)
    elevated = np.empty((dimension, num_nodes + 1), order="F")

    bezier._curve.elevate_nodes_curve(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &elevated[0, 0],
    )

    return elevated


def get_curvature(double[::1, :] nodes, double[::1, :] tangent_vec, double s):
    cdef int num_nodes
    cdef double curvature

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes = np.shape(nodes)
    # NOTE: We don't check that ``np.shape(tangent_vec) == (2, 1)``.

    bezier._curve.get_curvature(
        &num_nodes,
        &nodes[0, 0],
        &tangent_vec[0, 0],
        &s,
        &curvature,
    )

    return curvature


def reduce_pseudo_inverse(double[::1, :] nodes):
    cdef int num_nodes, dimension
    cdef bool_t not_implemented
    cdef ndarray_t[double, ndim=2, mode="fortran"] reduced

    dimension, num_nodes = np.shape(nodes)

    reduced = np.empty((dimension, num_nodes - 1), order="F")

    bezier._curve.reduce_pseudo_inverse(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &reduced[0, 0],
        &not_implemented,
    )

    if not_implemented:
        # NOTE: This import at runtime is expensive, but we don't mind it
        #       because the exception is intended to halt the program.
        from bezier._helpers import UnsupportedDegree
        raise UnsupportedDegree(num_nodes - 1, supported=(1, 2, 3, 4))

    return reduced


def full_reduce(double[::1, :] nodes):
    cdef int num_nodes, dimension
    cdef int num_reduced_nodes
    cdef ndarray_t[double, ndim=2, mode="fortran"] reduced
    cdef bool_t not_implemented

    dimension, num_nodes = np.shape(nodes)
    reduced = np.empty((dimension, num_nodes), order="F")

    bezier._curve.full_reduce(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &num_reduced_nodes,
        &reduced[0, 0],
        &not_implemented,
    )

    if not_implemented:
        # NOTE: This import at runtime is expensive, but we don't mind it
        #       because the exception is intended to halt the program.
        from bezier._helpers import UnsupportedDegree
        raise UnsupportedDegree(num_nodes - 1, supported=(0, 1, 2, 3, 4))

    if num_reduced_nodes == num_nodes:
        if isinstance(nodes.base, np.ndarray):
            return nodes.base
        else:
            return np.asarray(nodes)
    else:
        return reduced[:, :num_reduced_nodes]


def compute_length(double[::1, :] nodes):
    cdef int num_nodes, dimension
    cdef double length
    cdef int error_val

    dimension, num_nodes = np.shape(nodes)

    bezier._curve.compute_length(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &length,
        &error_val,
    )

    if error_val == 6:
        err_msg = DQAGSE_ERR_MSGS[5]
        raise ValueError(err_msg)
    elif error_val != 0:
        if 0 <= error_val - 1 < len(DQAGSE_ERR_MSGS):
            err_msg = DQAGSE_ERR_MSGS[error_val - 1]
        else:
            err_msg = "Unknown error: {!r}.".format(error_val)
        warnings.warn(err_msg, UserWarning)

    return length

#######################################
# Section: ``curve_intersection.f90`` #
#######################################

def newton_refine_curve_intersect(
        double s, double[::1, :] nodes1, double t, double[::1, :] nodes2):
    cdef int num_nodes1, num_nodes2
    cdef double new_s, new_t
    cdef bezier._status.Status status

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes1 = np.shape(nodes1)
    # NOTE: We don't check that there are 2 rows.
    _, num_nodes2 = np.shape(nodes2)

    bezier._curve_intersection.newton_refine_curve_intersect(
        &s,
        &num_nodes1,
        &nodes1[0, 0],
        &t,
        &num_nodes2,
        &nodes2[0, 0],
        &new_s,
        &new_t,
        &status,
    )

    if status == bezier._status.Status.SINGULAR:
        raise ValueError("Jacobian is singular.")

    return new_s, new_t


def bbox_intersect(double[::1, :] nodes1, double[::1, :] nodes2):
    cdef int num_nodes1, num_nodes2
    cdef BoxIntersectionType enum_val

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes1 = np.shape(nodes1)
    # NOTE: We don't check that there are 2 rows.
    _, num_nodes2 = np.shape(nodes2)

    bezier._curve_intersection.bbox_intersect(
        &num_nodes1,
        &nodes1[0, 0],
        &num_nodes2,
        &nodes2[0, 0],
        &enum_val,
    )

    return enum_val


def reset_curves_workspace(int workspace_size):
    global CURVES_WORKSPACE
    CURVES_WORKSPACE = np.empty((2, workspace_size), order="F")


def curves_workspace_size():
    global CURVES_WORKSPACE
    cdef int intersections_size

    # NOTE: We don't check that there are 2 rows.
    _, intersections_size = np.shape(CURVES_WORKSPACE)
    return intersections_size


def curve_intersections(
        double[::1, :] nodes_first, double[::1, :] nodes_second,
        bint allow_resize=True):
    global CURVES_WORKSPACE
    cdef int num_nodes_first, num_nodes_second
    cdef int intersections_size, num_intersections
    cdef bezier._status.Status status
    cdef ndarray_t[double, ndim=2, mode="fortran"] intersections
    cdef bool_t coincident

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes_first = np.shape(nodes_first)
    _, num_nodes_second = np.shape(nodes_second)
    # NOTE: We don't check that there are 2 rows.
    _, intersections_size = np.shape(CURVES_WORKSPACE)

    bezier._curve_intersection.curve_intersections(
        &num_nodes_first,
        &nodes_first[0, 0],
        &num_nodes_second,
        &nodes_second[0, 0],
        &intersections_size,
        &CURVES_WORKSPACE[0, 0],
        &num_intersections,
        &coincident,
        &status,
    )

    if status == bezier._status.Status.SUCCESS:
        intersections = np.empty((2, num_intersections), order="F")
        intersections[:, :] = CURVES_WORKSPACE[:, :num_intersections]
        return intersections, coincident
    elif status == bezier._status.Status.NO_CONVERGE:
        raise ValueError(SUBDIVISION_NO_CONVERGE)
    elif status == bezier._status.Status.INSUFFICIENT_SPACE:
        if allow_resize:
            reset_curves_workspace(num_intersections)
            return curve_intersections(
                nodes_first, nodes_second, allow_resize=False)
        else:
            msg = TOO_SMALL_TEMPLATE.format(
                num_intersections, intersections_size)
            raise ValueError(msg)
    elif status == bezier._status.Status.BAD_MULTIPLICITY:
        raise NotImplementedError(NEWTON_NO_CONVERGE)
    else:
        # NOTE: If ``status`` isn't one of the enum values, then it is the
        #       number of candidate intersections.
        raise NotImplementedError(TOO_MANY_TEMPLATE.format(status))


def free_curve_intersections_workspace():
    bezier._curve_intersection.free_curve_intersections_workspace()

############################
# Section: ``helpers.f90`` #
############################

def cross_product(double[::1] vec0, double[::1] vec1):
    cdef double result

    bezier._helpers.cross_product(
        &vec0[0],
        &vec1[0],
        &result,
    )

    return result


def bbox(double[::1, :] nodes):
    cdef int num_nodes
    cdef double left, right, bottom, top

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes = np.shape(nodes)

    bezier._helpers.bbox(
        &num_nodes,
        &nodes[0, 0],
        &left,
        &right,
        &bottom,
        &top,
    )

    return left, right, bottom, top


def wiggle_interval(double value):
    cdef double result
    cdef bool_t success

    bezier._helpers.wiggle_interval(
        &value,
        &result,
        &success,
    )

    return result, success


def contains_nd(double[::1, :] nodes, double[::1] point):
    cdef int num_nodes, dimension
    cdef bool_t predicate

    dimension, num_nodes = np.shape(nodes)
    if np.shape(point) != (dimension,):
        msg = "Point {} was expected to have shape ({},)".format(
            np.asarray(point), dimension)
        raise ValueError(msg)

    bezier._helpers.contains_nd(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &point[0],
        &predicate,
    )

    return predicate


def vector_close(double[::1] vec1, double[::1] vec2, double eps=EPS):
    cdef int num_values

    # NOTE: We don't check that ``np.shape(vec1) == np.shape(vec2)``.
    num_values, = np.shape(vec1)

    return bezier._helpers.vector_close(
        &num_values,
        &vec1[0],
        &vec2[0],
        &eps,
    )


def in_interval(double value, double start, double end):
    return bezier._helpers.in_interval(
        &value,
        &start,
        &end,
    )


def simple_convex_hull(double[::1, :] points):
    cdef int num_points
    cdef int polygon_size
    cdef ndarray_t[double, ndim=2, mode="fortran"] polygon

    # NOTE: We don't check that there are 2 rows.
    _, num_points = np.shape(points)
    polygon = np.empty((2, num_points), order="F")

    bezier._helpers.simple_convex_hull(
        &num_points,
        &points[0, 0],
        &polygon_size,
        &polygon[0, 0],
    )

    return polygon[:, :polygon_size]


def polygon_collide(double[::1, :] polygon1, double[::1, :] polygon2):
    cdef int polygon_size1, polygon_size2
    cdef bool_t collision

    # NOTE: We don't check that there are 2 rows.
    _, polygon_size1 = np.shape(polygon1)
    _, polygon_size2 = np.shape(polygon2)

    bezier._helpers.polygon_collide(
        &polygon_size1,
        &polygon1[0, 0],
        &polygon_size2,
        &polygon2[0, 0],
        &collision,
    )

    return collision

############################
# Section: ``surface.f90`` #
############################

def de_casteljau_one_round(
        double[::1, :] nodes, int degree,
        double lambda1, double lambda2, double lambda3):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] new_nodes

    dimension, num_nodes = np.shape(nodes)
    new_nodes = np.empty((dimension, num_nodes - degree - 1), order="F")

    bezier._surface.de_casteljau_one_round(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &lambda1,
        &lambda2,
        &lambda3,
        &new_nodes[0, 0],
    )

    return new_nodes


def evaluate_barycentric(
        double[::1, :] nodes, int degree,
        double lambda1, double lambda2, double lambda3):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] point

    dimension, num_nodes = np.shape(nodes)
    point = np.empty((dimension, 1), order="F")

    bezier._surface.evaluate_barycentric(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &lambda1,
        &lambda2,
        &lambda3,
        &point[0, 0],
    )

    return point


def evaluate_barycentric_multi(
        double[::1, :] nodes, int degree,
        double[::1, :] param_vals, int dimension):
    cdef int num_nodes, num_vals
    cdef ndarray_t[double, ndim=2, mode="fortran"] evaluated

    # NOTE: We don't check that there are ``dimension`` rows.
    _, num_nodes = np.shape(nodes)
    # NOTE: We don't check that there are 3 columns.
    num_vals, _ = np.shape(param_vals)
    evaluated = np.empty((dimension, num_vals), order="F")

    bezier._surface.evaluate_barycentric_multi(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &num_vals,
        &param_vals[0, 0],
        &evaluated[0, 0],
    )

    return evaluated


def evaluate_cartesian_multi(
        double[::1, :] nodes, int degree,
        double[::1, :] param_vals, int dimension):
    cdef int num_nodes, num_vals
    cdef ndarray_t[double, ndim=2, mode="fortran"] evaluated

    # NOTE: We don't check that there are ``dimension`` rows.
    _, num_nodes = np.shape(nodes)
    # NOTE: We don't check that there are 2 columns.
    num_vals, _ = np.shape(param_vals)
    evaluated = np.empty((dimension, num_vals), order="F")

    bezier._surface.evaluate_cartesian_multi(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &num_vals,
        &param_vals[0, 0],
        &evaluated[0, 0],
    )

    return evaluated


def jacobian_both(double[::1, :] nodes, int degree, int dimension):
    cdef int num_nodes
    cdef ndarray_t[double, ndim=2, mode="fortran"] new_nodes

    # NOTE: We don't check that there are ``dimension`` rows.
    _, num_nodes = np.shape(nodes)
    new_nodes = np.empty((2 * dimension, num_nodes - degree - 1), order="F")

    bezier._surface.jacobian_both(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &new_nodes[0, 0],
    )

    return new_nodes


def jacobian_det(double[::1, :] nodes, int degree, double[::1, :] st_vals):
    cdef int num_nodes, num_vals
    cdef ndarray_t[double, ndim=1, mode="fortran"] evaluated

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes = np.shape(nodes)
    # NOTE: We don't check that there are 2 columns.
    num_vals, _ = np.shape(st_vals)

    evaluated = np.empty((num_vals,), order="F")

    bezier._surface.jacobian_det(
        &num_nodes,
        &nodes[0, 0],
        &degree,
        &num_vals,
        &st_vals[0, 0],
        &evaluated[0],
    )

    return evaluated


def specialize_surface(
        double[::1, :] nodes, int degree,
        double[::1] weights_a, double[::1] weights_b, double[::1] weights_c):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] specialized

    dimension, num_nodes = np.shape(nodes)
    specialized = np.empty((dimension, num_nodes), order="F")

    bezier._surface.specialize_surface(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &weights_a[0],
        &weights_b[0],
        &weights_c[0],
        &specialized[0, 0],
    )

    return specialized


def subdivide_nodes_surface(double[::1, :] nodes, int degree):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] nodes_a
    cdef ndarray_t[double, ndim=2, mode="fortran"] nodes_b
    cdef ndarray_t[double, ndim=2, mode="fortran"] nodes_c
    cdef ndarray_t[double, ndim=2, mode="fortran"] nodes_d

    dimension, num_nodes = np.shape(nodes)
    nodes_a = np.empty((dimension, num_nodes), order="F")
    nodes_b = np.empty((dimension, num_nodes), order="F")
    nodes_c = np.empty((dimension, num_nodes), order="F")
    nodes_d = np.empty((dimension, num_nodes), order="F")

    bezier._surface.subdivide_nodes_surface(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &nodes_a[0, 0],
        &nodes_b[0, 0],
        &nodes_c[0, 0],
        &nodes_d[0, 0],
    )

    return nodes_a, nodes_b, nodes_c, nodes_d


def compute_edge_nodes(double[::1, :] nodes, int degree):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] nodes1
    cdef ndarray_t[double, ndim=2, mode="fortran"] nodes2
    cdef ndarray_t[double, ndim=2, mode="fortran"] nodes3

    dimension, num_nodes = np.shape(nodes)
    nodes1 = np.empty((dimension, degree + 1), order="F")
    nodes2 = np.empty((dimension, degree + 1), order="F")
    nodes3 = np.empty((dimension, degree + 1), order="F")

    bezier._surface.compute_edge_nodes(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &nodes1[0, 0],
        &nodes2[0, 0],
        &nodes3[0, 0],
    )

    return nodes1, nodes2, nodes3


def compute_area(tuple edges):
    cdef int num_edges
    cdef int[::1] sizes
    cdef double** nodes_pointers
    cdef int i
    cdef double[::1, :] edge_nodes
    cdef int num_nodes
    cdef double area
    cdef bool_t unused_not_implemented

    # First, create the workspace that will convey our data.
    num_edges = len(edges)
    sizes = np.empty(num_edges, dtype=np.intc)
    nodes_pointers = <double **>malloc(num_edges * sizeof(double *))

    # Then, populate the size and pointer information.
    for i, edge_nodes in enumerate(edges):
        # NOTE: We don't check that there are 2 rows.
        _, num_nodes = np.shape(edge_nodes)

        if num_nodes - 1 > 4:
            free(nodes_pointers)

            # NOTE: This import at runtime is expensive, but we don't mind it
            #       because the exception is intended to halt the program.
            from bezier._helpers import UnsupportedDegree
            raise UnsupportedDegree(num_nodes - 1, supported=(1, 2, 3, 4))

        sizes[i] = num_nodes
        nodes_pointers[i] = &edge_nodes[0, 0]

    # Pass along the pointers to the ABI (i.e. the Fortran layer).
    # This assumes that ``unused_not_implemented`` will be ``False``
    # since we already check the supported degrees above.
    bezier._surface.compute_area(
        &num_edges,
        &sizes[0],
        nodes_pointers,
        &area,
        &unused_not_implemented,
    )

    free(nodes_pointers)

    return area

#########################################
# Section: ``surface_intersection.f90`` #
#########################################

def newton_refine_surface(
        double[::1, :] nodes, int degree,
        double x_val, double y_val, double s, double t):
    cdef int num_nodes
    cdef double updated_s, updated_t

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes = np.shape(nodes)

    bezier._surface_intersection.newton_refine_surface(
        &num_nodes,
        &nodes[0, 0],
        &degree,
        &x_val,
        &y_val,
        &s,
        &t,
        &updated_s,
        &updated_t,
    )

    return updated_s, updated_t


def locate_point_surface(
        double[::1, :] nodes, int degree, double x_val, double y_val):
    cdef int num_nodes
    cdef double s_val, t_val

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes = np.shape(nodes)

    bezier._surface_intersection.locate_point_surface(
        &num_nodes,
        &nodes[0, 0],
        &degree,
        &x_val,
        &y_val,
        &s_val,
        &t_val,
    )

    if s_val == LOCATE_MISS:
        return None
    else:
        return s_val, t_val


def reset_surface_workspaces(int segment_ends_size=-1, int segments_size=-1):
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    if segment_ends_size != -1:
        SEGMENT_ENDS_WORKSPACE = np.empty(segment_ends_size, dtype=np.intc)
    if segments_size != -1:
        SEGMENTS_WORKSPACE = np.empty(segments_size, dtype=SEGMENT_DTYPE)


def surface_workspace_sizes():
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    cdef int segment_ends_size
    cdef int segments_size

    segment_ends_size, = np.shape(SEGMENT_ENDS_WORKSPACE)
    segments_size, = np.shape(SEGMENTS_WORKSPACE)
    return segment_ends_size, segments_size


def _surface_intersections_success(
        double[::1, :] nodes1, int degree1,
        double[::1, :] nodes2, int degree2,
        int num_intersected):
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    cdef int begin_index, end_index
    cdef size_t i, j
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode="fortran"] edge_nodes1
    cdef ndarray_t[double, ndim=2, mode="fortran"] edge_nodes2
    cdef ndarray_t[double, ndim=2, mode="fortran"] edge_nodes3
    cdef ndarray_t[double, ndim=2, mode="fortran"] edge_nodes4
    cdef ndarray_t[double, ndim=2, mode="fortran"] edge_nodes5
    cdef ndarray_t[double, ndim=2, mode="fortran"] edge_nodes6

    curved_polygons = []
    for i in range(num_intersected):
        if i == 0:
            begin_index = 0
        else:
            begin_index = SEGMENT_ENDS_WORKSPACE[i - 1]
        end_index = SEGMENT_ENDS_WORKSPACE[i]
        # NOTE: We switch from 1-based to 0-based indexing since
        #       the data moves from Fortran to Python.
        triples = tuple(
            (
                SEGMENTS_WORKSPACE[j].edge_index - 1,
                SEGMENTS_WORKSPACE[j].start,
                SEGMENTS_WORKSPACE[j].end,
            )
            for j in range(begin_index, end_index)
        )
        curved_polygons.append(triples)

    if not curved_polygons:
        return [], None, ()

    # NOTE: We compute the nodes for each of the six edges. This is a
    #       "wasted" computation / storage.
    dimension, num_nodes = np.shape(nodes1)
    edge_nodes1 = np.empty((dimension, degree1 + 1), order="F")
    edge_nodes2 = np.empty((dimension, degree1 + 1), order="F")
    edge_nodes3 = np.empty((dimension, degree1 + 1), order="F")
    bezier._surface.compute_edge_nodes(
        &num_nodes,
        &dimension,
        &nodes1[0, 0],
        &degree1,
        &edge_nodes1[0, 0],
        &edge_nodes2[0, 0],
        &edge_nodes3[0, 0],
    )

    dimension, num_nodes = np.shape(nodes2)
    edge_nodes4 = np.empty((dimension, degree2 + 1), order="F")
    edge_nodes5 = np.empty((dimension, degree2 + 1), order="F")
    edge_nodes6 = np.empty((dimension, degree2 + 1), order="F")
    bezier._surface.compute_edge_nodes(
        &num_nodes,
        &dimension,
        &nodes2[0, 0],
        &degree2,
        &edge_nodes4[0, 0],
        &edge_nodes5[0, 0],
        &edge_nodes6[0, 0],
    )

    all_edge_nodes = (
        edge_nodes1,
        edge_nodes2,
        edge_nodes3,
        edge_nodes4,
        edge_nodes5,
        edge_nodes6,
    )
    return curved_polygons, None, all_edge_nodes


def _surface_intersections_resize(
        double[::1, :] nodes1, int degree1,
        double[::1, :] nodes2, int degree2,
        int segment_ends_size, int segments_size,
        int num_intersected, int resizes_allowed):
    global SEGMENT_ENDS_WORKSPACE
    cdef int num_segments

    if resizes_allowed > 0:
        if num_intersected > segment_ends_size:
            reset_surface_workspaces(segment_ends_size=num_intersected)
        else:
            num_segments = SEGMENT_ENDS_WORKSPACE[num_intersected - 1]
            reset_surface_workspaces(segments_size=num_segments)

        return surface_intersections(
            nodes1, degree1, nodes2, degree2,
            resizes_allowed=resizes_allowed - 1)
    else:
        if num_intersected > segment_ends_size:
            msg = SEGMENT_ENDS_TOO_SMALL.format(
                num_intersected, segment_ends_size)
            raise ValueError(msg)
        else:
            num_segments = SEGMENT_ENDS_WORKSPACE[num_intersected - 1]
            msg = SEGMENTS_TOO_SMALL.format(
                num_segments, segments_size)
            raise ValueError(msg)


def surface_intersections(
        double[::1, :] nodes1, int degree1,
        double[::1, :] nodes2, int degree2,
        bint verify=True, int resizes_allowed=2):
    # NOTE: ``verify`` is unused. It is only provided as compatibility
    #       with the Python version.
    global SEGMENT_ENDS_WORKSPACE
    global SEGMENTS_WORKSPACE
    cdef int num_nodes1, num_nodes2
    cdef int segment_ends_size
    cdef int segments_size
    cdef int num_intersected
    cdef SurfaceContained contained
    cdef bezier._status.Status status

    # NOTE: We don't check that there are 2 rows.
    _, num_nodes1 = np.shape(nodes1)
    _, num_nodes2 = np.shape(nodes2)

    segment_ends_size, segments_size = surface_workspace_sizes()

    bezier._surface_intersection.surface_intersections(
        &num_nodes1,
        &nodes1[0, 0],
        &degree1,
        &num_nodes2,
        &nodes2[0, 0],
        &degree2,
        &segment_ends_size,
        &SEGMENT_ENDS_WORKSPACE[0],
        &segments_size,
        &SEGMENTS_WORKSPACE[0],
        &num_intersected,
        &contained,
        &status,
    )

    if status == bezier._status.Status.SUCCESS:
        if contained == bezier._surface_intersection.SurfaceContained.FIRST:
            return None, True, ()
        elif contained == bezier._surface_intersection.SurfaceContained.SECOND:
            return None, False, ()
        else:
            # Assumes, but does not check, that ``contained`` is equal to
            # ``bezier._surface_intersection.NEITHER``.
            return _surface_intersections_success(
                nodes1, degree1, nodes2, degree2, num_intersected)
    elif status == bezier._status.Status.INSUFFICIENT_SPACE:
        return _surface_intersections_resize(
            nodes1, degree1, nodes2, degree2,
            segment_ends_size, segments_size,
            num_intersected, resizes_allowed)
    elif status == bezier._status.Status.NO_CONVERGE:
        raise ValueError(SUBDIVISION_NO_CONVERGE)
    elif status == bezier._status.Status.BAD_MULTIPLICITY:
        raise NotImplementedError(NEWTON_NO_CONVERGE)
    elif status == bezier._status.Status.EDGE_END:
        # NOTE: This text is identical (or should be) to the exception
        #       in the Python ``classify_intersection()``.
        raise ValueError("Intersection occurs at the end of an edge")
    elif status == bezier._status.Status.SAME_CURVATURE:
        # NOTE: This text is identical (or should be) to the exception(s)
        #       in the Python ``classify_tangent_intersection()``.
        raise NotImplementedError("Tangent curves have same curvature.")
    elif status == bezier._status.Status.BAD_INTERIOR:
        # NOTE: This assumes that the Fortran ``interior_combine()`` has
        #       failed to return to the start node after ``MAX_EDGES = 10``
        #       edges have been added.
        raise RuntimeError("Unexpected number of edges")
    elif status == bezier._status.Status.UNKNOWN:
        raise RuntimeError("Unknown error has occured.")
    else:
        # NOTE: If ``status`` isn't one of the enum values, then it is the
        #       number of candidate intersections.
        raise NotImplementedError(TOO_MANY_TEMPLATE.format(status))


def free_surface_intersections_workspace():
    bezier._surface_intersection.free_surface_intersections_workspace()


def _type_info():
    return (
        SEGMENT_DTYPE.isnative,
        SEGMENT_DTYPE.itemsize,
        SEGMENT_DTYPE.num,
        sizeof(CurvedPolygonSegment),
    )
