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

"""Pure Python helper methods for :mod:`bezier.triangle`.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import collections
import itertools

import numpy as np

from bezier.hazmat import algebraic_intersection
from bezier.hazmat import geometric_intersection
from bezier.hazmat import helpers as _py_helpers
from bezier.hazmat import intersection_helpers
from bezier.hazmat import triangle_helpers


MAX_LOCATE_SUBDIVISIONS = 20
LOCATE_EPS = 0.5**47
INTERSECTION_T = geometric_intersection.BoxIntersectionType.INTERSECTION
CLASSIFICATION_T = intersection_helpers.IntersectionClassification
UNUSED_T = CLASSIFICATION_T.COINCIDENT_UNUSED
ACCEPTABLE_CLASSIFICATIONS = (
    CLASSIFICATION_T.FIRST,
    CLASSIFICATION_T.SECOND,
    CLASSIFICATION_T.COINCIDENT,
)
TANGENT_CLASSIFICATIONS = (
    CLASSIFICATION_T.TANGENT_FIRST,
    CLASSIFICATION_T.TANGENT_SECOND,
)
BAD_SEGMENT_PARAMS = (
    "Segment start should be strictly less than end and both "
    "should be in [0, 1]."
)
SEGMENTS_SAME_EDGE = "Consecutive segments lie on the same edge."


def newton_refine_solve(jac_both, x_val, triangle_x, y_val, triangle_y):
    r"""Helper for :func:`newton_refine`.

    We have a system:

    .. code-block:: rest

       [A C][ds] = [E]
       [B D][dt]   [F]

    This is not a typo, ``A->B->C->D`` matches the Fortran-ordering of the data
    in ``jac_both``. We solve directly rather than using a linear algebra
    utility:

    .. code-block:: rest

       ds = (D E - C F) / (A D - B C)
       dt = (A F - B E) / (A D - B C)

    Args:
        jac_both (numpy.ndarray): A ``4 x 1`` matrix of entries in a Jacobian.
        x_val (float): An ``x``-value we are trying to reach.
        triangle_x (float): The actual ``x``-value we are currently at.
        y_val (float): An ``y``-value we are trying to reach.
        triangle_y (float): The actual ``y``-value we are currently at.

    Returns:
        Tuple[float, float]: The pair of values the solve the
        linear system.
    """
    a_val, b_val, c_val, d_val = jac_both[:, 0]
    #       and
    e_val = x_val - triangle_x
    f_val = y_val - triangle_y
    # Now solve:
    denom = a_val * d_val - b_val * c_val
    delta_s = (d_val * e_val - c_val * f_val) / denom
    delta_t = (a_val * f_val - b_val * e_val) / denom
    return delta_s, delta_t


def newton_refine(nodes, degree, x_val, y_val, s, t):
    r"""Refine a solution to :math:`B(s, t) = p` using Newton's method.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Computes updates via

    .. math::

       \left[\begin{array}{c}
           0 \\ 0 \end{array}\right] \approx
           \left(B\left(s_{\ast}, t_{\ast}\right) -
           \left[\begin{array}{c} x \\ y \end{array}\right]\right) +
           \left[\begin{array}{c c}
               B_s\left(s_{\ast}, t_{\ast}\right) &
               B_t\left(s_{\ast}, t_{\ast}\right) \end{array}\right]
           \left[\begin{array}{c}
               \Delta s \\ \Delta t \end{array}\right]

    For example, (with weights
    :math:`\lambda_1 = 1 - s - t, \lambda_2 = s, \lambda_3 = t`)
    consider the triangle

    .. math::

       B(s, t) =
           \left[\begin{array}{c} 0 \\ 0 \end{array}\right] \lambda_1^2 +
           \left[\begin{array}{c} 1 \\ 0 \end{array}\right]
               2 \lambda_1 \lambda_2 +
           \left[\begin{array}{c} 2 \\ 0 \end{array}\right] \lambda_2^2 +
           \left[\begin{array}{c} 2 \\ 1 \end{array}\right]
               2 \lambda_1 \lambda_3 +
           \left[\begin{array}{c} 2 \\ 2 \end{array}\right]
               2 \lambda_2 \lambda_1 +
           \left[\begin{array}{c} 0 \\ 2 \end{array}\right] \lambda_3^2

    and the point
    :math:`B\left(\frac{1}{4}, \frac{1}{2}\right) =
    \frac{1}{4} \left[\begin{array}{c} 5 \\ 5 \end{array}\right]`.

    Starting from the **wrong** point
    :math:`s = \frac{1}{2}, t = \frac{1}{4}`, we have

    .. math::

       \begin{align*}
       \left[\begin{array}{c} x \\ y \end{array}\right] -
          B\left(\frac{1}{2}, \frac{1}{4}\right) &= \frac{1}{4}
          \left[\begin{array}{c} -1 \\ 2 \end{array}\right] \\
       DB\left(\frac{1}{2}, \frac{1}{4}\right) &= \frac{1}{2}
           \left[\begin{array}{c c} 3 & 2 \\ 1 & 6 \end{array}\right] \\
       \Longrightarrow \left[\begin{array}{c}
           \Delta s \\ \Delta t \end{array}\right] &= \frac{1}{32}
           \left[\begin{array}{c} -10 \\ 7 \end{array}\right]
       \end{align*}

    .. image:: ../../images/newton_refine_triangle.png
       :align: center

    .. testsetup:: newton-refine-triangle

       from bezier.hazmat.triangle_intersection import newton_refine

    .. doctest:: newton-refine-triangle

       >>> import bezier
       >>> import numpy as np
       >>> nodes = np.asfortranarray([
       ...     [0.0, 1.0, 2.0, 2.0, 2.0, 0.0],
       ...     [0.0, 0.0, 0.0, 1.0, 2.0, 2.0],
       ... ])
       >>> triangle = bezier.Triangle(nodes, degree=2)
       >>> bool(triangle.is_valid)
       True
       >>> (x_val,), (y_val,) = triangle.evaluate_cartesian(0.25, 0.5)
       >>> float(x_val), float(y_val)
       (1.25, 1.25)
       >>> s, t = 0.5, 0.25
       >>> new_s, new_t = newton_refine(nodes, 2, x_val, y_val, s, t)
       >>> float(32 * (new_s - s))
       -10.0
       >>> float(32 * (new_t - t))
       7.0

    .. testcleanup:: newton-refine-triangle

       import make_images
       make_images.newton_refine_triangle(
           triangle, x_val, y_val, s, t, new_s, new_t)

    Args:
        nodes (numpy.ndarray): Array of nodes in a triangle.
        degree (int): The degree of the triangle.
        x_val (float): The :math:`x`-coordinate of a point
            on the triangle.
        y_val (float): The :math:`y`-coordinate of a point
            on the triangle.
        s (float): Approximate :math:`s`-value to be refined.
        t (float): Approximate :math:`t`-value to be refined.

    Returns:
        Tuple[float, float]: The refined :math:`s` and :math:`t` values.
    """
    lambda1 = 1.0 - s - t
    (triangle_x,), (triangle_y,) = triangle_helpers.evaluate_barycentric(
        nodes, degree, lambda1, s, t
    )
    if triangle_x == x_val and triangle_y == y_val:
        # No refinement is needed.
        return s, t

    # NOTE: This function assumes ``dimension==2`` (i.e. since ``x, y``).
    jac_nodes = triangle_helpers.jacobian_both(nodes, degree, 2)
    # The degree of the jacobian is one less.
    jac_both = triangle_helpers.evaluate_barycentric(
        jac_nodes, degree - 1, lambda1, s, t
    )
    # The first row of the jacobian matrix is B_s (i.e. the
    # top-most values in ``jac_both``).
    delta_s, delta_t = newton_refine_solve(
        jac_both, x_val, triangle_x, y_val, triangle_y
    )
    return s + delta_s, t + delta_t


def update_locate_candidates(candidate, next_candidates, x_val, y_val, degree):
    """Update list of candidate triangles during geometric search for a point.

    .. note::

       This is used **only** as a helper for :func:`locate_point`.

    Checks if the point ``(x_val, y_val)`` is contained in the ``candidate``
    triangle. If not, this function does nothing. If the point is contained,
    the four subdivided triangles from ``candidate`` are added to
    ``next_candidates``.

    Args:
        candidate (Tuple[float, float, float, numpy.ndarray]): A 4-tuple
            describing a triangle and its centroid / width. Contains

            * Three times centroid ``x``-value
            * Three times centroid ``y``-value
            * "Width" of parameter space for the triangle
            * Control points for the triangle
        next_candidates (list): List of "candidate" sub-triangles that may
            contain the point being located.
        x_val (float): The ``x``-coordinate being located.
        y_val (float): The ``y``-coordinate being located.
        degree (int): The degree of the triangle.
    """
    centroid_x, centroid_y, width, candidate_nodes = candidate
    point = np.asfortranarray([x_val, y_val])
    if not _py_helpers.contains_nd(candidate_nodes, point):
        return

    nodes_a, nodes_b, nodes_c, nodes_d = triangle_helpers.subdivide_nodes(
        candidate_nodes, degree
    )
    half_width = 0.5 * width
    next_candidates.extend(
        (
            (
                centroid_x - half_width,
                centroid_y - half_width,
                half_width,
                nodes_a,
            ),
            (centroid_x, centroid_y, -half_width, nodes_b),
            (centroid_x + width, centroid_y - half_width, half_width, nodes_c),
            (centroid_x - half_width, centroid_y + width, half_width, nodes_d),
        )
    )


def mean_centroid(candidates):
    """Take the mean of all centroids in set of reference triangles.

    .. note::

       This is used **only** as a helper for :func:`locate_point`.

    Args:
        candidates (List[Tuple[float, float, float, numpy.ndarray]]): List of
            4-tuples, each of which has been produced by :func:`locate_point`.
            Each 4-tuple contains

            * Three times centroid ``x``-value
            * Three times centroid ``y``-value
            * "Width" of a parameter space for a triangle
            * Control points for a triangle

            We only use the first two values, which are triple the desired
            value so that we can put off division by three until summing in
            our average. We don't use the other two values, they are just an
            artifact of the way ``candidates`` is constructed by the caller.

    Returns:
        Tuple[float, float]: The mean of all centroids.
    """
    sum_x = 0.0
    sum_y = 0.0
    for centroid_x, centroid_y, _, _ in candidates:
        sum_x += centroid_x
        sum_y += centroid_y
    denom = 3.0 * len(candidates)
    return sum_x / denom, sum_y / denom


def locate_point(nodes, degree, x_val, y_val):
    r"""Locate a point on a triangle.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Does so by recursively subdividing the triangle and rejecting
    sub-triangles with bounding boxes that don't contain the point.
    After the sub-triangles are sufficiently small, uses Newton's
    method to narrow in on the pre-image of the point.

    Args:
        nodes (numpy.ndarray): Control points for B |eacute| zier triangle
            (assumed to be two-dimensional).
        degree (int): The degree of the triangle.
        x_val (float): The :math:`x`-coordinate of a point
            on the triangle.
        y_val (float): The :math:`y`-coordinate of a point
            on the triangle.

    Returns:
        Optional[Tuple[float, float]]: The :math:`s` and :math:`t`
        values corresponding to ``x_val`` and ``y_val`` or
        :data:`None` if the point is not on the ``triangle``.
    """
    # We track the centroid rather than base_x/base_y/width (by storing triple
    # the centroid -- to avoid division by three until needed). We also need
    # to track the width (or rather, just the sign of the width).
    candidates = [(1.0, 1.0, 1.0, nodes)]
    for _ in range(MAX_LOCATE_SUBDIVISIONS + 1):
        next_candidates = []
        for candidate in candidates:
            update_locate_candidates(
                candidate, next_candidates, x_val, y_val, degree
            )
        candidates = next_candidates
    if not candidates:
        return None

    # We take the average of all centroids from the candidates
    # that may contain the point.
    s_approx, t_approx = mean_centroid(candidates)
    s, t = newton_refine(nodes, degree, x_val, y_val, s_approx, t_approx)
    actual = triangle_helpers.evaluate_barycentric(
        nodes, degree, 1.0 - s - t, s, t
    )
    expected = np.asfortranarray([x_val, y_val])
    if not _py_helpers.vector_close(
        actual.ravel(order="F"), expected, eps=LOCATE_EPS
    ):
        s, t = newton_refine(nodes, degree, x_val, y_val, s, t)
    return s, t


def same_intersection(intersection1, intersection2, wiggle=0.5**40):
    """Check if two intersections are close to machine precision.

    .. note::

       This is a helper used only by :func:`verify_duplicates`, which in turn
       is only used by :func:`generic_intersect`.

    Args:
        intersection1 (Intersection): The first intersection.
        intersection2 (Intersection): The second intersection.
        wiggle (Optional[float]): The amount of relative error allowed
            in parameter values.

    Returns:
        bool: Indicates if the two intersections are the same to
        machine precision.
    """
    if intersection1.index_first != intersection2.index_first:
        return False

    if intersection1.index_second != intersection2.index_second:
        return False

    return np.allclose(
        [intersection1.s, intersection1.t],
        [intersection2.s, intersection2.t],
        atol=0.0,
        rtol=wiggle,
    )


def verify_duplicates(duplicates, uniques):
    """Verify that a set of intersections had expected duplicates.

    .. note::

       This is a helper used only by :func:`generic_intersect`.

    Args:
        duplicates (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of intersections corresponding to duplicates that were
            filtered out.
        uniques (List[~bezier.hazmat.intersection_helpers.Intersection]): List
            of "final" intersections with duplicates filtered out.

    Raises:
        ValueError: If the ``uniques`` are not actually all unique.
        ValueError: If one of the ``duplicates`` does not correspond to
            an intersection in ``uniques``.
        ValueError: If a duplicate occurs only once but does not have
            exactly one of ``s`` and ``t`` equal to ``0.0``.
        ValueError: If a duplicate occurs three times but does not have
            exactly both ``s == t == 0.0``.
        ValueError: If a duplicate occurs a number other than one or three
            times.
    """
    for uniq1, uniq2 in itertools.combinations(uniques, 2):
        if same_intersection(uniq1, uniq2):
            raise ValueError("Non-unique intersection")

    counter = collections.Counter()
    for dupe in duplicates:
        matches = []
        for index, uniq in enumerate(uniques):
            if same_intersection(dupe, uniq):
                matches.append(index)
        if len(matches) != 1:
            raise ValueError("Duplicate not among uniques", dupe)

        matched = matches[0]
        counter[matched] += 1
    for index, count in counter.items():
        uniq = uniques[index]
        if count == 1:
            if (uniq.s, uniq.t).count(0.0) != 1:
                raise ValueError("Count == 1 should be a single corner", uniq)

        elif count == 3:
            if (uniq.s, uniq.t) != (0.0, 0.0):
                raise ValueError("Count == 3 should be a double corner", uniq)

        else:
            raise ValueError("Unexpected duplicate count", count)


def verify_edge_segments(edge_infos):
    """Verify that the edge segments in an intersection are valid.

    .. note::

       This is a helper used only by :func:`generic_intersect`.

    Args:
        edge_infos (Optional[list]): List of "edge info" lists. Each list
            represents a curved polygon and contains 3-tuples of edge index,
            start and end (see the output of :func:`.ends_to_curve`).

    Raises:
        ValueError: If two consecutive edge segments lie on the same edge
            index.
        ValueError: If the start and end parameter are "invalid" (they should
            be between 0 and 1 and start should be strictly less than end).
    """
    if edge_infos is None:
        return

    for edge_info in edge_infos:
        num_segments = len(edge_info)
        for index in range(-1, num_segments - 1):
            index1, start1, end1 = edge_info[index]
            # First, verify the start and end parameters for the current
            # segment.
            if not 0.0 <= start1 < end1 <= 1.0:
                raise ValueError(BAD_SEGMENT_PARAMS, edge_info[index])

            # Then, verify that the indices are not the same.
            index2, _, _ = edge_info[index + 1]
            if index1 == index2:
                raise ValueError(
                    SEGMENTS_SAME_EDGE, edge_info[index], edge_info[index + 1]
                )


def add_edge_end_unused(intersection, duplicates, intersections):
    """Add intersection that is classified as "unused" and on an edge end.

    This is a helper for
    :func:`~.hazmat.triangle_intersection.add_intersection` for intersections
    classified as :attr:`~.IntersectionClassification.COINCIDENT_UNUSED`. It
    assumes that

    * ``intersection`` will have at least one of ``s == 0.0`` or ``t == 0.0``
    * A "misclassified" intersection in ``intersections`` that matches
      ``intersection`` will be the "same" if it matches both ``index_first``
      and ``index_second`` and if it matches the start index exactly

    Args:
        intersection (Intersection): An intersection to be added.
        duplicates (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of duplicate intersections.
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of "accepted" (i.e. non-duplicate) intersections.
    """
    found = None
    for other in intersections:
        if (
            intersection.index_first == other.index_first
            and intersection.index_second == other.index_second
        ):
            if intersection.s == 0.0 and other.s == 0.0:
                found = other
                break

            if intersection.t == 0.0 and other.t == 0.0:
                found = other
                break

    if found is not None:
        intersections.remove(found)
        duplicates.append(found)
    intersections.append(intersection)


def check_unused(intersection, duplicates, intersections):
    """Check if a "valid" ``intersection`` is already in ``intersections``.

    This assumes that

    * ``intersection`` will have at least one of ``s == 0.0`` or ``t == 0.0``
    * At least one of the intersections in ``intersections`` is classified as
      :attr:`~.IntersectionClassification.COINCIDENT_UNUSED`.

    Args:
        intersection (Intersection): An intersection to be added.
        duplicates (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of duplicate intersections.
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of "accepted" (i.e. non-duplicate) intersections.

    Returns:
        bool: Indicates if the ``intersection`` is a duplicate.
    """
    for other in intersections:
        if (
            other.interior_curve == UNUSED_T
            and intersection.index_first == other.index_first
            and intersection.index_second == other.index_second
        ):
            if intersection.s == 0.0 and other.s == 0.0:
                duplicates.append(intersection)
                return True

            if intersection.t == 0.0 and other.t == 0.0:
                duplicates.append(intersection)
                return True

    return False


def add_intersection(  # pylint: disable=too-many-arguments
    index1,
    s,
    index2,
    t,
    interior_curve,
    edge_nodes1,
    edge_nodes2,
    duplicates,
    intersections,
):
    """Create an :class:`.Intersection` and append.

    The intersection will be classified as either a duplicate or a valid
    intersection and appended to one of ``duplicates`` or ``intersections``
    depending on that classification.

    Args:
        index1 (int): The index (among 0, 1, 2) of the first edge in the
            intersection.
        s (float): The parameter along the first curve of the intersection.
        index2 (int): The index (among 0, 1, 2) of the second edge in the
            intersection.
        t (float): The parameter along the second curve of the intersection.
        interior_curve (Optional[ \
            ~bezier.hazmat.intersection_helpers.IntersectionClassification]):
            The classification of the intersection, if known. If :data:`None`,
            the classification will be computed.
        edge_nodes1 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the first triangle being intersected.
        edge_nodes2 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the second triangle being intersected.
        duplicates (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of duplicate intersections.
        intersections (List[~bezier.hazmat.intersection_helpers.Intersection]):
            List of "accepted" (i.e. non-duplicate) intersections.
    """
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-locals
    edge_end, is_corner, intersection_args = triangle_helpers.handle_ends(
        index1, s, index2, t
    )
    if edge_end:
        intersection = intersection_helpers.Intersection(*intersection_args)
        intersection.interior_curve = interior_curve
        if interior_curve == UNUSED_T:
            add_edge_end_unused(intersection, duplicates, intersections)
        else:
            duplicates.append(intersection)
    else:
        intersection = intersection_helpers.Intersection(index1, s, index2, t)
        if is_corner:
            is_duplicate = check_unused(
                intersection, duplicates, intersections
            )
            if is_duplicate:
                return

        # Classify the intersection.
        if interior_curve is None:
            interior_curve = triangle_helpers.classify_intersection(
                intersection, edge_nodes1, edge_nodes2
            )
        intersection.interior_curve = interior_curve
        intersections.append(intersection)


def classify_coincident(st_vals, coincident):
    """Determine if coincident parameters are "unused".

    .. note::

       This is a helper for :func:`triangle_intersections`.

    In the case that ``coincident`` is :data:`True`, then we'll have two
    sets of parameters :math:`(s_1, t_1)` and :math:`(s_2, t_2)`.

    If one of :math:`s_1 < s_2` or :math:`t_1 < t_2` is not satisfied, the
    coincident segments will be moving in opposite directions, hence don't
    define an interior of an intersection.

    .. warning::

       In the "coincident" case, this assumes, but doesn't check, that
       ``st_vals`` is ``2 x 2``.

    Args:
        st_vals (numpy.ndarray): ``2 X N`` array of intersection parameters.
        coincident (bool): Flag indicating if the intersections are the
            endpoints of coincident segments of two curves.

    Returns:
        Optional[ \
        ~bezier.hazmat.intersection_helpers.IntersectionClassification]:
        The classification of the intersections.
    """
    if not coincident:
        return None

    if st_vals[0, 0] >= st_vals[0, 1] or st_vals[1, 0] >= st_vals[1, 1]:
        return UNUSED_T

    else:
        return CLASSIFICATION_T.COINCIDENT


def should_use(intersection):
    """Check if an intersection can be used as part of a curved polygon.

    Will return :data:`True` if the intersection is classified as
    :attr:`~.IntersectionClassification.FIRST`,
    :attr:`~.IntersectionClassification.SECOND` or
    :attr:`~.IntersectionClassification.COINCIDENT` or if the intersection
    is classified is a corner / edge end which is classified as
    :attr:`~.IntersectionClassification.TANGENT_FIRST` or
    :attr:`~.IntersectionClassification.TANGENT_SECOND`.

    Args:
        intersection (Intersection): An intersection to be added.

    Returns:
        bool: Indicating if the intersection will be used.
    """
    if intersection.interior_curve in ACCEPTABLE_CLASSIFICATIONS:
        return True

    if intersection.interior_curve in TANGENT_CLASSIFICATIONS:
        return intersection.s == 0.0 or intersection.t == 0.0

    return False


def triangle_intersections(edge_nodes1, edge_nodes2, all_intersections):
    """Find all intersections among edges of two triangles.

    This treats intersections which have ``s == 1.0`` or ``t == 1.0``
    as duplicates. The duplicates may be checked by the caller, e.g.
    by :func:`verify_duplicates`.

    Args:
        edge_nodes1 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the first triangle being intersected.
        edge_nodes2 (Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]): The
            nodes of the three edges of the second triangle being intersected.
        all_intersections (Callable[[numpy.ndarray, numpy.ndarray], \
            Tuple[numpy.ndarray, bool]]): A helper that intersects
            B |eacute| zier curves. Takes the nodes of each curve as input and
            returns an array (``2 x N``) of intersections.

    Returns:
        Tuple[list, list, list, set]: 4-tuple of

        * The actual "unique" :class:`.Intersection`-s
        * Duplicate :class:`.Intersection`-s encountered (these will be
          corner intersections)
        * Intersections that won't be used, such as a tangent intersection
          along an edge
        * All the intersection classifications encountered
    """
    # NOTE: There is no corresponding "enable", but the disable only applies
    #       in this lexical scope.
    # pylint: disable=too-many-locals
    intersections = []
    duplicates = []
    for index1, nodes1 in enumerate(edge_nodes1):
        for index2, nodes2 in enumerate(edge_nodes2):
            st_vals, coincident = all_intersections(nodes1, nodes2)
            interior_curve = classify_coincident(st_vals, coincident)
            for s, t in st_vals.T:
                add_intersection(
                    index1,
                    s,
                    index2,
                    t,
                    interior_curve,
                    edge_nodes1,
                    edge_nodes2,
                    duplicates,
                    intersections,
                )
    all_types = set()
    to_keep = []
    unused = []
    for intersection in intersections:
        all_types.add(intersection.interior_curve)
        # Only keep the intersections which are "acceptable".
        if should_use(intersection):
            to_keep.append(intersection)
        else:
            unused.append(intersection)
    return to_keep, duplicates, unused, all_types


def generic_intersect(
    nodes1, degree1, nodes2, degree2, verify, all_intersections
):
    """Find all intersections among edges of two triangles.

    This treats intersections which have ``s == 1.0`` or ``t == 1.0``
    as duplicates. The duplicates will be checked by :func:`verify_duplicates`
    if ``verify`` is :data:`True`.

    Args:
        nodes1 (numpy.ndarray): The nodes defining the first triangle in
            the intersection (assumed in :math:`\\mathbf{R}^2`).
        degree1 (int): The degree of the triangle given by ``nodes1``.
        nodes2 (numpy.ndarray): The nodes defining the second triangle in
            the intersection (assumed in :math:`\\mathbf{R}^2`).
        degree2 (int): The degree of the triangle given by ``nodes2``.
        verify (Optional[bool]): Indicates if duplicate intersections
            should be checked.
        all_intersections (Callable[[numpy.ndarray, numpy.ndarray], \
            Tuple[numpy.ndarray, bool]]): A helper that intersects
            B |eacute| zier curves. Takes the nodes of each curve as input and
            returns an array (``2 x N``) of intersections.

    Returns:
        Tuple[Optional[list], Optional[bool], tuple]: 3-tuple of

        * List of "edge info" lists. Each list represents a curved polygon
          and contains 3-tuples of edge index, start and end (see the
          output of :func:`.ends_to_curve`).
        * "Contained" boolean. If not :data:`None`, indicates
          that one of the triangles is contained in the other.
        * The nodes of three edges of the first triangle being intersected
          followed by the nodes of the three edges of the second.
    """
    bbox_int = geometric_intersection.bbox_intersect(nodes1, nodes2)
    if bbox_int != INTERSECTION_T:
        return [], None, ()

    # We need **all** edges (as nodes).
    edge_nodes1 = triangle_helpers.compute_edge_nodes(nodes1, degree1)
    edge_nodes2 = triangle_helpers.compute_edge_nodes(nodes2, degree2)
    # Run through **all** pairs of edges.
    intersections, duplicates, unused, all_types = triangle_intersections(
        edge_nodes1, edge_nodes2, all_intersections
    )
    edge_infos, contained = triangle_helpers.combine_intersections(
        intersections, nodes1, degree1, nodes2, degree2, all_types
    )
    # Verify the duplicates and intersection if need be.
    if verify:
        verify_duplicates(duplicates, intersections + unused)
        verify_edge_segments(edge_infos)
    if edge_infos is None or edge_infos == []:
        return edge_infos, contained, ()

    return edge_infos, contained, edge_nodes1 + edge_nodes2


def geometric_intersect(nodes1, degree1, nodes2, degree2, verify):
    r"""Find all intersections among edges of two triangles.

    .. note::

       There is also a Fortran implementation of this function, which
       will be used if it can be built.

    Uses :func:`generic_intersect` with the
    :attr:`~.IntersectionStrategy.GEOMETRIC` intersection strategy.

    Args:
        nodes1 (numpy.ndarray): The nodes defining the first triangle in
            the intersection (assumed in :math:`\mathbf{R}^2`).
        degree1 (int): The degree of the triangle given by ``nodes1``.
        nodes2 (numpy.ndarray): The nodes defining the second triangle in
            the intersection (assumed in :math:`\mathbf{R}^2`).
        degree2 (int): The degree of the triangle given by ``nodes2``.
        verify (Optional[bool]): Indicates if duplicate intersections
            should be checked.

    Returns:
        Tuple[Optional[list], Optional[bool], tuple]: 3-tuple of

        * List of "edge info" lists. Each list represents a curved polygon
          and contains 3-tuples of edge index, start and end (see the
          output of :func:`.ends_to_curve`).
        * "Contained" boolean. If not :data:`None`, indicates
          that one of the triangles is contained in the other.
        * The nodes of three edges of the first triangle being intersected
          followed by the nodes of the three edges of the second.
    """
    all_intersections = geometric_intersection.all_intersections
    return generic_intersect(
        nodes1, degree1, nodes2, degree2, verify, all_intersections
    )


def algebraic_intersect(nodes1, degree1, nodes2, degree2, verify):
    r"""Find all intersections among edges of two triangles.

    Uses :func:`generic_intersect` with the
    :attr:`~.IntersectionStrategy.ALGEBRAIC` intersection strategy.

    Args:
        nodes1 (numpy.ndarray): The nodes defining the first triangle in
            the intersection (assumed in :math:`\mathbf{R}^2`).
        degree1 (int): The degree of the triangle given by ``nodes1``.
        nodes2 (numpy.ndarray): The nodes defining the second triangle in
            the intersection (assumed in :math:`\mathbf{R}^2`).
        degree2 (int): The degree of the triangle given by ``nodes2``.
        verify (Optional[bool]): Indicates if duplicate intersections
            should be checked.

    Returns:
        Tuple[Optional[list], Optional[bool], tuple]: 3-tuple of

        * List of "edge info" lists. Each list represents a curved polygon
          and contains 3-tuples of edge index, start and end (see the
          output of :func:`.ends_to_curve`).
        * "Contained" boolean. If not :data:`None`, indicates
          that one of the triangles is contained in the other.
        * The nodes of three edges of the first triangle being intersected
          followed by the nodes of the three edges of the second.
    """
    all_intersections = algebraic_intersection.all_intersections
    return generic_intersect(
        nodes1, degree1, nodes2, degree2, verify, all_intersections
    )
