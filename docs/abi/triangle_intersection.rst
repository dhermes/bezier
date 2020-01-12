############################
triangle_intersection module
############################

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

This is a collection of procedures and types for computing intersections
between two B |eacute| zier triangles in :math:`\mathbf{R}^2`. The region(s)
of intersection (if non-empty) will be curved polygons defined by
the B |eacute| zier curve segments that form the exterior.

.. note::

    In most of the procedures below both the number of nodes :math:`N` and
    the degree :math:`d` of a B |eacute| zier triangle are provided. It is
    redundant to require both as arguments since :math:`N = \binom{d + 2}{2}`.
    However, both are provided as arguments to avoid unnecessary
    re-computation, i.e. we expect the caller to know both :math:`N` and
    :math:`d`.

**********
Procedures
**********

.. c:function:: void BEZ_locate_point_triangle(const int *num_nodes, \
                                              const double *nodes, \
                                              const int *degree, \
                                              const double *x_val, \
                                              const double *y_val, \
                                              double *s_val, \
                                              double *t_val)

   This solves the inverse problem :math:`B(s, t) = (x, y)` (if it can be
   solved). Does so by subdividing the triangle until the sub-triangles are
   sufficiently small, then using Newton's method to narrow in on the
   pre-image of the point.

   This assumes the triangle is "valid", i.e. it has positive Jacobian
   throughout the unit triangle. (If this were not true, then multiple
   inputs could map to the same output.)

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: const int*
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`2 \times N` array. This should be laid out in Fortran order, with
      :math:`2 N` total values.
   :type nodes: const double*
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: const int*
   :param x_val:
      **[Input]** The :math:`x`\-value of the point being located.
   :type x_val: const double*
   :param y_val:
      **[Input]** The :math:`y`\-value of the point being located.
   :type y_val: const double*
   :param double* s_val:
      **[Output]** The first parameter :math:`s` of the solution. If
      :math:`(x, y)` can't be located on the triangle, then ``s_val = -1.0``.
   :param double* t_val:
      **[Output]** The second parameter :math:`t` of the solution. If
      :math:`(x, y)` can't be located on the triangle, then this value
      is undefined.

   **Signature:**

   .. code-block:: c

      void
      BEZ_locate_point_triangle(const int *num_nodes,
                               const double *nodes,
                               const int *degree,
                               const double *x_val,
                               const double *y_val,
                               double *s_val,
                               double *t_val);

.. c:function:: void BEZ_newton_refine_triangle(const int *num_nodes, \
                                               const double *nodes, \
                                               const int *degree, \
                                               const double *x_val, \
                                               const double *y_val, \
                                               const double *s, \
                                               const double *t, \
                                               double *updated_s, \
                                               double *updated_t)

   This refines a solution to :math:`B(s, t) = (x, y) = p` using Newton's
   method. Given a current approximation :math:`(s_n, t_n)` for a solution,
   this produces the updated approximation via

   .. math::

      \left[\begin{array}{c} s_{n + 1} \\ t_{n + 1} \end{array}\right] =
      \left[\begin{array}{c} s_n \\ t_n \end{array}\right] -
      DB(s_n, t_n)^{-1} \left(B(s_n, t_n) - p\right).

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: const int*
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`2 \times N` array. This should be laid out in Fortran order, with
      :math:`2 N` total values.
   :type nodes: const double*
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: const int*
   :param x_val:
      **[Input]** The :math:`x`\-value of the point :math:`p`.
   :type x_val: const double*
   :param y_val:
      **[Input]** The :math:`y`\-value of the point :math:`p`.
   :type y_val: const double*
   :param s:
      **[Input]** The first parameter :math:`s_n` of the current approximation
      of a solution.
   :type s: const double*
   :param t:
      **[Input]** The second parameter :math:`t_n` of the current approximation
      of a solution.
   :type t: const double*
   :param double* updated_s:
      **[Output]** The first parameter :math:`s_{n + 1}` of the updated
      approximation.
   :param double* updated_t:
      **[Output]** The second parameter :math:`t_{n + 1}` of the updated
      approximation.

   **Signature:**

   .. code-block:: c

      void
      BEZ_newton_refine_triangle(const int *num_nodes,
                                const double *nodes,
                                const int *degree,
                                const double *x_val,
                                const double *y_val,
                                const double *s,
                                const double *t,
                                double *updated_s,
                                double *updated_t);

.. c:function:: void BEZ_triangle_intersections(const int *num_nodes1, \
                                               const double *nodes1, \
                                               const int *degree1, \
                                               const int *num_nodes2, \
                                               const double *nodes2, \
                                               const int *degree2, \
                                               const int *segment_ends_size, \
                                               int *segment_ends, \
                                               const int *segments_size, \
                                               CurvedPolygonSegment *segments, \
                                               int *num_intersected, \
                                               TriangleContained *contained, \
                                               Status *status)

   Compute the intersection of two B |eacute| zier triangles. This will
   first compute all intersection points between edges of the first and
   second triangle (nine edge pairs in total). Then, it will classify each
   point according to which triangle is "interior" at that point. Finally,
   it will form a loop of intersection points using the classifications
   until all intersections have been used or discarded.

   .. tip::

      If the ``status`` returned is :c:data:`INSUFFICIENT_SPACE` that means
      either

      * ``segment_ends_size`` is smaller than ``num_intersected``
        so ``segment_ends`` needs to be resized to at least as large as
        ``num_intersected``.
      * ``segments_size`` is smaller than the number of segments. The number
        of segments will be the last index in the list of edge indices:
        ``segment_ends[num_intersected - 1]``. In this case ``segments``
        needs to be resized.

      This means a successful invocation of this procedure may take three
      attempts. To avoid false starts occurring on a regular basis, keep a
      static workspace around that will continue to grow as resizing is
      needed, but will never shrink.

   :param num_nodes1:
      **[Input]** The number of nodes :math:`N_1` in the control net of the
      first B |eacute| zier triangle.
   :type num_nodes1: const int*
   :param nodes1:
      **[Input]** The actual control net of the first B |eacute| zier triangle
      as a :math:`2 \times N_1` array. This should be laid out in Fortran
      order, with :math:`2 N_1` total values.
   :type nodes1: const double*
   :param degree1:
      **[Input]** The degree :math:`d_1` of the first B |eacute| zier triangle.
   :type degree1: const int*
   :param num_nodes2:
      **[Input]** The number of nodes :math:`N_2` in the control net of the
      second B |eacute| zier triangle.
   :type num_nodes2: const int*
   :param nodes2:
      **[Input]** The actual control net of the second B |eacute| zier triangle
      as a :math:`2 \times N_2` array. This should be laid out in Fortran
      order, with :math:`2 N_2` total values.
   :type nodes2: const double*
   :param degree2:
      **[Input]** The degree :math:`d_2` of the second B |eacute| zier triangle.
   :type degree1: const int*
   :param segment_ends_size:
      **[Input]** The size of ``segment_ends``, which must be pre-allocated by
      the caller.
   :type segment_ends_size: const int*
   :param int* segment_ends:
      **[Output]** An array (pre-allocated by the caller) of the end indices
      for each group of segments in ``segments``. For example, if the triangles
      intersect in two distinct curved polygons, the first of which has four
      sides and the second of which has three, then the first two values in
      ``segment_ends`` will be ``[4, 7]`` and ``num_intersected`` will be
      ``2``.
   :param segments_size:
      **[Input]** The size of ``segments``, which must be pre-allocated by
      the caller.
   :type segments_size: const int*
   :param CurvedPolygonSegment* segments:
      **[Output]** An array (pre-allocated by the caller) of the edge segments
      that make up the boundary of the curved polygon(s) that form the
      intersection of the two triangles.
   :param int* num_intersected:
      **[Output]** The number of curved polygons in the intersection of two
      triangles.
   :param TriangleContained* contained:
      **[Output]** Enum indicating if one triangle is **fully** contained in
      the other.
   :param Status* status:
      **[Output]** The status code for the procedure. Will be

      * :c:data:`SUCCESS` on success.
      * :c:data:`INSUFFICIENT_SPACE` if ``segment_ends_size`` is smaller than
        ``num_intersected`` **or** if ``segments_size`` is smaller than the
        number of segments.
      * :c:data:`UNKNOWN` if the intersection points are classified in an
        unexpected way (e.g. if there is both an ignored corner and a tangent
        intersection, but no other types).
      * :c:data:`NO_CONVERGE` if the two curves in an edge pair don't converge
        to approximately linear after being subdivided 20 times. (This error
        will occur via :c:func:`BEZ_curve_intersections`.)
      * An integer :math:`N_C \geq 64` to indicate that there were :math:`N_C`
        pairs of candidate segments during edge-edge intersection that had
        overlapping convex hulls. This is a sign of either round-off error
        in detecting that the edges are coincident curve segments on the same
        algebraic curve or that the intersection is a non-simple root. (This
        error will occur via :c:func:`BEZ_curve_intersections`.)
      * :c:data:`BAD_MULTIPLICITY` if the two curves in an edge pair have an
        intersection that doesn't converge to either a simple or double root
        via Newton's method. (This error will occur via
        :c:func:`BEZ_curve_intersections`.)
      * :c:data:`EDGE_END` If there is an attempt to add an intersection
        point with either the :math:`s` or :math:`t`\-parameter equal to 1
        (i.e. if the intersection is at the end of an edge). This should
        not occur because such intersections are "rotated" to the beginning
        of the neighboring edge before the boundary of the curved polygon
        is formed.
      * :c:data:`SAME_CURVATURE` if the two curves in an edge pair have
        identical curvature at a tangent intersection.
      * :c:data:`BAD_INTERIOR` if a curved polygon requires more than
        10 sides. This could be due to either a particular complex
        intersection, a programming error or round-off which causes an
        infinite loop of intersection points to be added without wrapping
        around back to the first intersection point.

   **Signature:**

   .. code-block:: c

      void
      BEZ_triangle_intersections(const int *num_nodes1,
                                const double *nodes1,
                                const int *degree1,
                                const int *num_nodes2,
                                const double *nodes2,
                                const int *degree2,
                                const int *segment_ends_size,
                                int *segment_ends,
                                const int *segments_size,
                                CurvedPolygonSegment *segments,
                                int *num_intersected,
                                TriangleContained *contained,
                                Status *status);

.. c:function:: void BEZ_free_triangle_intersections_workspace(void)

   This frees any long-lived workspace(s) used by ``libbezier`` throughout
   the life of a program. It should be called during clean-up for any code
   which invokes :c:func:`BEZ_triangle_intersections`.

   **Signature:**

   .. code-block:: c

      void
      BEZ_free_triangle_intersections_workspace(void);

*****
Types
*****

.. c:type:: CurvedPolygonSegment

   Describes an edge of a :class:`.CurvedPolygon` formed when intersecting
   two curved B |eacute| zier triangles. The edges of the intersection need
   not be an entire edge of one of the triangles. For example, an edge
   :math:`E(s)` may be restricted to
   :math:`E\left(\left[\frac{1}{4}, \frac{7}{8}\right]\right)`.

   .. c:type:: double start

      The start parameter of the segment. In the restriction
      :math:`E\left(\left[\frac{1}{4}, \frac{7}{8}\right]\right)`, the
      ``start`` would be ``0.25``.

   .. c:type:: double end

      The end parameter of the segment. In the restriction
      :math:`E\left(\left[\frac{1}{4}, \frac{7}{8}\right]\right)`, the
      ``end`` would be ``0.875``.

   .. c:type:: int edge_index

      An index describing which edge the segment falls on. The edges
      of the first triangle in the intersection are given index values
      of ``1``, ``2`` and ``3`` while those of the second triangle are
      ``4``, ``5`` and ``6``.

   In the header ``bezier/triangle_intersection.h``, this is defined as

   .. code-block:: c

      typedef struct CurvedPolygonSegment {
        double start;
        double end;
        int edge_index;
      } CurvedPolygonSegment;

.. c:type:: TriangleContained

   This enum is used to indicate if one triangle is contained in
   another when doing triangle-triangle intersection.

   .. c:var:: NEITHER

      (``0``)
      Indicates that neither triangle is contained in the other. This
      could mean the triangles are disjoint or that they intersect
      in a way other than full containment.

   .. c:var:: FIRST

      (``1``)
      Indicates that the first triangle (arguments will be ordered) is
      fully contained in the second. This allows for points of tangency,
      shared corners or shared segments along an edge.

   .. c:var:: SECOND

      (``2``)
      Indicates that the second triangle (arguments will be ordered) is
      fully contained in the first. This allows for points of tangency,
      shared corners or shared segments along an edge.
