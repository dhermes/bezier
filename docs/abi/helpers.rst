##############
helpers module
##############

This is a collection of generic utility procedures that help
with various geometric computations.

**********
Procedures
**********

.. c:function:: void BEZ_bbox(const int *num_nodes, \
                              const double *nodes, \
                              double *left, \
                              double *right, \
                              double *bottom, \
                              double *top)

   Computes a rectangular bounding box
   :math:`\left[m_x, M_x\right] \times \left[m_y, M_y\right]`
   for a set of points in :math:`\mathbf{R}^2`.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` we are bounding.
   :type num_nodes: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual nodes as a :math:`2 \times N` array. This should
      be laid out in Fortran order, with :math:`2 N` total values.
   :type nodes: :c:expr:`const double*`
   :param left:
      **[Output]** The minimal :math:`x`\-value :math:`m_x`.
   :type left: :c:expr:`double*`
   :param right:
      **[Output]** The maximal :math:`x`\-value :math:`M_x`.
   :type right: :c:expr:`double*`
   :param bottom:
      **[Output]** The minimal :math:`y`\-value :math:`m_y`.
   :type bottom: :c:expr:`double*`
   :param top:
      **[Output]** The maximal :math:`y`\-value :math:`M_y`.
   :type top: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_bbox(const int *num_nodes,
               const double *nodes,
               double *left,
               double *right,
               double *bottom,
               double *top);

.. c:function:: void BEZ_contains_nd(const int *num_nodes, \
                                     const int *dimension, \
                                     const double *nodes, \
                                     const double *point, \
                                     bool *predicate)

   Checks if a point :math:`p` is contained in the bounding box
   for a set of points.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` define the bounding box.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the nodes and point lie in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual nodes as a :math:`D \times N` array. This should
      be laid out in Fortran order, with :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param point:
      **[Input]** The point :math:`p` as an array containing :math:`D` values.
   :type point: :c:expr:`const double*`
   :param predicate:
      **[Output]** Flag indicating if the point lies inside the box.
   :type predicate: :c:expr:`bool*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_contains_nd(const int *num_nodes,
                      const int *dimension,
                      const double *nodes,
                      const double *point,
                      bool *predicate);

.. c:function:: void BEZ_cross_product(const double *vec0, \
                                       const double *vec1, \
                                       double *result)

   Computes the cross-product of two vectors :math:`v_1, v_2` in
   :math:`\mathbf{R}^2`. This is done as if they were embedded in
   :math:`\mathbf{R}^3` and the result is the resulting :math:`z`\-component
   :math:`x_1 y_2 - x_2 y_1`.

   :param vec0:
      **[Input]** The first vector :math:`v_1` in :math:`\mathbf{R}^2`.
   :type vec0: :c:expr:`const double*`
   :param vec1:
      **[Input]** The second vector :math:`v_2` in :math:`\mathbf{R}^2`.
   :type vec1: :c:expr:`const double*`
   :param result:
      **[Output]** The cross-product.
   :type result: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_cross_product(const double *vec0,
                        const double *vec1,
                        double *result);

.. c:function:: bool BEZ_in_interval(const double *value, \
                                     const double *start, \
                                     const double *end)

   Checks if a value :math:`v` is in an interval :math:`\left[s, e\right]`.

   :param value:
      **[Input]** The value :math:`v`.
   :type value: :c:expr:`const double*`
   :param start:
      **[Input]** The start :math:`s` of the interval
      :math:`\left[s, e\right]`.
   :type start: :c:expr:`const double*`
   :param end:
      **[Input]** The end :math:`e` of the interval :math:`\left[s, e\right]`.
   :type end: :c:expr:`const double*`
   :returns: Flag indicating if :math:`v \in \left[s, e\right]`.
   :rtype: bool

   **Signature:**

   .. code-block:: c

      bool
      BEZ_in_interval(const double *value,
                      const double *start,
                      const double *end);

.. c:function:: void BEZ_polygon_collide(const int *polygon_size1, \
                                         const double *polygon1, \
                                         const int *polygon_size2, \
                                         const double *polygon2, \
                                         bool *collision)

   Determines if two polygons collide.

   :param polygon_size1:
      **[Input]** The number of sides :math:`N_1` in the first polygon.
   :type polygon_size1: :c:expr:`const int*`
   :param polygon1:
      **[Input]** The nodes of the first polygon as a :math:`2 \times N_1`
      array. This should be laid out in Fortran order.
   :type polygon1: :c:expr:`const double*`
   :param polygon_size2:
      **[Input]** The number of sides :math:`N_2` in the second polygon.
   :type polygon_size2: :c:expr:`const int*`
   :param polygon2:
      **[Input]** The nodes of the second polygon as a :math:`2 \times N_2`
      array. This should be laid out in Fortran order.
   :type polygon2: :c:expr:`const double*`
   :param collision:
      **[Output]** Flag indicating if the polygons collide.
   :type collision: :c:expr:`bool*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_polygon_collide(const int *polygon_size1,
                          const double *polygon1,
                          const int *polygon_size2,
                          const double *polygon2,
                          bool *collision);

.. c:function:: void BEZ_simple_convex_hull(const int *num_points, \
                                            const double *points, \
                                            int *polygon_size, \
                                            double *polygon)

   Computes the convex hull of a set of points.

   :param num_points:
      **[Input]** The number of points :math:`N`.
   :type num_points: :c:expr:`const int*`
   :param points:
      **[Input]** The points being considered, as a :math:`2 \times N`
      array. This should be laid out in Fortran order.
   :type points: :c:expr:`const double*`
   :param polygon_size:
      **[Output]** The number of sides :math:`M` in the convex hull. This
      will be at most :math:`N`.
   :type polygon_size: :c:expr:`int*`
   :param polygon:
      **[Output]** The nodes in the convex hull, as a :math:`2 \times N`
      array laid out in Fortran order. This must be allocated by the caller
      and must be size :math:`N` to account for the extreme case.
   :type polygon: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_simple_convex_hull(const int *num_points,
                             const double *points,
                             int *polygon_size,
                             double *polygon);

.. c:function:: bool BEZ_vector_close(const int *num_values, \
                                      const double *vec1, \
                                      const double *vec2, \
                                      const double *eps)

   Determines if two vectors are close to machine precision.

   :param num_values:
      **[Input]** The dimension :math:`D` such that the vectors lie in
      :math:`\mathbf{R}^D`.
   :type num_values: :c:expr:`const int*`
   :param vec1:
      **[Input]** The first vector :math:`v_1`, as an array of :math:`D`
      values.
   :type vec1: :c:expr:`const double*`
   :param vec2:
      **[Input]** The second vector :math:`v_2`, as an array of :math:`D`
      values.
   :type vec2: :c:expr:`const double*`
   :param eps:
      **[Input]** The tolerance :math:`\varepsilon` used when comparing
      :math:`\left\lVert v_1 - v_2 \right\rVert` to
      :math:`\left\lVert v_1 \right\rVert` and
      :math:`\left\lVert v_2 \right\rVert`.
   :type eps: :c:expr:`const double*`
   :returns:
      Flag indicating if :math:`v_1` and :math:`v_2` are close to the desired
      precision.
   :rtype: bool

   **Signature:**

   .. code-block:: c

      bool
      BEZ_vector_close(const int *num_values,
                       const double *vec1,
                       const double *vec2,
                       const double *eps);

.. c:function:: void BEZ_wiggle_interval(const double *value, \
                                         double *result, \
                                         bool *success)

   Round a value :math:`v` into the unit interval if it is sufficiently
   close. The real line will be broken into five intervals and handled
   differently on each interval:

   * :math:`v \in \left(-\infty, -2^{-44}\right]` will not be rounded
     and will set ``success`` to ``FALSE``.
   * :math:`v \in \left(-2^{-44}, 2^{-44}\right)` will be rounded to
     ``0.0``.
   * :math:`v \in \left[2^{-44}, 1 - 2^{-44}\right]` will be left
     untouched (i.e. they are safely in the unit interval).
   * :math:`v \in \left(1 - 2^{-44}, 1 + 2^{-44}\right)` will be rounded to
     ``1.0``.
   * :math:`v \in \left[1 + 2^{-44}, \infty\right)` will not be rounded
     and will set ``success`` to ``FALSE``.

   :param value:
      **[Input]** The value :math:`v` to be rounded.
   :type value: :c:expr:`const double*`
   :param result:
      **[Output]** The rounded version of :math:`v`. If ``success`` is
      ``FALSE`` this is undefined.
   :type result: :c:expr:`double*`
   :param success:
      **[Output]** Flag indicating if :math:`v` was in the unit interval or
      sufficiently close to it.
   :type success: :c:expr:`bool*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_wiggle_interval(const double *value,
                          double *result,
                          bool *success);
