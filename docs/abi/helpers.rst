##############
helpers module
##############

This is a collection of generic utility procedures that help
with various geometric computations.

**********
Procedures
**********

.. c:function:: void bbox(int *num_nodes, \
                          double *nodes, \
                          double *left, \
                          double *right, \
                          double *bottom, \
                          double *top)

   Computes a rectangular bounding box
   :math:`\left[m_x, M_x\right] \times \left[m_y, M_y\right]`
   for a set of points in :math:`\mathbf{R}^2`.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` we are bounding.
   :param double* nodes:
      **[Input]** The actual nodes as a :math:`2 \times N` array. This should
      be laid out in Fortran order, with :math:`2 N` total values.
   :param double* left:
      **[Output]** The minimal :math:`x`\-value :math:`m_x`.
   :param double* right:
      **[Output]** The maximal :math:`x`\-value :math:`M_x`.
   :param double* bottom:
      **[Output]** The minimal :math:`y`\-value :math:`m_y`.
   :param double* top:
      **[Output]** The maximal :math:`y`\-value :math:`M_y`.

   **Signature:**

   .. code-block:: c

      void
      bbox(int *num_nodes,
           double *nodes,
           double *left,
           double *right,
           double *bottom,
           double *top);

.. c:function:: void contains_nd(int *num_nodes, \
                                 int *dimension, \
                                 double *nodes, \
                                 double *point, \
                                 bool *predicate)

   Checks if a point :math:`p` is contained in the bounding box
   for a set of points.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` define the bounding box.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the nodes and point lie in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual nodes as a :math:`D \times N` array. This should
      be laid out in Fortran order, with :math:`D N` total values.
   :param double* point:
      **[Input]** The point :math:`p` as an array containing :math:`D` values.
   :param bool* predicate:
      **[Output]** Flag indicating if the point lies inside the box.

   **Signature:**

   .. code-block:: c

      void
      contains_nd(int *num_nodes,
                  int *dimension,
                  double *nodes,
                  double *point,
                  bool *predicate);

.. c:function:: void cross_product(double *vec0, \
                                   double *vec1, \
                                   double *result)

   Computes the cross-product of two vectors :math:`v_1, v_2` in
   :math:`\mathbf{R}^2`. This is done as if they were embedded in
   :math:`\mathbf{R}^3` and the result is the resulting :math:`z`\-component
   :math:`x_1 y_2 - x_2 y_1`.

   :param double* vec0:
      **[Input]** The first vector :math:`v_1` in :math:`\mathbf{R}^2`.
   :param double* vec1:
      **[Input]** The second vector :math:`v_2` in :math:`\mathbf{R}^2`.
   :param double* result:
      **[Output]** The cross-product.

   **Signature:**

   .. code-block:: c

      void
      cross_product(double *vec0,
                    double *vec1,
                    double *result);

.. c:function:: bool in_interval(double *value, \
                                 double *start, \
                                 double *end)

   Checks if a value :math:`v` is in an interval :math:`\left[s, e\right]`.

   :param double* value:
      **[Input]** The value :math:`v`.
   :param double* start:
      **[Input]** The start :math:`s` of the interval
      :math:`\left[s, e\right]`.
   :param double* end:
      **[Input]** The end :math:`e` of the interval :math:`\left[s, e\right]`.
   :returns: Flag indicating if :math:`v \in \left[s, e\right]`.
   :rtype: bool

   **Signature:**

   .. code-block:: c

      bool
      in_interval(double *value,
                  double *start,
                  double *end);

.. c:function:: void polygon_collide(int *polygon_size1, \
                                     double *polygon1, \
                                     int *polygon_size2, \
                                     double *polygon2, \
                                     bool *collision)

   Determines if two polygons collide.

   :param int* polygon_size1:
      **[Input]** The number of sides :math:`N_1` in the first polygon.
   :param double* polygon1:
      **[Input]** The nodes of the first polygon as a :math:`2 \times N_1`
      array. This should be laid out in Fortran order.
   :param int* polygon_size2:
      **[Input]** The number of sides :math:`N_2` in the second polygon.
   :param double* polygon2:
      **[Input]** The nodes of the second polygon as a :math:`2 \times N_2`
      array. This should be laid out in Fortran order.
   :param bool* collision:
      **[Output]** Flag indicating if the polygons collide.

   **Signature:**

   .. code-block:: c

      void
      polygon_collide(int *polygon_size1,
                      double *polygon1,
                      int *polygon_size2,
                      double *polygon2,
                      bool *collision);

.. c:function:: void simple_convex_hull(int *num_points, \
                                        double *points, \
                                        int *polygon_size, \
                                        double *polygon)

   Computes the convex hull of a set of points.

   :param int* num_points:
      **[Input]** The number of points :math:`N`.
   :param double* points:
      **[Input]** The points being considered, as a :math:`2 \times N`
      array. This should be laid out in Fortran order.
   :param int* polygon_size:
      **[Output]** The number of sides :math:`M` in the convex hull. This
      will be at most :math:`N`.
   :param double* polygon:
      **[Output]** The nodes in the convex hull, as a :math:`2 \times N`
      array laid out in Fortran order. This must be allocated by the caller
      and must be size :math:`N` to account for the extreme case.

   **Signature:**

   .. code-block:: c

      void
      simple_convex_hull(int *num_points,
                         double *points,
                         int *polygon_size,
                         double *polygon);

.. c:function:: bool vector_close(int *num_values, \
                                  double *vec1, \
                                  double *vec2, \
                                  double *eps)

   Determines if two vectors are close to machine precision.

   :param int* num_values:
      **[Input]** The dimension :math:`D` such that the vectors lie in
      :math:`\mathbf{R}^D`.
   :param double* vec1:
      **[Input]** The first vector :math:`v_1`, as an array of :math:`D`
      values.
   :param double* vec2:
      **[Input]** The second vector :math:`v_2`, as an array of :math:`D`
      values.
   :param double* eps:
      **[Input]** The tolerance :math:`\varepsilon` used when comparing
      :math:`\left\lVert v_1 - v_2 \right\rVert` to
      :math:`\left\lVert v_1 \right\rVert` and
      :math:`\left\lVert v_2 \right\rVert`.
   :returns:
      Flag indicating if :math:`v_1` and :math:`v_2` are close to the desired
      precision.
   :rtype: bool

   **Signature:**

   .. code-block:: c

      bool
      vector_close(int *num_values,
                   double *vec1,
                   double *vec2,
                   double *eps);

.. c:function:: void wiggle_interval(double *value, \
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

   :param double* value:
      **[Input]** The value :math:`v` to be rounded.
   :param double* result:
      **[Output]** The rounded version of :math:`v`. If ``success`` is
      ``FALSE`` this is undefined.
   :param bool* success:
      **[Output]** Flag indicating if :math:`v` was in the unit interval or
      sufficiently close to it.

   **Signature:**

   .. code-block:: c

      void
      wiggle_interval(double *value,
                      double *result,
                      bool *success);
