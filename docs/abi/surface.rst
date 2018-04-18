##############
surface module
##############

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

This is a collection of procedures for performing computations on a
B |eacute| zier `surface`_.

.. _surface: https://en.wikipedia.org/wiki/B%C3%A9zier_triangle

.. note::

    In most of the procedures below both the number of nodes :math:`N` and
    the degree :math:`d` of a B |eacute| zier surface are provided. It is
    redundant to require both as arguments since :math:`N = \binom{d + 2}{2}`.
    However, both are provided as arguments to avoid unnecessary
    re-computation, i.e. we expect the caller to know both :math:`N` and
    :math:`d`.

**********
Procedures
**********

.. c:function:: void compute_area(int *num_edges, \
                                  int *sizes, \
                                  double **nodes_pointers, \
                                  double *area, \
                                  bool *not_implemented)

   This computes the area of a curved polygon in :math:`\mathbf{R}^2` via
   Green's theorem. In order to do this, it assumes that the edges

   * form a closed loop
   * have no intersections with another edge (including self)
   * are oriented in the direction of the exterior

   :param int* num_edges:
      **[Input]** The number of edges :math:`N` that bound the curved polygon.
   :param int* sizes:
      **[Input]** An array of the size (i.e. the number of nodes) of each
      of the :math:`N` edges.
   :param double** nodes_pointers:
      **[Input]** An array of :math:`N` pointers. Pointer :math:`j` points to
      the start of the edge :math:`E_j` along the boundary of the curved
      polygon.
   :param double* area:
      **[Output]** The computed area of the curved polygon. If
      ``not_implemented`` is ``TRUE``, then this is undefined.
   :param bool* not_implemented:
      **[Output]** Indicates if a Green's theorem implementation exists for
      each of the edges. (Currently, the only degrees supported are 1, 2,
      3 and 4.)

   **Signature:**

   .. code-block:: c

      void
      compute_area(int *num_edges,
                   int *sizes,
                   double **nodes_pointers,
                   double *area,
                   bool *not_implemented);

.. c:function:: void compute_edge_nodes(int *num_nodes, \
                                        int *dimension, \
                                        double *nodes, \
                                        int *degree, \
                                        double *nodes1, \
                                        double *nodes2, \
                                        double *nodes3)

   Extracts the edge nodes from the control net of a B |eacute| zier surface.
   For example, if the control net of a quadratic surface is:

   .. math::

      \left[\begin{array}{c c c c c c}
          v_{2,0,0} & v_{1,1,0} &
          v_{0,2,0} & v_{1,0,1} &
          v_{0,1,1} & v_{0,0,2} \end{array}\right] =
      \left[\begin{array}{c c c c c c}
          a & b & c & d & e & f \end{array}\right]

   then the edges are

   .. math::

      \begin{align*}
      E_1 &= \left[\begin{array}{c c c}
          a & b & c \end{array}\right] \\
      E_2 &= \left[\begin{array}{c c c}
          c & e & f \end{array}\right] \\
      E_3 &= \left[\begin{array}{c c c}
          f & d & a \end{array}\right].
      \end{align*}

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the surface lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param double* nodes1:
      **[Output]** The control points of the first edge B |eacute| zier curve
      as a :math:`D \times (d + 1)` array, laid out in Fortran order.
   :param double* nodes2:
      **[Output]** The control points of the second edge B |eacute| zier curve
      as a :math:`D \times (d + 1)` array, laid out in Fortran order.
   :param double* nodes3:
      **[Output]** The control points of the third edge B |eacute| zier curve
      as a :math:`D \times (d + 1)` array, laid out in Fortran order.

   **Signature:**

   .. code-block:: c

      void
      compute_edge_nodes(int *num_nodes,
                         int *dimension,
                         double *nodes,
                         int *degree,
                         double *nodes1,
                         double *nodes2,
                         double *nodes3);

.. c:function:: void de_casteljau_one_round(int *num_nodes, \
                                            int *dimension, \
                                            double *nodes, \
                                            int *degree, \
                                            double *lambda1, \
                                            double *lambda2, \
                                            double *lambda3, \
                                            double *new_nodes)

   This performs a single round of the de Casteljau algorithm for evaluation
   in barycentric coordinates :math:`B(\lambda_1, \lambda_2, \lambda_3)`. This
   reduces the control net :math:`v_{i, j, k}^d` to a lower degree control net

   .. math::

      v_{i, j, k}^{d - 1} = \lambda_1 v_{i + 1, j, k}^d +
          \lambda_2 v_{i, j + 1, k}^d + \lambda_3 v_{i, j, k + 1}^d.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the surface lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param double* lambda1:
      **[Input]** The first barycentric parameter along the reference triangle.
   :param double* lambda2:
      **[Input]** The second barycentric parameter along the reference
      triangle.
   :param double* lambda3:
      **[Input]** The third barycentric parameter along the reference triangle.
   :param double* new_nodes:
      **[Output]** The newly-formed degree :math:`d - 1` control net. This will
      be a :math:`D \times (N - d - 1)` array.

   **Signature:**

   .. code-block:: c

      void
      de_casteljau_one_round(int *num_nodes,
                             int *dimension,
                             double *nodes,
                             int *degree,
                             double *lambda1,
                             double *lambda2,
                             double *lambda3,
                             double *new_nodes);

.. c:function:: void evaluate_barycentric(int *num_nodes, \
                                          int *dimension, \
                                          double *nodes, \
                                          int *degree, \
                                          double *lambda1, \
                                          double *lambda2, \
                                          double *lambda3, \
                                          double *point)

   Evaluates a single point on a B |eacute| zier surface, with input
   in barycentric coordinates: :math:`B(\lambda_1, \lambda_2, \lambda_3)`.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the surface lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param double* lambda1:
      **[Input]** The first barycentric parameter along the reference triangle.
   :param double* lambda2:
      **[Input]** The second barycentric parameter along the reference
      triangle.
   :param double* lambda3:
      **[Input]** The third barycentric parameter along the reference triangle.
   :param double* point:
      **[Output]** A :math:`D \times 1` array, will contain
      :math:`B(\lambda_1, \lambda_2, \lambda_3)`.

   **Signature:**

   .. code-block:: c

      void
      evaluate_barycentric(int *num_nodes,
                           int *dimension,
                           double *nodes,
                           int *degree,
                           double *lambda1,
                           double *lambda2,
                           double *lambda3,
                           double *point);

.. c:function:: void evaluate_barycentric_multi(int *num_nodes, \
                                                int *dimension, \
                                                double *nodes, \
                                                int *degree, \
                                                int *num_vals, \
                                                double *param_vals, \
                                                double *evaluated)

   Evaluates many points on a B |eacute| zier surface, with input
   in barycentric coordinates:
   :math:`\left\{B(\lambda_{1,j}, \lambda_{2,j}, \lambda_{3,j})\right\}_j`.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the surface lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param int* num_vals:
      **[Input]** The number of points :math:`k` where :math:`B` is
      being evaluated.
   :param double* param_vals:
      **[Input]** A :math:`k \times 3` array of :math:`k` triples of
      barycentric coordinates, laid out in Fortran order. This way, the
      first column contains all :math:`\lambda_1` values in contiguous order,
      and similarly for the other columns.
   :param double* evaluated:
      **[Output]** A :math:`D \times k` array of all evaluated points on the
      surface. Column :math:`j` will contain
      :math:`B(\lambda_{1,j}, \lambda_{2,j}, \lambda_{3,j})`.

   **Signature:**

   .. code-block:: c

      void
      evaluate_barycentric_multi(int *num_nodes,
                                 int *dimension,
                                 double *nodes,
                                 int *degree,
                                 int *num_vals,
                                 double *param_vals,
                                 double *evaluated);

.. c:function:: void evaluate_cartesian_multi(int *num_nodes, \
                                              int *dimension, \
                                              double *nodes, \
                                              int *degree, \
                                              int *num_vals, \
                                              double *param_vals, \
                                              double *evaluated)

   Evaluates many points on a B |eacute| zier surface, with input
   in cartesian coordinates:
   :math:`\left\{B(s_j, t_j)\right\}_j`. Each input :math:`(s, t)` is
   equivalent to the barycentric input :math:`\lambda_1 = 1 - s - t`,
   :math:`\lambda_2 = s` and :math:`\lambda_3 = t`.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the surface lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param int* num_vals:
      **[Input]** The number of points :math:`k` where :math:`B` is
      being evaluated.
   :param double* param_vals:
      **[Input]** A :math:`k \times 2` array of :math:`k` pairs of
      cartesian coordinates, laid out in Fortran order. This way, the
      first column contains all :math:`s`\-values in contiguous order,
      and similarly for the other column.
   :param double* evaluated:
      **[Output]** A :math:`D \times k` array of all evaluated points on the
      surface. Column :math:`j` will contain
      :math:`B(s_j, t_j)`.

   **Signature:**

   .. code-block:: c

      void
      evaluate_cartesian_multi(int *num_nodes,
                               int *dimension,
                               double *nodes,
                               int *degree,
                               int *num_vals,
                               double *param_vals,
                               double *evaluated);

.. c:function:: void jacobian_both(int *num_nodes, \
                                   int *dimension, \
                                   double *nodes, \
                                   int *degree, \
                                   double *new_nodes)

   Computes control nets for both cartesian partial derivatives of a
   B |eacute| ezier surface :math:`B_s(s, t)` and :math:`B_t(s, t)`. Taking
   a single (partial) derivative lowers the degree by 1.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the surface lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param double* new_nodes:
      **[Output]** The combined control nets :math:`B_s` and :math:`B_t` as
      a :math:`(2D) \times (N - d - 1)` array, laid out in Fortran order. The
      first :math:`D` columns contain the control net of :math:`B_s` and
      final :math:`D` columns contain the control net of :math:`B_t`.

   **Signature:**

   .. code-block:: c

      void
      jacobian_both(int *num_nodes,
                    int *dimension,
                    double *nodes,
                    int *degree,
                    double *new_nodes);

.. c:function:: void jacobian_det(int *num_nodes, \
                                  double *nodes, \
                                  int *degree, \
                                  int *num_vals, \
                                  double *param_vals, \
                                  double *evaluated)

   Computes :math:`\det(DB)` at a many points :math:`(s_j, t_j)`. This is only
   well-defined if :math:`\det(DB)` has two rows, hence the surface must lie
   in :math:`\mathbf{R}^2`.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`2 \times N` array. This should be laid out in Fortran order, with
      :math:`2 N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param int* num_vals:
      **[Input]** The number of points :math:`k` where :math:`\det(DB)` is
      being evaluated.
   :param double* param_vals:
      **[Input]** A :math:`k \times 2` array of :math:`k` pairs of
      cartesian coordinates, laid out in Fortran order. This way, the
      first column contains all :math:`s`\-values in contiguous order,
      and similarly for the other column.
   :param double* evaluated:
      **[Output]** A :math:`k` array of all evaluated determinants. The
      surface. The :math:`j`\-th value will be :math:`\det(DB(s_j, t_j))`.

   **Signature:**

   .. code-block:: c

      void
      jacobian_det(int *num_nodes,
                   double *nodes,
                   int *degree,
                   int *num_vals,
                   double *param_vals,
                   double *evaluated);

.. c:function:: void specialize_surface(int *num_nodes, \
                                        int *dimension, \
                                        double *nodes, \
                                        int *degree, \
                                        double *weights_a, \
                                        double *weights_b, \
                                        double *weights_c, \
                                        double *specialized)

   Changes the control net for a B |eacute| zier surface by specializing
   from the original triangle :math:`(0, 0), (1, 0), (0, 1)` to a new
   triangle :math:`p_1, p_2, p_3`.

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the surface lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param double* weights_a:
      **[Input]** A 3-array containing the barycentric weights for the first
      node :math:`p_1` in the new triangle.
   :param double* weights_b:
      **[Input]** A 3-array containing the barycentric weights for the second
      node :math:`p_2` in the new triangle.
   :param double* weights_c:
      **[Input]** A 3-array containing the barycentric weights for the third
      node :math:`p_3` in the new triangle.
   :param double* specialized:
      **[Input]** The control net of the newly formed B |eacute| zier surface
      as a :math:`D \times N` array.

   **Signature:**

   .. code-block:: c

      void
      specialize_surface(int *num_nodes,
                         int *dimension,
                         double *nodes,
                         int *degree,
                         double *weights_a,
                         double *weights_b,
                         double *weights_c,
                         double *specialized);

.. c:function:: void subdivide_nodes_surface(int *num_nodes, \
                                             int *dimension, \
                                             double *nodes, \
                                             int *degree, \
                                             double *nodes_a, \
                                             double *nodes_b, \
                                             double *nodes_c, \
                                             double *nodes_d)

   Subdivides a B |eacute| surface into four sub-surfaces that cover
   the original surface. See :meth:`.Surface.subdivide` for more
   details

   :param int* num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier surface.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the surface lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control net of the B |eacute| zier surface as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :param int* degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier surface.
   :param double* nodes_a:
      **[Output]** The control net of the lower left sub-surface as a
      :math:`D \times N` array, laid out in Fortran order.
   :param double* nodes_b:
      **[Output]** The control net of the central sub-surface as a
      :math:`D \times N` array, laid out in Fortran order.
   :param double* nodes_c:
      **[Output]** The control net of the lower right sub-surface as a
      :math:`D \times N` array, laid out in Fortran order.
   :param double* nodes_d:
      **[Output]** The control net of the upper left sub-surface as a
      :math:`D \times N` array, laid out in Fortran order.

   **Signature:**

   .. code-block:: c

      void
      subdivide_nodes_surface(int *num_nodes,
                              int *dimension,
                              double *nodes,
                              int *degree,
                              double *nodes_a,
                              double *nodes_b,
                              double *nodes_c,
                              double *nodes_d);
