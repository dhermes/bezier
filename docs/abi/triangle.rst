###############
triangle module
###############

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

This is a collection of procedures for performing computations on a
B |eacute| zier `triangle`_.

.. _triangle: https://en.wikipedia.org/wiki/B%C3%A9zier_triangle

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

.. c:function:: void BEZ_compute_area(const int *num_edges, \
                                      const int *sizes, \
                                      const double* const* nodes_pointers, \
                                      double *area, \
                                      bool *not_implemented)

   This computes the area of a curved polygon in :math:`\mathbf{R}^2` via
   Green's theorem. In order to do this, it assumes that the edges

   * form a closed loop
   * have no intersections with another edge (including self)
   * are oriented in the direction of the exterior

   :param num_edges:
      **[Input]** The number of edges :math:`N` that bound the curved polygon.
   :type num_edges: :c:expr:`const int*`
   :param sizes:
      **[Input]** An array of the size (i.e. the number of nodes) of each
      of the :math:`N` edges.
   :type sizes: :c:expr:`const int*`
   :param nodes_pointers:
      **[Input]** An array of :math:`N` pointers. Pointer :math:`j` points to
      the start of the edge :math:`E_j` along the boundary of the curved
      polygon. (This is a ``const`` array of ``const`` pointers, i.e. the
      ``double**`` pointer array cannot be modified and the underlying
      ``double*`` arrays pointed to by each element cannot be modified.)
   :type nodes_pointers: :c:expr:`const double* const*`
   :param area:
      **[Output]** The computed area of the curved polygon. If
      ``not_implemented`` is ``TRUE``, then this is undefined.
   :type area: :c:expr:`double*`
   :param not_implemented:
      **[Output]** Indicates if a Green's theorem implementation exists for
      each of the edges. (Currently, the only degrees supported are 1, 2,
      3 and 4.)
   :type not_implemented: :c:expr:`bool*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_compute_area(const int *num_edges,
                       const int *sizes,
                       const double* const* nodes_pointers,
                       double *area,
                       bool *not_implemented);

.. c:function:: void BEZ_compute_edge_nodes(const int *num_nodes, \
                                            const int *dimension, \
                                            const double *nodes, \
                                            const int *degree, \
                                            double *nodes1, \
                                            double *nodes2, \
                                            double *nodes3)

   Extracts the edge nodes from the control net of a B |eacute| zier triangle.
   For example, if the control net of a quadratic triangle is:

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

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the triangle lies in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param nodes1:
      **[Output]** The control points of the first edge B |eacute| zier curve
      as a :math:`D \times (d + 1)` array, laid out in Fortran order.
   :type nodes1: :c:expr:`double*`
   :param nodes2:
      **[Output]** The control points of the second edge B |eacute| zier curve
      as a :math:`D \times (d + 1)` array, laid out in Fortran order.
   :type nodes2: :c:expr:`double*`
   :param nodes3:
      **[Output]** The control points of the third edge B |eacute| zier curve
      as a :math:`D \times (d + 1)` array, laid out in Fortran order.
   :type nodes3: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_compute_edge_nodes(const int *num_nodes,
                             const int *dimension,
                             const double *nodes,
                             const int *degree,
                             double *nodes1,
                             double *nodes2,
                             double *nodes3);

.. c:function:: void BEZ_de_casteljau_one_round(const int *num_nodes, \
                                                const int *dimension, \
                                                const double *nodes, \
                                                const int *degree, \
                                                const double *lambda1, \
                                                const double *lambda2, \
                                                const double *lambda3, \
                                                double *new_nodes)

   This performs a single round of the de Casteljau algorithm for evaluation
   in barycentric coordinates :math:`B(\lambda_1, \lambda_2, \lambda_3)`. This
   reduces the control net :math:`v_{i, j, k}^d` to a lower degree control net

   .. math::

      v_{i, j, k}^{d - 1} = \lambda_1 v_{i + 1, j, k}^d +
          \lambda_2 v_{i, j + 1, k}^d + \lambda_3 v_{i, j, k + 1}^d.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the triangle lies in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param lambda1:
      **[Input]** The first barycentric parameter along the reference triangle.
   :type lambda1: :c:expr:`const double*`
   :param lambda2:
      **[Input]** The second barycentric parameter along the reference
      triangle.
   :type lambda2: :c:expr:`const double*`
   :param lambda3:
      **[Input]** The third barycentric parameter along the reference triangle.
   :type lambda3: :c:expr:`const double*`
   :param new_nodes:
      **[Output]** The newly-formed degree :math:`d - 1` control net. This will
      be a :math:`D \times (N - d - 1)` array.
   :type new_nodes: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_de_casteljau_one_round(const int *num_nodes,
                                 const int *dimension,
                                 const double *nodes,
                                 const int *degree,
                                 const double *lambda1,
                                 const double *lambda2,
                                 const double *lambda3,
                                 double *new_nodes);

.. c:function:: void BEZ_evaluate_barycentric(const int *num_nodes, \
                                              const int *dimension, \
                                              const double *nodes, \
                                              const int *degree, \
                                              const double *lambda1, \
                                              const double *lambda2, \
                                              const double *lambda3, \
                                              double *point)

   Evaluates a single point on a B |eacute| zier triangle, with input
   in barycentric coordinates: :math:`B(\lambda_1, \lambda_2, \lambda_3)`.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the triangle lies in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param lambda1:
      **[Input]** The first barycentric parameter along the reference triangle.
   :type lambda1: :c:expr:`const double*`
   :param lambda2:
      **[Input]** The second barycentric parameter along the reference
      triangle.
   :type lambda2: :c:expr:`const double*`
   :param lambda3:
      **[Input]** The third barycentric parameter along the reference triangle.
   :type lambda3: :c:expr:`const double*`
   :param point:
      **[Output]** A :math:`D \times 1` array, will contain
      :math:`B(\lambda_1, \lambda_2, \lambda_3)`.
   :type point: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_evaluate_barycentric(const int *num_nodes,
                               const int *dimension,
                               const double *nodes,
                               const int *degree,
                               const double *lambda1,
                               const double *lambda2,
                               const double *lambda3,
                               double *point);

.. c:function:: void BEZ_evaluate_barycentric_multi(const int *num_nodes, \
                                                    const int *dimension, \
                                                    const double *nodes, \
                                                    const int *degree, \
                                                    const int *num_vals, \
                                                    const double *param_vals, \
                                                    double *evaluated)

   Evaluates many points on a B |eacute| zier triangle, with input
   in barycentric coordinates:
   :math:`\left\{B(\lambda_{1,j}, \lambda_{2,j}, \lambda_{3,j})\right\}_j`.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the triangle lies in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param num_vals:
      **[Input]** The number of points :math:`k` where :math:`B` is
      being evaluated.
   :type num_vals: :c:expr:`const int*`
   :param param_vals:
      **[Input]** A :math:`k \times 3` array of :math:`k` triples of
      barycentric coordinates, laid out in Fortran order. This way, the
      first column contains all :math:`\lambda_1` values in contiguous order,
      and similarly for the other columns.
   :type param_vals: :c:expr:`const double*`
   :param evaluated:
      **[Output]** A :math:`D \times k` array of all evaluated points on the
      triangle. Column :math:`j` will contain
      :math:`B(\lambda_{1,j}, \lambda_{2,j}, \lambda_{3,j})`.
   :type evaluated: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_evaluate_barycentric_multi(const int *num_nodes,
                                     const int *dimension,
                                     const double *nodes,
                                     const int *degree,
                                     const int *num_vals,
                                     const double *param_vals,
                                     double *evaluated);

.. c:function:: void BEZ_evaluate_cartesian_multi(const int *num_nodes, \
                                                  const int *dimension, \
                                                  const double *nodes, \
                                                  const int *degree, \
                                                  const int *num_vals, \
                                                  const double *param_vals, \
                                                  double *evaluated)

   Evaluates many points on a B |eacute| zier triangle, with input
   in cartesian coordinates:
   :math:`\left\{B(s_j, t_j)\right\}_j`. Each input :math:`(s, t)` is
   equivalent to the barycentric input :math:`\lambda_1 = 1 - s - t`,
   :math:`\lambda_2 = s` and :math:`\lambda_3 = t`.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the triangle lies in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param num_vals:
      **[Input]** The number of points :math:`k` where :math:`B` is
      being evaluated.
   :type num_vals: :c:expr:`const int*`
   :param param_vals:
      **[Input]** A :math:`k \times 2` array of :math:`k` pairs of
      cartesian coordinates, laid out in Fortran order. This way, the
      first column contains all :math:`s`\-values in contiguous order,
      and similarly for the other column.
   :type param_vals: :c:expr:`const double*`
   :param evaluated:
      **[Output]** A :math:`D \times k` array of all evaluated points on the
      triangle. Column :math:`j` will contain
      :math:`B(s_j, t_j)`.
   :type evaluated: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_evaluate_cartesian_multi(const int *num_nodes,
                                   const int *dimension,
                                   const double *nodes,
                                   const int *degree,
                                   const int *num_vals,
                                   const double *param_vals,
                                   double *evaluated);

.. c:function:: void BEZ_jacobian_both(const int *num_nodes, \
                                       const int *dimension, \
                                       const double *nodes, \
                                       const int *degree, \
                                       double *new_nodes)

   Computes control nets for both cartesian partial derivatives of a
   B |eacute| zier triangle :math:`B_s(s, t)` and :math:`B_t(s, t)`. Taking
   a single (partial) derivative lowers the degree by 1.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the triangle lies in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param new_nodes:
      **[Output]** The combined control nets :math:`B_s` and :math:`B_t` as
      a :math:`(2D) \times (N - d - 1)` array, laid out in Fortran order. The
      first :math:`D` columns contain the control net of :math:`B_s` and
      final :math:`D` columns contain the control net of :math:`B_t`.
   :type new_nodes: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_jacobian_both(const int *num_nodes,
                        const int *dimension,
                        const double *nodes,
                        const int *degree,
                        double *new_nodes);

.. c:function:: void BEZ_jacobian_det(const int *num_nodes, \
                                      const double *nodes, \
                                      const int *degree, \
                                      const int *num_vals, \
                                      const double *param_vals, \
                                      double *evaluated)

   Computes :math:`\det(DB)` at many points :math:`(s_j, t_j)`. This is only
   well-defined if :math:`\det(DB)` has two rows, hence the triangle must lie
   in :math:`\mathbf{R}^2`.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`2 \times N` array. This should be laid out in Fortran order, with
      :math:`2 N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param num_vals:
      **[Input]** The number of points :math:`k` where :math:`\det(DB)` is
      being evaluated.
   :type num_vals: :c:expr:`const int*`
   :param param_vals:
      **[Input]** A :math:`k \times 2` array of :math:`k` pairs of
      cartesian coordinates, laid out in Fortran order. This way, the
      first column contains all :math:`s`\-values in contiguous order,
      and similarly for the other column.
   :type param_vals: :c:expr:`const double*`
   :param evaluated:
      **[Output]** A :math:`k` array of all evaluated determinants. The
      :math:`j`\-th value will be :math:`\det(DB(s_j, t_j))`.
   :type evaluated: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_jacobian_det(const int *num_nodes,
                       const double *nodes,
                       const int *degree,
                       const int *num_vals,
                       const double *param_vals,
                       double *evaluated);

.. c:function:: void BEZ_specialize_triangle(const int *num_nodes, \
                                             const int *dimension, \
                                             const double *nodes, \
                                             const int *degree, \
                                             const double *weights_a, \
                                             const double *weights_b, \
                                             const double *weights_c, \
                                             double *specialized)

   Changes the control net for a B |eacute| zier triangle by specializing
   from the original triangle :math:`(0, 0), (1, 0), (0, 1)` to a new
   triangle :math:`p_1, p_2, p_3`.

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the triangle lies in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param weights_a:
      **[Input]** A 3-array containing the barycentric weights for the first
      node :math:`p_1` in the new triangle.
   :type weights_a: :c:expr:`const double*`
   :param weights_b:
      **[Input]** A 3-array containing the barycentric weights for the second
      node :math:`p_2` in the new triangle.
   :type weights_b: :c:expr:`const double*`
   :param weights_c:
      **[Input]** A 3-array containing the barycentric weights for the third
      node :math:`p_3` in the new triangle.
   :type weights_c: :c:expr:`const double*`
   :param specialized:
      **[Output]** The control net of the newly formed B |eacute| zier triangle
      as a :math:`D \times N` array.
   :type specialized: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_specialize_triangle(const int *num_nodes,
                              const int *dimension,
                              const double *nodes,
                              const int *degree,
                              const double *weights_a,
                              const double *weights_b,
                              const double *weights_c,
                              double *specialized);

.. c:function:: void BEZ_subdivide_nodes_triangle(const int *num_nodes, \
                                                  const int *dimension, \
                                                  const double *nodes, \
                                                  const int *degree, \
                                                  double *nodes_a, \
                                                  double *nodes_b, \
                                                  double *nodes_c, \
                                                  double *nodes_d)

   Subdivides a B |eacute| zier triangle into four sub-triangles that cover
   the original triangle. See :meth:`.Triangle.subdivide` for more
   details

   :param num_nodes:
      **[Input]** The number of nodes :math:`N` in the control net of the
      B |eacute| zier triangle.
   :type num_nodes: :c:expr:`const int*`
   :param dimension:
      **[Input]** The dimension :math:`D` such that the triangle lies in
      :math:`\mathbf{R}^D`.
   :type dimension: :c:expr:`const int*`
   :param nodes:
      **[Input]** The actual control net of the B |eacute| zier triangle as a
      :math:`D \times N` array. This should be laid out in Fortran order, with
      :math:`D N` total values.
   :type nodes: :c:expr:`const double*`
   :param degree:
      **[Input]** The degree :math:`d` of the B |eacute| zier triangle.
   :type degree: :c:expr:`const int*`
   :param nodes_a:
      **[Output]** The control net of the lower left sub-triangle as a
      :math:`D \times N` array, laid out in Fortran order.
   :type nodes_a: :c:expr:`double*`
   :param nodes_b:
      **[Output]** The control net of the central sub-triangle as a
      :math:`D \times N` array, laid out in Fortran order.
   :type nodes_b: :c:expr:`double*`
   :param nodes_c:
      **[Output]** The control net of the lower right sub-triangle as a
      :math:`D \times N` array, laid out in Fortran order.
   :type nodes_c: :c:expr:`double*`
   :param nodes_d:
      **[Output]** The control net of the upper left sub-triangle as a
      :math:`D \times N` array, laid out in Fortran order.
   :type nodes_d: :c:expr:`double*`

   **Signature:**

   .. code-block:: c

      void
      BEZ_subdivide_nodes_triangle(const int *num_nodes,
                                   const int *dimension,
                                   const double *nodes,
                                   const int *degree,
                                   double *nodes_a,
                                   double *nodes_b,
                                   double *nodes_c,
                                   double *nodes_d);
