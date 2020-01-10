############
curve module
############

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

This is a collection of procedures for performing computations on a
B |eacute| zier `curve`_.

.. _curve: https://en.wikipedia.org/wiki/B%C3%A9zier_curve

**********
Procedures
**********

.. c:function:: void BEZ_compute_length(int *num_nodes, \
                                        int *dimension, \
                                        double *nodes, \
                                        double *length, \
                                        int *error_val)

   Computes the length of a B |eacute| zier curve via

   .. math::

      \ell = \int_0^1 \left\lVert B'(s) \right\rVert \, ds.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param double* length:
      **[Output]** The computed length :math:`\ell`.
   :param int* error_val:
      **[Output]** An error status passed along from ``dqagse`` (a QUADPACK
      procedure).

   **Signature:**

   .. code-block:: c

      void
      BEZ_compute_length(int *num_nodes,
                         int *dimension,
                         double *nodes,
                         double *length,
                         int *error_val);

   **Example:**

   .. literalinclude:: example_compute_length.c
      :language: c
      :dedent: 4
      :lines: 18-29

   Consider the line segment :math:`B(s) = \left[\begin{array}{c} 3s \\ 4s
   \end{array}\right]`, we can verify the length:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_compute_length.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Length: 5.000000
      Error value: 0

.. c:function:: void BEZ_elevate_nodes_curve(int *num_nodes, \
                                             int *dimension, \
                                             double *nodes, \
                                             double *elevated)

   Degree-elevate a B |eacute| zier curve. Does so by producing
   control points of a higher degree that define the exact same curve.

   See :meth:`.Curve.elevate` for more details.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param double* elevated:
      **[Output]** The control points of the degree-elevated curve as a
      :math:`D \times (N + 1)` array, laid out in Fortran order.

   **Signature:**

   .. code-block:: c

      void
      BEZ_elevate_nodes_curve(int *num_nodes,
                              int *dimension,
                              double *nodes,
                              double *elevated);

   **Example:**

   After elevating :math:`B(s) = \left[\begin{array}{c} 0 \\ 0
   \end{array}\right] (1 - s)^2 + \frac{1}{2} \left[\begin{array}{c} 3 \\ 3
   \end{array}\right] 2 (1 - s) s + \left[\begin{array}{c} 3 \\ 0
   \end{array}\right] s^2`:

   .. literalinclude:: example_elevate_nodes_curve.c
      :language: c
      :dedent: 4
      :lines: 18-31

   we have :math:`B(s) = \left[\begin{array}{c} 0 \\ 0
   \end{array}\right] (1 - s)^3 + \left[\begin{array}{c} 1 \\ 1
   \end{array}\right] 3 (1 - s)^2 s + \left[\begin{array}{c} 2 \\ 1
   \end{array}\right] 3 (1 - s) s^2 + \left[\begin{array}{c} 3 \\ 0
   \end{array}\right] s^3`:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_elevate_nodes_curve.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Elevated:
      0.000000, 1.000000, 2.000000, 3.000000
      0.000000, 1.000000, 1.000000, 0.000000

   .. image:: ../images/curve_elevate.png
      :align: center

.. c:function:: void BEZ_evaluate_curve_barycentric(int *num_nodes, \
                                                    int *dimension, \
                                                    double *nodes, \
                                                    int *num_vals, \
                                                    double *lambda1, \
                                                    double *lambda2, \
                                                    double *evaluated)

   For a B |eacute| zier curve with control points :math:`p_0, \ldots, p_d`,
   this evaluates the quantity

   .. math::

      Q(\lambda_1, \lambda_2) =
          \sum_{j = 0}^d \binom{d}{j} \lambda_1^{d - j} \lambda_2^j p_j.

   The typical case is barycentric, i.e. :math:`\lambda_1 + \lambda_2 = 1`, but
   this is not required.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param int* num_vals:
      **[Input]** The number of values :math:`k` where the quantity will be
      evaluated.
   :param double* lambda1:
      **[Input]** An array of :math:`k` values used for the first parameter
      :math:`\lambda_1`.
   :param double* lambda2:
      **[Input]** An array of :math:`k` values used for the second parameter
      :math:`\lambda_2`.
   :param double* evaluated:
      **[Output]** The evaluated quantites as a :math:`D \times k` array, laid
      out in Fortran order. Column :math:`j` of ``evaluated`` will contain
      :math:`Q\left(\lambda_1\left[j\right], \lambda_2\left[j\right]\right)`.

   **Signature:**

   .. code-block:: c

      void
      BEZ_evaluate_curve_barycentric(int *num_nodes,
                                     int *dimension,
                                     double *nodes,
                                     int *num_vals,
                                     double *lambda1,
                                     double *lambda2,
                                     double *evaluated);

   **Example:**

   For the curve :math:`B(s) = \left[\begin{array}{c} 0 \\ 1
   \end{array}\right] (1 - s)^2 + \left[\begin{array}{c} 2 \\ 1
   \end{array}\right] 2 (1 - s) s + \left[\begin{array}{c} 3 \\ 3
   \end{array}\right] s^2 = \left[\begin{array}{c} s(4 - s) \\ 2s^2 + 1
   \end{array}\right]`:

   .. literalinclude:: example_evaluate_curve_barycentric.c
      :language: c
      :dedent: 4
      :lines: 18-35

   we have

   .. math::

      \begin{align*}
        Q\left(\frac{1}{4}, \frac{3}{4}\right) &= \frac{1}{16} \left[
          \begin{array}{c} 39 \\ 34 \end{array}\right] \\
        Q\left(\frac{1}{2}, \frac{1}{4}\right) &= \frac{1}{16} \left[
          \begin{array}{c} 11 \\ 11 \end{array}\right] \\
        Q\left(0, \frac{1}{2}\right) &= \frac{1}{4} \left[
          \begin{array}{c} 3 \\ 3 \end{array}\right] \\
        Q\left(1, \frac{1}{4}\right) &= \frac{1}{16} \left[
          \begin{array}{c} 19 \\ 27 \end{array}\right]
      \end{align*}

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_evaluate_curve_barycentric.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Evaluated:
      2.437500, 0.687500, 0.750000, 1.187500
      2.125000, 0.687500, 0.750000, 1.687500

.. c:function:: void BEZ_evaluate_hodograph(double *s, \
                                            int *num_nodes, \
                                            int *dimension, \
                                            double *nodes, \
                                            double *hodograph)

   Evaluates the hodograph (or derivative) of a B |eacute| zier curve
   function :math:`B'(s)`.

   :param double* s:
      **[Input]** The parameter :math:`s` where the hodograph is being
      computed.
   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param double* hodograph:
      **[Output]** The hodograph :math:`B'(s)` as a :math:`D \times 1` array.

   **Signature:**

   .. code-block:: c

      void
      BEZ_evaluate_hodograph(double *s,
                             int *num_nodes,
                             int *dimension,
                             double *nodes,
                             double *hodograph);

   **Example:**

   For the curve :math:`B(s) = \left[\begin{array}{c} 1 \\ 0
   \end{array}\right] (1 - s)^3 + \left[\begin{array}{c} 1 \\ 1
   \end{array}\right] 3 (1 - s)^2 s + \left[\begin{array}{c} 2 \\ 0
   \end{array}\right] 3 (1 - s) s^2 + \left[\begin{array}{c} 2 \\ 1
   \end{array}\right] s^3`:

   .. literalinclude:: example_evaluate_hodograph.c
      :language: c
      :dedent: 4
      :lines: 18-28

   we have :math:`B'\left(\frac{1}{8}\right) = \frac{1}{32} \left[
   \begin{array}{c} 21 \\ 54 \end{array}\right]`:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_evaluate_hodograph.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Hodograph:
      0.656250
      1.687500

.. c:function:: void BEZ_evaluate_multi(int *num_nodes, \
                                        int *dimension, \
                                        double *nodes, \
                                        int *num_vals, \
                                        double *s_vals, \
                                        double *evaluated)

   Evaluate a B |eacute| zier curve function :math:`B(s_j)` at
   multiple values :math:`\left\{s_j\right\}_j`.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param int* num_vals:
      **[Input]** The number of values :math:`k` where the :math:`B(s)` will be
      evaluated.
   :param double* s_vals:
      **[Input]** An array of :math:`k` values :math:`s_j`.
   :param double* evaluated:
      **[Output]** The evaluated points as a :math:`D \times k` array, laid
      out in Fortran order. Column :math:`j` of ``evaluated`` will contain
      :math:`B\left(s_j\right)`.

   **Signature:**

   .. code-block:: c

      void
      BEZ_evaluate_multi(int *num_nodes,
                         int *dimension,
                         double *nodes,
                         int *num_vals,
                         double *s_vals,
                         double *evaluated);

   **Example:**

   For the curve :math:`B(s) = \left[\begin{array}{c} 1 \\ 0
   \end{array}\right] (1 - s)^3 + \left[\begin{array}{c} 1 \\ 1
   \end{array}\right] 3 (1 - s)^2 s + \left[\begin{array}{c} 2 \\ 0
   \end{array}\right] 3 (1 - s) s^2 + \left[\begin{array}{c} 2 \\ 1
   \end{array}\right] s^3`:

   .. literalinclude:: example_evaluate_multi.c
      :language: c
      :dedent: 4
      :lines: 18-32

   we have :math:`B\left(0\right) = \left[\begin{array}{c}
   1 \\ 0 \end{array}\right], B\left(\frac{1}{2}\right) =
   \frac{1}{2} \left[\begin{array}{c} 3 \\ 1 \end{array}\right]` and
   :math:`B\left(1\right) = \left[\begin{array}{c} 2 \\ 1 \end{array}\right]`:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_evaluate_multi.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Evaluated:
      1.000000, 1.500000, 2.000000
      0.000000, 0.500000, 1.000000

.. c:function:: void BEZ_full_reduce(int *num_nodes, \
                                     int *dimension, \
                                     double *nodes, \
                                     int *num_reduced_nodes, \
                                     double *reduced, \
                                     bool *not_implemented)

   Perform a "full" degree reduction. Does so by using
   :c:func:`BEZ_reduce_pseudo_inverse` continually until the degree of
   the curve can no longer be reduced.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param int* num_reduced_nodes:
      **[Output]** The number of control points :math:`M` of the fully reduced
      curve.
   :param double* reduced:
      **[Output]** The control points of the fully reduced curve as a
      :math:`D \times N` array. The first :math:`M` columns will contain the
      reduced nodes. ``reduced`` must be allocated by the caller and since
      :math:`M = N` occurs when no reduction is possible, this array must be
      :math:`D \times N`.
   :param bool* not_implemented:
      **[Output]** Indicates if degree-reduction has been implemented for the
      current curve's degree. (Currently, the only degrees supported are 1,
      2, 3 and  4.)

   **Signature:**

   .. code-block:: c

      void
      BEZ_full_reduce(int *num_nodes,
                      int *dimension,
                      double *nodes,
                      int *num_reduced_nodes,
                      double *reduced,
                      bool *not_implemented);

   **Example:**

   When taking a curve that is degree-elevated from linear to quartic:

   .. literalinclude:: example_full_reduce.c
      :language: c
      :dedent: 4
      :lines: 18-34

   this procedure reduces it to the line
   :math:`B(s) = \left[\begin{array}{c} 1 \\ 3
   \end{array}\right] (1 - s) + \left[\begin{array}{c} 2 \\ 5
   \end{array}\right] s = \left[\begin{array}{c} 1 + s \\ 3 + 2s
   \end{array}\right]`:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_full_reduce.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Number of reduced nodes: 2
      Reduced:
      1.000000, 2.000000
      3.000000, 5.000000
      Not implemented: FALSE

.. c:function:: void BEZ_get_curvature(int *num_nodes, \
                                       double *nodes, \
                                       double *tangent_vec, \
                                       double *s, \
                                       double *curvature)

   Get the signed curvature of a B |eacute| zier curve at a point. See
   :func:`._py_curve_helpers.get_curvature` for more details.

   .. note::

      This **only** computes curvature for plane curves (i.e. curves
      in :math:`\mathbf{R}^2`). An equivalent notion of curvature exists for
      space curves, but support for that is not implemented here.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`2 \times N` array. This should be laid out in Fortran order,
      with :math:`2 N` total values.
   :param double* tangent_vec:
      **[Input]** The hodograph :math:`B'(s)` as a :math:`2 \times 1` array.
      Note that this could be computed once :math:`s` and :math:`B` are known,
      but this allows the caller to re-use an already computed tangent vector.
   :param double* s:
      **[Input]** The parameter :math:`s` where the curvature is being
      computed.
   :param double* curvature:
      **[Output]** The signed curvature :math:`\kappa`.

   **Signature:**

   .. code-block:: c

      void
      BEZ_get_curvature(int *num_nodes,
                        double *nodes,
                        double *tangent_vec,
                        double *s,
                        double *curvature);

   **Example:**

   .. literalinclude:: example_get_curvature.c
      :language: c
      :dedent: 4
      :lines: 18-28

   .. image:: ../images/get_curvature.png
      :align: center

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_get_curvature.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Curvature: -12.000000

.. c:function:: void BEZ_locate_point_curve(int *num_nodes, \
                                            int *dimension, \
                                            double *nodes, \
                                            double *point, \
                                            double *s_approx)

   This solves the inverse problem :math:`B(s) = p` (if it can be
   solved). Does so by subdividing the curve until the segments are
   sufficiently small, then using Newton's method to narrow in on the
   pre-image of the point.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param double* point:
      **[Input]** The point :math:`p` as a :math:`D \times 1` array.
   :param double* s_approx:
      **[Output]** The parameter :math:`s` of the solution. If
      :math:`p` can't be located on the curve, then ``s_approx = -1.0``.
      If there are **multiple** parameters that satisfy :math:`B(s) = p`
      (indicating that :math:`B(s)` has a self-crossing) then
      ``s_approx = -2.0``.

   **Signature:**

   .. code-block:: c

      void
      BEZ_locate_point_curve(int *num_nodes,
                             int *dimension,
                             double *nodes,
                             double *point,
                             double *s_approx);

   **Example:**

   For :math:`B(s) = \left[\begin{array}{c} 0 \\ 2
   \end{array}\right] (1 - s)^3 + \left[\begin{array}{c} -1 \\ 0
   \end{array}\right] 3 (1 - s)^2 s + \left[\begin{array}{c} 1 \\ 1
   \end{array}\right] 3 (1 - s) s^2 + \frac{1}{8} \left[\begin{array}{c}
   -6 \\ 13 \end{array}\right] s^3`:

   .. literalinclude:: example_locate_point_curve.c
      :language: c
      :dedent: 4
      :lines: 18-34

   We can locate the point :math:`B\left(\frac{1}{2}\right) = \frac{1}{64}
   \left[\begin{array}{c} -6 \\ 53 \end{array}\right]` but find that
   :math:`\frac{1}{2} \left[\begin{array}{c} 0 \\ 3 \end{array}\right]` is
   not on the curve and that

   .. math::

      B\left(\frac{3 - \sqrt{5}}{6}\right) =
          B\left(\frac{3 + \sqrt{5}}{6}\right) = \frac{1}{8} \left[
          \begin{array}{c} -2 \\ 11 \end{array}\right]

   is a self-crossing:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_locate_point_curve.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      When B(s) = [-0.093750, 0.828125]; s =  0.500000
      When B(s) = [ 0.000000, 1.500000]; s = -1.000000
      When B(s) = [-0.250000, 1.375000]; s = -2.000000

   .. image:: ../images/curve_locate.png
      :align: center

.. c:function:: void BEZ_newton_refine_curve(int *num_nodes, \
                                             int *dimension, \
                                             double *nodes, \
                                             double *point, \
                                             double *s, \
                                             double *updated_s)

   This refines a solution to :math:`B(s) = p` using Newton's
   method. Given a current approximation :math:`s_n` for a solution,
   this produces the updated approximation via

   .. math::

      s_{n + 1} = s_n - \frac{B'(s_n)^T \left[B(s_n) - p\right]}{
          B'(s_n)^T B'(s_n)}.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param double* point:
      **[Input]** The point :math:`p` as a :math:`D \times 1` array.
   :param double* s:
      **[Input]** The parameter :math:`s_n` of the current approximation
      of a solution.
   :param double* updated_s:
      **[Output]** The parameter :math:`s_{n + 1}` of the updated
      approximation.

   **Signature:**

   .. code-block:: c

      void
      BEZ_newton_refine_curve(int *num_nodes,
                              int *dimension,
                              double *nodes,
                              double *point,
                              double *s,
                              double *updated_s);

   **Example:**

   When trying to locate :math:`B\left(\frac{1}{4}\right) = \frac{1}{16}
   \left[\begin{array}{c} 9 \\ 13 \end{array}\right]` on the curve
   :math:`B(s) = \left[\begin{array}{c} 0 \\ 0
   \end{array}\right] (1 - s)^2 + \left[\begin{array}{c} 1 \\ 2
   \end{array}\right] 2 (1 - s) s + \left[\begin{array}{c} 3 \\ 1
   \end{array}\right] s^2`, starting at :math:`s = \frac{3}{4}`:

   .. literalinclude:: example_newton_refine_curve.c
      :language: c
      :dedent: 4
      :lines: 18-29

   we expect a Newton update :math:`\Delta s = -\frac{2}{5}`, which produces
   a new parameter value :math:`s = \frac{7}{20}`:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_newton_refine_curve.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Updated s: 0.350000

   .. image:: ../images/newton_refine_curve.png
      :align: center

.. c:function:: void BEZ_reduce_pseudo_inverse(int *num_nodes, \
                                               int *dimension, \
                                               double *nodes, \
                                               double *reduced, \
                                               bool *not_implemented)

   Perform a pseudo inverse to :c:func:`BEZ_elevate_nodes_curve`. If an
   inverse can be found, i.e. if a curve can be degree-reduced, then
   this will produce the equivalent curve of lower degree. If no
   inverse can be found, then this will produce the "best" answer in
   the least squares sense.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param double* reduced:
      **[Output]** The control points of the degree-(pseudo)reduced curve
      :math:`D \times (N - 1)` array, laid out in Fortran order.
   :param bool* not_implemented:
      **[Output]** Indicates if degree-reduction has been implemented for the
      current curve's degree. (Currently, the only degrees supported are 1,
      2, 3 and  4.)

   **Signature:**

   .. code-block:: c

      void
      BEZ_reduce_pseudo_inverse(int *num_nodes,
                                int *dimension,
                                double *nodes,
                                double *reduced,
                                bool *not_implemented);

   **Example:**

   After reducing :math:`B(s) = \left[\begin{array}{c} -3 \\ 3
   \end{array}\right] (1 - s)^3 + \left[\begin{array}{c} 0 \\ 2
   \end{array}\right] 3 (1 - s)^2 s + \left[\begin{array}{c} 1 \\ 3
   \end{array}\right] 3 (1 - s) s^2 + \left[\begin{array}{c} 0 \\ 6
   \end{array}\right] s^3`:

   .. literalinclude:: example_reduce_pseudo_inverse.c
      :language: c
      :dedent: 4
      :lines: 18-32

   we get the valid quadratic representation of :math:`B(s) =
   \left[\begin{array}{c} 3(1 - s)(2s - 1) \\ 3(2s^2 - s + 1)
   \end{array}\right]`:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_reduce_pseudo_inverse.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Reduced:
      -3.000000, 1.500000, 0.000000
       3.000000, 1.500000, 6.000000
      Not implemented: FALSE

   .. image:: ../images/curve_reduce.png
      :align: center

.. c:function:: void BEZ_specialize_curve(int *num_nodes, \
                                          int *dimension, \
                                          double *nodes, \
                                          double *start, \
                                          double *end, \
                                          double *new_nodes)

   Specialize a B |eacute| zier curve to an interval
   :math:`\left[a, b\right]`. This produces the control points
   for the curve given by :math:`B\left(a + (b - a) s\right)`.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param double* start:
      **[Input]** The start :math:`a` of the specialized interval.
   :param double* end:
      **[Input]** The end :math:`b` of the specialized interval.
   :param double* new_nodes:
      **[Output]** The control points of the specialized curve, as a
      :math:`D \times N` array, laid out in Fortran order.

   **Signature:**

   .. code-block:: c

      void
      BEZ_specialize_curve(int *num_nodes,
                           int *dimension,
                           double *nodes,
                           double *start,
                           double *end,
                           double *new_nodes);

   **Example:**

   When we specialize the curve :math:`B(s) = \left[\begin{array}{c} 0 \\ 0
   \end{array}\right] (1 - s)^2 + \frac{1}{2} \left[\begin{array}{c} 1 \\ 2
   \end{array}\right] 2 (1 - s) s + \left[\begin{array}{c} 1 \\ 0
   \end{array}\right] s^2 = \left[\begin{array}{c} s \\ 2s(1 - s)
   \end{array}\right]` to the interval :math:`\left[-\frac{1}{4},
   \frac{3}{4}\right]`:

   .. literalinclude:: example_specialize_curve.c
      :language: c
      :dedent: 4
      :lines: 18-31

   we get the specialized curve :math:`S(t) = \frac{1}{8} \left[
   \begin{array}{c} -2 \\ -5 \end{array}\right] (1 - s)^2 + \frac{1}{8}
   \left[\begin{array}{c} 2 \\ 7 \end{array}\right] 2 (1 - s) s + \frac{1}{8}
   \left[\begin{array}{c} 6 \\ 3 \end{array}\right] s^2 = \frac{1}{8}
   \left[\begin{array}{c} 2(4t - 1) \\ (4t - 1)(5 - 4t) \end{array}\right]`,
   which still lies on :math:`y = 2x(1 - x)`:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_specialize_curve.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      New Nodes:
      -0.250000, 0.250000, 0.750000
      -0.625000, 0.875000, 0.375000

   .. image:: ../images/curve_specialize.png
      :align: center

.. c:function:: void BEZ_subdivide_nodes_curve(int *num_nodes, \
                                               int *dimension, \
                                               double *nodes, \
                                               double *left_nodes, \
                                               double *right_nodes)

   Split a B |eacute| zier curve into two halves
   :math:`B\left(\left[0, \frac{1}{2}\right]\right)` and
   :math:`B\left(\left[\frac{1}{2}, 1\right]\right)`.

   :param int* num_nodes:
      **[Input]** The number of control points :math:`N` of a
      B |eacute| zier curve.
   :param int* dimension:
      **[Input]** The dimension :math:`D` such that the curve lies in
      :math:`\mathbf{R}^D`.
   :param double* nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :param double* left_nodes:
      **[Output]** The control points of the left half curve
      :math:`B\left(\left[0, \frac{1}{2}\right]\right)` as a
      :math:`D \times N` array, laid out in Fortran order.
   :param double* right_nodes:
      **[Output]** The control points of the right half curve
      :math:`B\left(\left[\frac{1}{2}, 1\right]\right)` as a
      :math:`D \times N` array, laid out in Fortran order.

   **Signature:**

   .. code-block:: c

      void
      BEZ_subdivide_nodes_curve(int *num_nodes,
                                int *dimension,
                                double *nodes,
                                double *left_nodes,
                                double *right_nodes);

   **Example:**

   For example, subdividing the curve :math:`B(s) = \left[\begin{array}{c}
   0 \\ 0 \end{array}\right] (1 - s)^2 + \frac{1}{4} \left[\begin{array}{c}
   5 \\ 12 \end{array}\right] 2 (1 - s) s + \left[\begin{array}{c} 2 \\ 1
   \end{array}\right] s^2`:

   .. literalinclude:: example_subdivide_nodes_curve.c
      :language: c
      :dedent: 4
      :lines: 18-34

   yields:

   .. code-block:: console

      $ gcc \
      >   -o example \
      >   example_subdivide_nodes_curve.c \
      >   -I .../src/fortran/include \
      >   -L .../site-packages/bezier/lib \
      >   -lbezier \
      >   -lm -lgfortran
      $ ./example
      Left Nodes:
      0.000000, 0.625000, 1.125000
      0.000000, 1.500000, 1.750000
      Right Nodes:
      1.125000, 1.625000, 2.000000
      1.750000, 2.000000, 1.000000

   .. image:: ../images/curve_subdivide.png
      :align: center
