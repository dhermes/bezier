#############
status module
#############

.. c:type:: Status

   This enum contains all status codes that can be returned from
   procedures that can fail for **some** inputs.

   .. c:var:: SUCCESS

      (``0``)
      Procedure exited with no error.

   .. c:var:: BAD_MULTIPLICITY

      (``1``)
      An iterative method (e.g. Newton's method) failed to converge because it
      encountered a solution with an unsupported multiplicity (i.e. it was a
      triple root or higher). The multiplicity is detected by the observed
      rate of convergence.

   .. c:var:: NO_CONVERGE

      (``2``)
      An iterative method failed to converge in the maximum allowed
      number of iterations.

   .. c:var:: INSUFFICIENT_SPACE

      (``3``)
      All allocation must be done by the caller, hence any variable-size
      output must be already allocated. This will be used when the caller has
      not allocated enough space for the output values.

   .. c:var:: SAME_CURVATURE

      (``4``)
      Classification of a curve-curve intersection failed due to two curves
      having identical curvature at a tangent intersection.

   .. c:var:: BAD_INTERIOR

      (``5``)
      Caused by a failure during the process of triangle-triangle intersection.
      Occurs when the corners and edge-edge intersections can't be converted
      into the curved polygon(s) that make up the intersection of the two
      triangles.

   .. c:var:: EDGE_END

      (``6``)
      A triangle-triangle intersection point occurs at the **end** of an
      edge (only intersections at the beginning of an edge should be used).

   .. c:var:: SINGULAR

      (``7``)
      Signifies that an attempt was made to solve a linear system that was
      singular.

   .. c:var:: UNKNOWN

      (``999``)
      Signifies a block of code reached an "impossible" state. Either the
      code has violated some mathematical invariant or the author of that
      block of code has misunderstood the possible states of the system.

   **Example:**

   Consider curves which intersect at the point
   :math:`B_1\left(\frac{3}{2}\right) = B_2(-2)`:

   .. math::

      \begin{align*}
      B_1(s) &= \left[\begin{array}{c} 0 \\ 0 \end{array}\right] (1 - s)^2
          + \left[\begin{array}{c} 2 \\ 4 \end{array}\right] 2s(1 - s)
          + \left[\begin{array}{c} 4 \\ 0 \end{array}\right] s^2
          = \left[\begin{array}{c} 4s \\ 8s(1 - s) \end{array}\right] \\
      B_2(t) &= \left[\begin{array}{c} 2 \\ 0 \end{array}\right] (1 - t)
          + \left[\begin{array}{c} 0 \\ 3 \end{array}\right] t
          = \left[\begin{array}{c} 2(1 - t) \\ 3t \end{array}\right]
          \end{align*}

   When trying to use Newton's method to find the root of

   .. math::

      F(s, t) = B_1(s) - B_2(t)

   nearest to :math:`s = \frac{7}{8}, t = -2`:

   .. literalinclude:: example_status.c
      :language: c
      :dedent: 4
      :lines: 18-34

   the method fails with the status :c:data:`SINGULAR` because the
   Jacobian

   .. math::

      DF\left(\frac{7}{8}, -2\right) = \left[\begin{array}{c c}
       4 &  2 \\
      -6 & -3
      \end{array}\right]

   is singular to numerical precision:

   .. testsetup:: example-status

      import tests.utils


      build_and_run_c = tests.utils.build_and_run_c

   .. testcode:: example-status
      :hide:

      build_and_run_c("example_status.c")

   .. testoutput:: example-status
      :options: +NORMALIZE_WHITESPACE
      :windows-skip:

      $ INCLUDE_DIR=.../libbezier-release/usr/include
      $ LIB_DIR=.../libbezier-release/usr/lib
      $ gcc \
      >     -o example \
      >     example_status.c \
      >     -I "${INCLUDE_DIR}" \
      >     -L "${LIB_DIR}" \
      >     -Wl,-rpath,"${LIB_DIR}" \
      >     -lbezier \
      >     -lm -lgfortran
      $ ./example
      Jacobian is singular.
