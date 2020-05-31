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

.. cpp:namespace:: bezier

.. cpp:function:: template <std::size_t N, std::size_t D> \
                  std::tuple<double, int> compute_length(const Matrix<D, N>& nodes)

   Computes the length of a B |eacute| zier curve via

   .. math::

      \ell = \int_0^1 \left\lVert B'(s) \right\rVert \, ds.

   :param nodes:
      **[Input]** The actual control points of the curve as a
      :math:`D \times N` array. This should be laid out in Fortran order,
      with :math:`D N` total values.
   :returns:
      A pair (``tuple``) of values

      * The computed length :math:`\ell`.
      * An error status passed along from ``dqagse`` (a QUADPACK procedure).

   **Example:**

   .. literalinclude:: example_compute_length.cpp
      :language: cpp
      :dedent: 4
      :lines: 20-27

   Consider the line segment :math:`B(s) = \left[\begin{array}{c} 3s \\ 4s
   \end{array}\right]`, we can verify the length:

   .. code-block:: console

      $ XTENSOR_INCLUDE_DIR=.../xtensor-release/usr/include
      $ XTL_INCLUDE_DIR=.../xtl-release/usr/include
      $ INCLUDE_DIR=.../libbezier-release/usr/include
      $ LIB_DIR=.../libbezier-release/usr/lib
      $ g++ \
      >     -o example \
      >     example_compute_length.cpp \
      >     -I "${INCLUDE_DIR}" \
      >     -I "${XTENSOR_INCLUDE_DIR}" \
      >     -I "${XTL_INCLUDE_DIR}" \
      >     -L "${LIB_DIR}" \
      >     -Wl,-rpath,"${LIB_DIR}" \
      >     -lbezier
      $ ./example
      Length: 5
      Error value: 0
