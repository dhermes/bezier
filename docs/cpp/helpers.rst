##############
helpers module
##############

This is a collection of generic utility procedures that help
with various geometric computations.

*****
Types
*****

.. cpp:namespace:: bezier

.. cpp:type:: template <std::size_t M, std::size_t N> \
              Matrix

   .. code-block:: cpp

      template <std::size_t M, std::size_t N>
      using Matrix = xt::xtensor_fixed<double, xt::xshape<M, N>, xt::layout_type::column_major>;
