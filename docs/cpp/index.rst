#######
C++ API
#######

In addition to the low-level binary ``libbezier`` :doc:`ABI <../abi/index>`,
there is a high-level C++ API. Providing a high-level API for matrices (i.e.
two-dimensional arrays) requires specifying compatible container type(s) as
part of the API. This library (currently) relies on the ``xtensor`` `library`_
to define a ``bezier::Matrix`` type that is contiguous, in Fortran order and
with shape defined at compile time:

.. code-block:: cpp

   template <std::size_t M, std::size_t N>
   using Matrix = xt::xtensor_fixed<double, xt::xshape<M, N>, xt::layout_type::column_major>;

.. _library: https://xtensor.readthedocs.io/en/latest/

The C++ API is just a header-only wrapper around the C ABI; the headers are
:doc:`installed <../abi/installation>` when installing ``libbezier``. The
corresponding header files for the C++ API match the module structure from
the C ABI headers:

.. toctree::
   :titlesonly:

   curve
   helpers
