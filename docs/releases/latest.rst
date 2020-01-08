Latest Release (``0.10.1.dev1``)
================================

|pypi| |docs|

Python Changes
--------------

Packaging
~~~~~~~~~

-  Explicit support for Python 3.8 has been added
   (`#161 <https://github.com/dhermes/bezier/issues/161>`__).

New Features
~~~~~~~~~~~~

-  Loosening type constraints in ``Curve``
   `constructor <https://bezier.readthedocs.io/en/latest/python/reference/bezier.curve.html#bezier.curve.Curve>`__
   and ``Surface``
   `constructor <https://bezier.readthedocs.io/en/latest/python/reference/bezier.surface.html#bezier.surface.Surface>`__;
   now any sequence type is accepted rather than **only** NumPy arrays
   (`68f7dc7 <https://github.com/dhermes/bezier/commit/68f7dc7c1f26bb678d09b4221fd917531fb79860>`__,
   `a8c68a3 <https://github.com/dhermes/bezier/commit/a8c68a3368a1edf90cd76cd6ff77ab698b6c3907>`__,
   `f5c7869 <https://github.com/dhermes/bezier/commit/f5c7869e86b196aca3db272a2e85413357864bc7>`__).
   Fixed `#146 <https://github.com/dhermes/bezier/issues/146>`__.

Internals
~~~~~~~~~

-  Re-factored non-public modules so that algorithms implemented in pure Python
   **only** invoke other algorithms written in pure Python
   (`#160 <https://github.com/dhermes/bezier/pull/160>`__). Previously
   these algorithms invoked the equivalent Fortran speedup if present for a
   given function. Fixed
   `#159 <https://github.com/dhermes/bezier/issues/159>`__.

Miscellany
~~~~~~~~~~

-  Moving ``*.f90`` Fortran files **out** of Python source tree
   (`#152 <https://github.com/dhermes/bezier/pull/152>`__).

Bug Fixes
---------

-  Explicitly handling length 0 curves (with an error) in the
   ``compute_length()`` Fortran
   `subroutine <https://bezier.readthedocs.io/en/latest/abi/curve.html#c.compute_length>`__
   that is used by the ``Curve.length``
   `property <https://bezier.readthedocs.io/en/latest/python/reference/bezier.curve.html#bezier.curve.Curve.length>`__
   (`a24368f <https://github.com/dhermes/bezier/commit/a24368fc690b2c6d6a676b9d569f25b5919c400d>`__).
   Fixed `#148 <https://github.com/dhermes/bezier/issues/148>`__.
-  Fixing high-degree error in ``Curve.evaluate()``
   `method <https://bezier.readthedocs.io/en/latest/python/reference/bezier.curve.html#bezier.curve.Curve.evaluate>`__,
   via the ``evaluate_curve_barycentric()`` Fortran
   `subroutine <https://bezier.readthedocs.io/en/latest/abi/curve.html#c.evaluate_curve_barycentric>`__
   (`5768824 <https://github.com/dhermes/bezier/commit/57688243b9264ca7ea48423f100e8f516ba2fa2f>`__).
   Fixed `#156 <https://github.com/dhermes/bezier/issues/156>`__. The code uses
   :math:`\binom{n}{k + 1} = \frac{n - k}{k + 1} \binom{n}{k}` to update the
   value and :math:`(30 - 14) \binom{30}{14}` overflows a 32-bit signed
   integer.

Documentation
-------------

-  Updating install instructions to show how to disable the binary extension
   via ``BEZIER_NO_EXTENSIONS``
   (`6262594 <https://github.com/dhermes/bezier/commit/626259493997a9d83924d100900189f32b87e6c5>`__).
   Fixed `#147 <https://github.com/dhermes/bezier/issues/147>`__.
-  Adding "Citation" section to landing page
   (`9885063 <https://github.com/dhermes/bezier/commit/9885063a2e3795e0bec35a4fc1574dc294d359e0>`__).
   Fixed `#150 <https://github.com/dhermes/bezier/issues/150>`__.

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.10.1.svg
   :target: https://pypi.org/project/bezier/0.10.1/
   :alt: PyPI link to release 0.10.1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.10.1
   :target: https://bezier.readthedocs.io/en/0.10.1/
   :alt: Documentation for release 0.10.1
