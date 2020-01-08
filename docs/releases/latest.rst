Latest Release (``0.10.1.dev1``)
================================

|pypi| |docs|

Bug Fixes
---------

-  Fixing high-degree error in ``Curve.evaluate()``
   `method <https://bezier.readthedocs.io/en/latest/python/reference/bezier.curve.html#bezier.curve.Curve.evaluate>`__,
   via the ``evaluate_curve_barycentric()`` Fortran
   `subroutine <https://bezier.readthedocs.io/en/latest/abi/curve.html#c.evaluate_curve_barycentric>`__
   (`5768824 <https://github.com/dhermes/bezier/commit/57688243b9264ca7ea48423f100e8f516ba2fa2f>`__).
   Fixed `#156 <https://github.com/dhermes/bezier/issues/156>`__. The code uses
   :math:`\binom{n}{k + 1} = \frac{n - k}{k + 1} \binom{n}{k}` to update the
   value and :math:`(30 - 14) \binom{30}{14}` caused the 32-bit signed integer
   to overflow.

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.10.1.svg
   :target: https://pypi.org/project/bezier/0.10.1/
   :alt: PyPI link to release 0.10.1
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.10.1
   :target: https://bezier.readthedocs.io/en/0.10.1/
   :alt: Documentation for release 0.10.1
