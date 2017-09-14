Algorithm Helpers
=================

In an attempt to thoroughly vet each algorithm used in this
library, each computation is split into small units that can
be tested independently.

Though many of these computational units aren't provided as part
of the public interface of this library, they are still interesting.
(Possibly) more importantly, it's useful to see these algorithms
at work.

In this document, these helper functions and objects are documented.
This is to help with the exposition of the computation and
**does not** imply that these are part of the stable public interface.

.. autofunction:: bezier._intersection_helpers._linearization_error
.. autofunction:: bezier._intersection_helpers._newton_refine
.. autofunction:: bezier._intersection_helpers._segment_intersection
.. autofunction:: bezier._intersection_helpers._parallel_different
.. autofunction:: bezier._curve_helpers._get_curvature
.. autofunction:: bezier._curve_helpers._newton_refine
.. autoclass:: bezier._intersection_helpers.Intersection
   :members:
.. autoclass:: bezier._intersection_helpers.Linearization
   :members:
.. autoclass:: bezier._surface_helpers.IntersectionClassification
   :members:
.. autofunction:: bezier._surface_helpers.newton_refine
.. autofunction:: bezier._surface_helpers.classify_intersection
.. autofunction:: bezier._surface_helpers._jacobian_det
.. autofunction:: bezier._surface_helpers.two_by_two_det
.. autofunction:: bezier._implicitization.bezier_roots
.. autofunction:: bezier._implicitization.lu_companion
.. autofunction:: bezier._implicitization.bezier_value_check

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
