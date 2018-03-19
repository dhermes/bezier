Latest Release (``0.8.0``)
==========================

|pypi| |docs|

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

New Features
------------

-  Adding support for surface-surface intersections that have
   coincident segments shared between each surface
   (`cfa2b93 <https://github.com/dhermes/bezier/commit/cfa2b93792695b87f11ece9da1959013ecf77678>`__,
   `0a9645c <https://github.com/dhermes/bezier/commit/0a9645c9a3f1df3274677ad3def3d934c590b642>`__).
   See cases:

   -  4: `10Q-18Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces10Q_and_18Q.png>`__
   -  5: `10Q-19Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces10Q_and_19Q.png>`__
   -  43: `29Q-42Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces29Q_and_42Q.png>`__
   -  44: `29Q-43Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces29Q_and_43Q.png>`__
   -  45: `10Q-44Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces10Q_and_44Q.png>`__
   -  46: `1Q-45Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces1Q_and_45Q.png>`__
   -  47: `1Q-2C <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces1Q_and_2C.png>`__
-  Adding support for curve-curve intersections that are also points of
   tangency. This was accomplished by five broad changes to the geometric
   intersection algorithm:

   -  Checking if almost-linear curves have disjoint bounding boxes
      **before** intersecting the linearized segments
      (`05f0343 <https://github.com/dhermes/bezier/commit/05f0343ca1962dbc5ab3b143b5c6fe20b87272d1>`__).
   -  Adding a "full" Newton iteration for finding ``B1(s) = B2(t)`` when known
      to be near a solution. In particular, this has **special** handling for
      tangencies, which cause a singular Jacobian and make convergence drop
      from quadratic to linear and stalls out convergence early
      (`13a5be5 <https://github.com/dhermes/bezier/commit/13a5be5d80d6a07a1a71326493baa06dbda70f13>`__,
      `4bac61a <https://github.com/dhermes/bezier/commit/4bac61a243b08002c4b0154d2b346cc356097eaf>`__).
   -  Changing how "bad" linearized segments are handled. After subdividing
      to approximately linear curve segments, there were two problems which
      are now being handled in the same way. If the line segments connecting
      the subdivided curve endpoints

      -  are parallel, then the algorithm failed with a ``PARALLEL`` status
      -  intersect outside of the unit interval (for either ``s`` or ``t``),
         the curve-curve candidate was rejected (a small amount, ``0.5^{16}``,
         of "wiggle" room was allowed outside of ``[0, 1]``).

      Now both cases are handled in the same way. First, the subdivided curve
      segments will have a convex hull check applied (which is more strict than
      a bounding box check). If their convex hulls do collide, they are
      treated as a normal intersection of curved segments
      (`4457f64 <https://github.com/dhermes/bezier/commit/4457f64eaf28bb9fb5c91a8740cd0d618fafc3da>`__,
      `fe453c3 <https://github.com/dhermes/bezier/commit/fe453c3839b19ce4a85dfd0b5ad78f71a0973daf>`__).
   -  Using the newly added "full" Newton's iteration for all intersections.
      Before, a single Newton step was applied after intersection the
      linearized segments
      (`d06430f <https://github.com/dhermes/bezier/commit/d06430fbb027eb9d62b6b724f70e62d0efb0732b>`__).
   -  Changing how a candidate pair of ``s-t`` parameters is added.
      (`c998445 <https://github.com/dhermes/bezier/commit/c998445026a5487c59af17c9cbdfc9a6cf4d72c0>`__).
      In the previous implementation, a pair was considered a duplicate
      only if there was a difference of at most 1
      `ULP <https://en.wikipedia.org/wiki/Unit_in_the_last_place>`__ from
      an existing intersection (though this could be toggled via
      ``set_similar_ulps()``). Now, the pair is "normalized" so that ``s``
      and ``t`` are away from ``0``. For example, if ``s < 2^{-10}`` then we
      use ``1 - s`` instead. (This is perfectly "appropriate" since evaluating
      a B |eacute| zier curve requires using both ``s`` and ``1 - s``, so both
      values are equally relevant.) Once normalized, a relative error threshold
      is used.

   Four curve-curve functional test cases have gone from failing to passing:

   -  11: `14-15 <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/curves14_and_15.png>`__
   -  31: `38-39 <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/curves38_and_39.png>`__
   -  43: `58-59 <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/curves58_and_59.png>`__
   -  44: `60-59 <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/curves60_and_59.png>`__

   and two surface-surface cases have as well:

   -  10: `20Q-21Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces20Q_and_21Q.png>`__
   -  42: `41Q-21Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces41Q_and_21Q.png>`__

   In order to support the aforementioned surface-surface cases, special
   support for "tangent corners" was added
   (`12b0de4 <https://github.com/dhermes/bezier/commit/12b0de4e4dae1d84e0681386fd312794ac8736ff>`__).

ABI Changes
-----------

Breaking Changes
~~~~~~~~~~~~~~~~

-  Removed ``BAD_TANGENT`` status enum
   (`b89b2b1 <https://github.com/dhermes/bezier/commit/b89b2b1de1726cdc9f508bd761f4c20e7d655321>`__).
   The case where that failure was issued has now been handled as an acceptable
   ``TANGENT_BOTH`` classification for surface-surface intersection points.
   (See the ``classify_intersection()``
   `function <http://bezier.readthedocs.io/en/0.8.0/algorithm-helpers.html#bezier._surface_helpers.classify_intersection>`__
   for an example.)
-  Adding ``BAD_INTERIOR`` status enum
   (`6348dc6 <https://github.com/dhermes/bezier/commit/6348dc63b5d11453fa8312997429448bbdad0a3f>`__).
   (This is a **breaking** change rather than additive because it re-uses
   the enum value of ``5`` previously used by ``BAD_TANGENT``.) This
   value is used by ``interior_combine()`` in the case that the
   curved polygon intersection(s) cannot be determined from the edge-edge
   intersections for a given surface-surface pair. See
   `#101 <https://github.com/dhermes/bezier/issues/101>`__.
-  Removed ``PARALLEL`` status enum
   (`fe453c3 <https://github.com/dhermes/bezier/commit/fe453c3839b19ce4a85dfd0b5ad78f71a0973daf>`__).
   Now when doing geometric curve-curve intersection, parallel linearized
   segments are handled by checking if the convex hulls collide and then
   (if they do) using a modifed Newton iteration to converge to a root.
-  Adding ``BAD_MULTIPLICITY`` status enum
   (`fe453c3 <https://github.com/dhermes/bezier/commit/fe453c3839b19ce4a85dfd0b5ad78f71a0973daf>`__).
   (This is a **breaking** change rather than additive because it re-uses
   the enum value of ``1`` previously used by ``PARALLEL``.) This is used
   when Newton's method fails to converge to either a simple intersection
   or a tangent intersection. Such failures to converge, when already starting
   near an intersection, may be caused by one of:

   -  The intersection was of multiplicity greater than 2
   -  The curves don't actually intersect, though they come very close
   -  Numerical issues caused the iteration to leave the region
      of convergence
-  Removed ``ulps_away()``
   (`c998445 <https://github.com/dhermes/bezier/commit/c998445026a5487c59af17c9cbdfc9a6cf4d72c0>`__).
-  Removed ``set_similar_ulps()`` and ``get_similar_ulps()``
   (`c998445 <https://github.com/dhermes/bezier/commit/c998445026a5487c59af17c9cbdfc9a6cf4d72c0>`__).

Surface Changes
~~~~~~~~~~~~~~~

-  Added ``SINGULAR`` status enum for cases when a linear system can't be
   solved due to a singular matrix
   (`4457f64 <https://github.com/dhermes/bezier/commit/4457f64eaf28bb9fb5c91a8740cd0d618fafc3da>`__).
-  Adding ``status`` as a return value in ``newton_refine_curve_intersect()``.
   This way, when the Jacobian is singular (which happens at points of
   tangency), the ``SINGULAR`` status can be returned
   (`4457f64 <https://github.com/dhermes/bezier/commit/4457f64eaf28bb9fb5c91a8740cd0d618fafc3da>`__).
   The old behavior would've resulted in a division by zero.

Non-Public API
~~~~~~~~~~~~~~

-  Adding custom linear solver for the ``2 x 2`` case
   (`a3fb476 <https://github.com/dhermes/bezier/commit/a3fb476cf9a82a34754bdd9b9881fbe857883d57>`__).
   This is modelled after ``dgesv`` from LAPACK.

Python Changes
--------------

-  (**Bug fix**) The ``0.7.0`` release broke ``Surface.plot()`` and
   ``CurvedPolygon.plot()`` (when the nodes were transposed, the plotting
   helpers were not correctly updated). The ``add_patch()`` helper was
   fixed to account for the changes in data layout
   (`80bfaaa <https://github.com/dhermes/bezier/commit/80bfaaa74219f9053585aa8970131018baa516d1>`__).
-  Added custom ``UnsupportedDegree``
   `exception <http://bezier.readthedocs.io/en/0.8.0/reference/bezier.html#bezier.UnsupportedDegree>`__
   to be used by routines that have implementations that are hard-coded for
   specific degrees
   (`87a1f21 <https://github.com/dhermes/bezier/commit/87a1f2171f6b810516544ff1691856d7fadfa12f>`__).
   See `#103 <https://github.com/dhermes/bezier/issues/103>`__.
-  Removed ``ulps_away()``
   (`c998445 <https://github.com/dhermes/bezier/commit/c998445026a5487c59af17c9cbdfc9a6cf4d72c0>`__).
-  Removed ``set_similar_ulps()`` and ``get_similar_ulps()``
   (`c998445 <https://github.com/dhermes/bezier/commit/c998445026a5487c59af17c9cbdfc9a6cf4d72c0>`__).

Non-Public API
~~~~~~~~~~~~~~

-  Returning ``coincident`` flag from curve-curve ``all_intersections``
   (`ebe6617 <https://github.com/dhermes/bezier/commit/ebe66178d0ab6f359ba206ded7b5d629d849955c>`__).
-  Adding a ``TANGENT_BOTH`` classification for surface-surface intersection
   points that are interior to both surfaces at the point of tangency
   (`b89b2b1 <https://github.com/dhermes/bezier/commit/b89b2b1de1726cdc9f508bd761f4c20e7d655321>`__).
   This previously failed with a :exc:`NotImplementedError`.
-  Added ``COINCIDENT`` classification for surface-surface intersection
   points that occur on a segment that is coincident on an edges of each
   surface
   (`8b1c59d <https://github.com/dhermes/bezier/commit/8b1c59d2b48281d38275af6c5b6e11c1699b92c6>`__).
   Such points previously failed classification because they were interpreted
   as being tangent and having the same curvature (because the segments
   are identical).
-  Added a ``COINCIDENT_UNUSED`` classification
   (`cfa2b93 <https://github.com/dhermes/bezier/commit/cfa2b93792695b87f11ece9da1959013ecf77678>`__)
   for cases where coincident segments are moving in opposite directions (i.e.
   the surfaces don't share a common interior). For example see case 44
   (`29Q-43Q <https://github.com/dhermes/bezier/blob/0.8.0/docs/images/surfaces29Q_and_43Q.png>`__).
-  Adding custom linear solver for the ``2 x 2`` case
   (`764e56d <https://github.com/dhermes/bezier/commit/764e56db5bb4987d31e3c9f5fbabbe6564d6e0c0>`__).
   This is modelled after ``dgesv`` from LAPACK.
-  Adding some support for B |eacute| zier clipping algorithm
   (`fbed62d <https://github.com/dhermes/bezier/commit/fbed62df305b8c2679ff260bba4f57d414e79a77>`__,
   `ada4ea3 <https://github.com/dhermes/bezier/commit/ada4ea34bf31cff5cc34491d6689f0f3a2b9f0a1>`__).
   See the original `paper <https://dx.doi.org/10.1016/0010-4485(90)90039-F>`__
   by Sederberg and Nishita for more information.

.. |pypi| image:: https://img.shields.io/pypi/v/bezier/0.8.0.svg
   :target: https://pypi.org/project/bezier/0.8.0/
   :alt: PyPI link to release 0.8.0
.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=0.8.0
   :target: https://bezier.readthedocs.io/en/0.8.0/
   :alt: Documentation for release 0.8.0
