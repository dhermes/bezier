# Changelog

## 0.5.0

PyPI: https://pypi.org/project/bezier/0.5.0/ <br>
Docs: https://bezier.readthedocs.io/en/0.5.0/

### Performance Optimizations

-   Change `wiggle_interval` to return `success` bool instead of raising an
    exception. This allows the implicitization approach to use it without
    having to use exceptions for flow-control. (Fixes [#22][0.5.0-5].)
-   Switching Fortran speedups from `f2py` to Cython (this is because `f2py`
    artificially limits the feature set of Fortran, i.e. user defined types)
-   Moving some more code to Fortran (e.g. `bbox_line_intersect()`
    [`3dcf640`][0.5.0-11])

### New Features

-   Making Fortran features available outside of Python (see
    [Native Libraries][0.5.0-1])
-   C headers for each Fortran module (via
    [`bezier.get_include()`][0.5.0-2])
    -   Cython `.pxd` declarations for all Fortran modules
    -   `libbezier` static library (via [`bezier.get_lib()`][0.5.0-3])
-   Implementing [`bezier_roots()`][0.5.0-13] polynomial root solver for
    polynomials written in Bernstein basis. ([`0dd6369`][0.5.0-12])

### Miscellany

-   Getting `bezier` [published][0.5.0-10] in the Journal of Open Source
    Science (JOSS). See [review][0.5.0-9]. ([`e6c4536`][0.5.0-7] and
    [`975ac6b`][0.5.0-8])
-   Updating error message for `locate()` methods and adding a note that
    `locate()` / `evaluate*()` are (essentially) inverses. H/T to
    [@pdknsk][0.5.0-22] [#36][0.5.0-4]
-   Using Fortran-contiguous arrays in
    `_check_non_simple()`. ([`b06c78e`][0.5.0-6])
-   Moving most of `Curve.subdivide()` and `Surface.subdivide()` logic into
    helpers. This is part of an effort to make all helpers take low-level
    data types rather than `Curve`s, `Surface`s, etc. ([`34515bd`][0.5.0-14]
    and [`1fc80e5`][0.5.0-15])
-   Split `speedup.f90` into submodules `curve.f90`, `surface.f90`,
    etc. ([`75349b7`][0.5.0-16], [`dfd6bba`][0.5.0-17], [`7096a9d`][0.5.0-18],
    [`c326c00`][0.5.0-19])
-   Adding `BEZIER_JOURNAL` option to `setup.py`. This stores a record of
    compiler commands invoked during installation. See
    [Native Libraries][0.5.0-1] for more details. ([`3d832e7`][0.5.0-20] and
    [`c64a97a`][0.5.0-21])

[0.5.0-1]: http://bezier.readthedocs.io/en/0.5.0/native-libraries.html
[0.5.0-2]: http://bezier.readthedocs.io/en/0.5.0/reference/bezier.html#bezier.get_include
[0.5.0-3]: http://bezier.readthedocs.io/en/0.5.0/reference/bezier.html#bezier.get_lib
[0.5.0-4]: https://github.com/dhermes/bezier/pull/36
[0.5.0-5]: https://github.com/dhermes/bezier/pull/22
[0.5.0-6]: https://github.com/dhermes/bezier/commit/b06c78e50d53bf673bcf0b71fa84b36c8df564d8
[0.5.0-7]: https://github.com/dhermes/bezier/commit/e6c45360f0c8412ae90d967463a14c49490d70ee
[0.5.0-8]: https://github.com/dhermes/bezier/commit/975ac6b1a4313db4dcdc17396d6d34561005939e
[0.5.0-9]: https://github.com/openjournals/joss-reviews/issues/267
[0.5.0-10]: http://joss.theoj.org/papers/10.21105/joss.00267
[0.5.0-11]: https://github.com/dhermes/bezier/commit/3dcf64090bb5874320dcde86eaf449e94278dd08
[0.5.0-12]: https://github.com/dhermes/bezier/commit/0dd6369b0f77e4c0cf8113f2d25812addc90482a
[0.5.0-13]: http://bezier.readthedocs.io/en/0.5.0/algorithm-helpers.html#bezier._implicitization.bezier_roots
[0.5.0-14]: https://github.com/dhermes/bezier/commit/34515bd6246f57fbb311b4089520a24e8237294a
[0.5.0-15]: https://github.com/dhermes/bezier/commit/1fc80e54ad1b45cb628af06e5a2100eeb9282865
[0.5.0-16]: https://github.com/dhermes/bezier/commit/75349b745063a9bbc623808b3f7bbf6b7641c008
[0.5.0-17]: https://github.com/dhermes/bezier/commit/dfd6bba303ac0a8492fac1f309086b685e52ab59
[0.5.0-18]: https://github.com/dhermes/bezier/commit/7096a9d646930378476e650c77d0652a48bf148a
[0.5.0-19]: https://github.com/dhermes/bezier/commit/c326c00a5c0ee74f9aa53c2b104ac6d4eb5c6794
[0.5.0-20]: https://github.com/dhermes/bezier/commit/3d832e78af2a951a642ff5860b9593abfa674ec3
[0.5.0-21]: https://github.com/dhermes/bezier/commit/c64a97aa5599220b927094a41de04b0c75bbec33
[0.5.0-22]: https://github.com/pdknsk

## 0.4.0

PyPI: https://pypi.org/project/bezier/0.4.0/ <br>
Docs: https://bezier.readthedocs.io/en/0.4.0/

### Performance Optimizations

-   [Adding][0.4.0-29] Fortran [speedups][0.4.0-30] for many crucial computation helpers
    including
    -   intersecting line segments
    -   (vectorized) Horner's method for evaluating a B&#xe9;zier curve at
        multiple parameters at once
    -   (vectorized) Horner's method for evaluating a B&#xe9;zier surface
    -   computing "linearization error" (how close a curve is to a line)
    -   specializing a B&#xe9;zier curve to a sub-interval
    -   using Newton's method to refine a curve-curve intersection
-   [Adding][0.4.0-10] `_verify` switch to [`Surface.locate()`][0.4.0-11] and
    [`Curve.intersect()`][0.4.0-14] to selectively disable overly defensive
    value checking. (Making sure to use this switch during "internal"
    computation.)
-   Making sure NumPy arrays are Fortan-contiguous as often as possible (e.g.
    snippets and source, via `np.asfortranarray()`). This is to avoid (and
    emphasize) a non-trivial overhead when passing a C-contiguous array to a
    Fortran function. ([`03a7242`][0.4.0-15], [`6064e4c`][0.4.0-16],
    [`f1804f4`][0.4.0-17])
-   Using Horner's method in `Curve.evaluate_multi()` and
    `Surface.evaluate_barycentric()`, rather than inferior (sometimes
    non-vectorized) approaches ([`dee8181`][0.4.0-18], [`2611e64`][0.4.0-19])
-   Made surface-surface intersection more resilient / lenient for corner
    intersections. For "nearby" intersections, parameter values can be
    rounded to `0` or `1`. ([`4a8458c`][0.4.0-25])

### New Features

-   [Adding][0.4.0-23] optional `strategy` argument (one of geometric or
    algebraic) to [`Surface.intersect()`][0.4.0-24]
    -   Added "algebraic" [`IntersectionStrategy`][0.4.0-20] via curve
        [implicitization][0.4.0-27] ([reference][0.4.0-28])
-   Adding [`Curve.reduce_()`][0.4.0-21] which acts as a partial inverse to
    [`Curve.elevate()`][0.4.0-31]. It is only a complete inverse when a curve
    is degree-elevated, otherwise it returns the "best" reduced form (in the
    least squares sense).

### Interface Changes

-   (**Breaking change**) [Removing][0.4.0-5] `show` keyword from
    [`Curve.plot()`][0.4.0-2], [`Surface.plot()`][0.4.0-3] and
    [`CurvedPolygon.plot()`][0.4.0-4]
-   [Adding][0.4.0-32] `color` keyword to [`Curve.plot()`][0.4.0-2]
-   [Adding][0.4.0-26] `alpha` keyword to [`Curve.plot()`][0.4.0-2]
-   (**Breaking change**) [Splitting][0.4.0-6] the
    [`Surface.evaluate_multi()`][0.4.0-7] method into
    [`Surface.evaluate_barycentric_multi()`][0.4.0-8] and
    [`Surface.evaluate_cartesian_multi()`][0.4.0-9]
-   [Adding][0.4.0-22] `__dict__` helpers on `Curve`, `CurvedPolygon` and
    `Surface`. These are `@property`s intended only for REPL use, since
    classes with `__slots__` no longer have a `__dict__` attribute.

### Miscellany

-   Adding [`IntersectionClassification`][0.4.0-1] to docs ([ref][0.4.0-5])
-   [Moving][0.4.0-12] most plotting into a dedicated module. More
    importantly, importing plotting helpers at **run-time** rather at
    **import time**. So if computational code never plots, it won't eat the
    import cost of `matplotlib`. [Removing][0.4.0-13] `matplotlib` as a
    dependency.

[0.4.0-1]: http://bezier.readthedocs.io/en/0.4.0/algorithm-helpers.html#bezier._surface_helpers.IntersectionClassification
[0.4.0-2]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.curve.html#bezier.curve.Curve.plot
[0.4.0-3]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.surface.html#bezier.surface.Surface.plot
[0.4.0-4]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.curved_polygon.html#bezier.curved_polygon.CurvedPolygon.plot
[0.4.0-5]: https://github.com/dhermes/bezier/commit/828f4238971b12a9d494ce38387cec855d063c91
[0.4.0-6]: https://github.com/dhermes/bezier/commit/cea88285b8c9002a57efd88e69b5bd2ef46e7ca7
[0.4.0-7]: http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.evaluate_multi
[0.4.0-8]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.surface.html#bezier.surface.Surface.evaluate_barycentric_multi
[0.4.0-9]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.surface.html#bezier.surface.Surface.evaluate_cartesian_multi
[0.4.0-10]: https://github.com/dhermes/bezier/commit/dcf40f4c9ed2167e96fc8f4675aeedcc2d811a0b
[0.4.0-11]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.surface.html#bezier.surface.Surface.locate
[0.4.0-12]: https://github.com/dhermes/bezier/commit/dc4d33cfcf7f9ac6e794b856dc6d76635d362922
[0.4.0-13]: https://github.com/dhermes/bezier/commit/064e2c5efe7fa6498d74a33798a363e2c8e0b83e
[0.4.0-14]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.curve.html#bezier.curve.Curve.intersect
[0.4.0-15]: https://github.com/dhermes/bezier/commit/03a72428c6f9d3bd3a1fac9b7f9afa615ce12d46
[0.4.0-16]: https://github.com/dhermes/bezier/commit/6064e4c314d8d717873d46e6ef35c0bbc9772728
[0.4.0-17]: https://github.com/dhermes/bezier/commit/f1804f442f190d0bc36782e940ee0b8a68c5ecd6
[0.4.0-18]: https://github.com/dhermes/bezier/commit/dee81813e34d5f69c52f48aa90f7c11eb4ddc3ec
[0.4.0-19]: https://github.com/dhermes/bezier/commit/2611e64a735e46317cce08a41270d61024705fd9
[0.4.0-20]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.curve.html#bezier.curve.IntersectionStrategy
[0.4.0-21]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.curve.html#bezier.curve.Curve.reduce_
[0.4.0-22]: https://github.com/dhermes/bezier/commit/f0fca088ac6f70c39f9f5af457c29e3c82f094b5
[0.4.0-23]: https://github.com/dhermes/bezier/commit/e72ca20f0f4ee0f6399b56805b30fe67a02aa04f
[0.4.0-24]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.surface.html#bezier.surface.Surface.intersect
[0.4.0-25]: https://github.com/dhermes/bezier/commit/4a8458c823d8acc185818f856889cff6f46300d3
[0.4.0-26]: https://github.com/dhermes/bezier/commit/dcbeefc25b7f5f9a1fa725dac04e81a43039f680
[0.4.0-27]: https://github.com/dhermes/bezier/commits/0.4.0/src/bezier/_implicitization.py
[0.4.0-28]: https://en.wikipedia.org/wiki/Resultant
[0.4.0-29]: https://github.com/dhermes/bezier/commits/0.4.0/src/bezier/speedup.f90
[0.4.0-30]: https://github.com/dhermes/bezier/blob/0.4.0/src/bezier/speedup.f90
[0.4.0-31]: http://bezier.readthedocs.io/en/0.4.0/reference/bezier.curve.html#bezier.curve.Curve.elevate
[0.4.0-32]: https://github.com/dhermes/bezier/commit/ce838a2aaef2281f06603d1c76324a3aa8289cf9

## 0.3.0

PyPI: https://pypi.org/project/bezier/0.3.0/ <br>
Docs: http://bezier.readthedocs.io/en/0.3.0/

### Performance Optimizations

-   Adding `__slots__` for all classes
-   Removing all usage of `@property` calls from internal callers (to avoid
    function call overhead)
-   Avoiding un-necessary data copying, e.g. `nodes[[0], :]` creates a copy but

    ``` python
    nodes[0, :].reshape((1, 2))
    ```

    does not ([more details](https://docs.scipy.org/doc/numpy-1.6.0/reference/arrays.indexing.html#advanced-indexing))
-   Adding `_verify` switches to selectively disable overly defensive value
    checking. Added to
    [`CurvedPolygon`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.curved_polygon.html#bezier.curved_polygon.CurvedPolygon)
    constructor,
    [`Surface.evaluate_barycentric()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.evaluate_barycentric),
    [`Surface.evaluate_cartesian()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.evaluate_cartesian),
    [`Surface.evaluate_multi()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.evaluate_multi) and
    [`Surface.intersect()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.intersect).
    Internal callers with already verified data now skip verification steps
-   [Bailing out early](https://github.com/dhermes/bezier/commit/db816eb5a748bb997adcc2d7d9008638e22a824c)
    if surface bounding boxes are disjoint in
    [`Surface.intersect()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.intersect)

### Breaking Changes

-   Requiring `degree` in
    [`Curve`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.curve.html#bezier.curve.Curve) and
    [`Surface`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface)
    constructors, but adding
    [`Curve.from_nodes()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.curve.html#bezier.curve.Curve.from_nodes) and
    [`Surface.from_nodes()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.from_nodes)
    factories to accept nodes only (computing the degree in the constructor
    every time is a waste of flops, especially if the caller knows the degree)
-   [Removing](https://github.com/dhermes/bezier/commit/3393b9010c26b55a9c29afc2702426bb179b85a1)
    public
    [`Curve.copy()`](http://bezier.readthedocs.io/en/0.2.1/reference/bezier.curve.html#bezier.curve.Curve.copy)
    and
    [`Surface.copy()`](http://bezier.readthedocs.io/en/0.2.1/reference/bezier.surface.html#bezier.surface.Surface.copy)
-   [Removing](https://github.com/dhermes/bezier/commit/3393b9010c26b55a9c29afc2702426bb179b85a1)
    custom equality checks for
    [`Curve`](http://bezier.readthedocs.io/en/0.2.1/reference/bezier.curve.html#bezier.curve.Curve.__eq__)
    and
    [`Surface`](http://bezier.readthedocs.io/en/0.2.1/reference/bezier.surface.html#bezier.surface.Surface.__eq__)
    objects. The previous implementation did not factor in all relevant values
-   Returning `1xD` arrays
    [instead of flattened](https://github.com/dhermes/bezier/commit/b5e5b327594c6143956ed98703f596ff82b7501a)
    `D`-dimensional 1D arrays from
    [`Curve.evaluate()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.curve.html#bezier.curve.Curve.evaluate),
    [`Surface.evaluate_barycentric()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.evaluate_barycentric),
    [`Surface.evaluate_cartesian()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.evaluate_cartesian),
    and related helpers
-   Renaming
    [`Intersection.left/right`](http://bezier.readthedocs.io/en/0.2.1/algorithm-helpers.html#bezier._intersection_helpers.Intersection.left)
    as
    [`first/second`](http://bezier.readthedocs.io/en/0.3.0/algorithm-helpers.html#bezier._intersection_helpers.Intersection.first)
    (They were poorly named originally, since "left" and "right" were in
    reference to where they were used **in code**, not geometry. This class
    is not part of the public interface, but it is documented.)

### Bug Fixes

-   Handling cases where one corner of a surface touches another but their
    interiors don't intersect (in
    [`Surface.intersect()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.intersect)).
    Adding `ignored_corner` classification to handle these curve-curve
    intersecions that don't contribute to a surface-surface intersection
-   Throwing exception in
    [`Curve.locate()`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.curve.html#bezier.curve.Curve.locate)
    when the subdivided intervals are very far apart ([#13][0.3.0-2])
-   Improving
    [`Surface.is_valid`](http://bezier.readthedocs.io/en/0.3.0/reference/bezier.surface.html#bezier.surface.Surface.is_valid)
    by considering the signs of the Jacobian determinant at corner
    nodes ([#12][0.3.0-1])

### Miscellany

-   Adding possible strategy to avoid linear convergence in
    [`newton_refine()`](http://bezier.readthedocs.io/en/0.3.0/algorithm-helpers.html#bezier._intersection_helpers.newton_refine)
-   Adding AppVeyor configuration to make sure there are no Windows
    issues, testing exclusively with `conda` install
-   Updating generated images with `matplotlib` 2.0

[0.3.0-1]: https://github.com/dhermes/bezier/issues/12
[0.3.0-2]: https://github.com/dhermes/bezier/issues/13

## 0.2.1

PyPI: https://pypi.org/project/bezier/0.2.1/ <br>
Docs: http://bezier.readthedocs.io/en/0.2.1/

- Added
  [`Curve.locate()`](http://bezier.readthedocs.io/en/0.2.1/reference/bezier.curve.html#bezier.curve.Curve.locate)
  and
  [`_curve_helpers.newton_refine()`](http://bezier.readthedocs.io/en/0.2.1/algorithm-helpers.html#bezier._curve_helpers.newton_refine)
  helper
- Adding optional `color` to
  [`Surface.plot()`](http://bezier.readthedocs.io/en/0.2.1/reference/bezier.surface.html#bezier.surface.Surface.plot)
- Adding
  [`Surface.elevate()`](http://bezier.readthedocs.io/en/0.2.1/reference/bezier.surface.html#bezier.surface.Surface.elevate)
  for degree elevation
- Fixing nodes defining the
  [self-intersecting curve](http://bezier.readthedocs.io/en/0.2.1/curve-curve-intersection.html#detecting-self-intersections)
  in `curve-curve-intersection` (no code in `bezier` was broken / fixed,
  just "bad" docs)
- Allow wiggle outside of `[0, 1]` when intersecting linearizations in
  `from_linearized()`
- Collapsing almost-same parameter values in `intersect_one_round()` (via
  `from_linearized()`). Previously checked for bitwise equality and relied
  on checking values at the boundary of a subdivided interval
- Adding non-public `bezier._plot_helpers` module

## 0.2.0

PyPI: https://pypi.org/project/bezier/0.2.0/ <br>
Docs: http://bezier.readthedocs.io/en/0.2.0/

- **Primary feature**:
  [`Surface.intersect()`](http://bezier.readthedocs.io/en/0.2.0/reference/bezier.surface.html#bezier.surface.Surface.intersect)
  added
- To support intersection, needed
  [`CurvedPolygon`](http://bezier.readthedocs.io/en/0.2.0/reference/bezier.curved_polygon.html#bezier.curved_polygon.CurvedPolygon),
  i.e. an object defined only by its curved sides (whereas a `Surface`
  may have interior control points)
- [`Curve.specialize`](http://bezier.readthedocs.io/en/0.2.0/reference/bezier.curve.html#bezier.curve.Curve.specialize)
  for chopping a `Curve` at arbitrary parameter values (this is also used
  in surface-surface intersection)
- Added images to most documented functions and methods to illustrate the
  concept at hand. For example
  [`classify_intersection`](http://bezier.readthedocs.io/en/0.2.0/algorithm-helpers.html#bezier._surface_helpers.classify_intersection)
  has **seven** images to enumerate all of the possible cases covered in
  the algorithm.
- Added
  [`Surface.locate()`](http://bezier.readthedocs.io/en/0.2.0/reference/bezier.surface.html#bezier.surface.Surface.locate),
  made possible by
  [`newton_refine`](http://bezier.readthedocs.io/en/0.2.0/algorithm-helpers.html#bezier._surface_helpers.newton_refine)
- Added
  [Algorithm Helpers](http://bezier.readthedocs.io/en/0.2.0/algorithm-helpers.html)
  doc to try to explain some of the core algorithms at work (not all are
  documented yet). Some of this content was previously documented in the
  `bezier.curve` module, but was moved. Since, documentation has been added
  for `get_curvature`, `newton_refine` (for surfaces), `classify_intersection`
  (to determine how two curves interact while intersecting) and for some
  helper classes.
- Added `Surface.base_x`, `Surface.base_y` and `Surface.width`
  [properties](http://bezier.readthedocs.io/en/0.2.0/reference/bezier.surface.html#bezier.surface.Surface.width)
  to allow tracking a sub-surface during the subdivision process (this is an
  analogue to the `Curve.start` and `Curve.end`
  [properties](http://bezier.readthedocs.io/en/0.2.0/reference/bezier.curve.html#bezier.curve.Curve.start))
- Added `Curve.edge_index`, `Curve.next_edge` and `Curve.previous_edge`
  [properties](http://bezier.readthedocs.io/en/0.2.0/reference/bezier.curve.html#bezier.curve.Curve.edge_index)
  to allow tracking when curves are actually the sides of a `Surface`

## 0.1.1

PyPI: https://pypi.org/project/bezier/0.1.1/ <br>
Docs: http://bezier.readthedocs.io/en/0.1.1/

Changes:

- Adding
  [`Curve.elevate()`](http://bezier.readthedocs.io/en/0.1.1/reference/bezier.curve.html#bezier.curve.Curve.elevate) for degree elevation
- Upgrading curve-curve intersection algorithm to ignore parallel line
  segments that don't meet (rather than throwing `NotImplementedError`)
- Making
  [`segment_intersection()`](http://bezier.readthedocs.io/en/0.1.1/reference/bezier.curve.html#bezier._intersection_helpers.segment_intersection)
  return a `success` bool instead of raising `NotImplementedError` on failure
- Updating docs for
  [`newton_refine()`](http://bezier.readthedocs.io/en/0.1.1/reference/bezier.curve.html#bezier._intersection_helpers.newton_refine)
  with two examples and making
  [`parallel_different()`](http://bezier.readthedocs.io/en/0.1.1/reference/bezier.curve.html#bezier._intersection_helpers.parallel_different)
  a publicly documented function (as a partner to `segment_intersection()`)
- Adding some more examples / failures to `curve-curve-intersection`
  [doc](http://bezier.readthedocs.io/en/0.1.1/curve-curve-intersection.html)

## 0.1.0

Second Release: https://pypi.org/project/bezier/0.1.0/ <br>
Docs: http://bezier.readthedocs.io/en/0.1.0/

Primary changes since previous release
([`0.0.1`](https://pypi.org/project/bezier/0.0.1/)) are related to
curve-curve intersection. See
[the intersection docs](http://bezier.readthedocs.io/en/0.1.0/curve-curve-intersection.html)
for examples and more information.

## 0.0.1

PyPI: https://pypi.org/project/bezier/0.0.1/
