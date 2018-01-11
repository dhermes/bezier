# Changelog

## 0.6.0

PyPI: https://pypi.org/project/bezier/0.6.0/ <br>
Docs: https://bezier.readthedocs.io/en/0.6.0/

### Performance Optimizations

- Added recommended performance flags for `gfortran` based on
  [recommendations][0.6.0-3] on `fortran90.org` ([`3877982`][0.6.0-22]).
  - Extensions can be compiled in debug mode by setting
    `DEBUG=True` ([`b62460b`][0.6.0-62]).
  - Setting `BEZIER_NO_EXTENSIONS=True` will build pure-Python modules
    only ([`3f6280c`][0.6.0-25])
- Added [QUADPACK][0.6.0-86] to use in
  `curve.f90::compute_length` ([`985a4c0`][0.6.0-53]).
- Implemented curve-curve intersection completely in Fortran (e.g.
  [`4a8f801`][0.6.0-28]) which resulted in a 10x speedup when called from
  Python. Also implemented surface-surface intersection completely
  in Fortran, resulting in a 3x speedup.

### Python Changes

#### New Features

- Added `CurvedPolygon._metadata` to track where edges originated, e.g.
  from a surface-surface intersection ([`871d23d`][0.6.0-45]). This is used
  for sanity checking in functional tests ([`e253da2`][0.6.0-78]).
- Made speedup checks specific to the module, not all four. I.e.
  `bezier._HAS_SPEEDUP` was dropped in favor of five members, e.g.
  `_HAS_CURVE_SPEEDUP` ([`d798f66`][0.6.0-73]).
- Added `bezier.__author__` and [`bezier.__version__`][0.6.0-87] attributes.
- Added [`bezier.get_dll()`][0.6.0-88] for Windows ([`699e39b`][0.6.0-34]).
- Added `bezier/__config__.py` that adds `libbezier` to `%PATH%` on
  Windows ([`8538af4`][0.6.0-43]).
- Fortran / Cython speedups added:
  - `_curve_speedup.pyx::subdivide_nodes`
  - `_curve_speedup.pyx::newton_refine`
  - `_curve_speedup.pyx::locate_point`
  - `_curve_speedup.pyx::elevate_nodes`
  - `_curve_speedup.pyx::get_curvature`
  - `_curve_speedup.pyx::reduce_pseudo_inverse`
  - `_curve_speedup.pyx::full_reduce`
  - `_curve_speedup.pyx::compute_length`
  - `_curve_intersection_speedup.pyx::all_intersections`
  - `_curve_intersection_speedup.pyx::free_curve_intersections_workspace`
  - `_helpers_speedup.pyx::contains_nd`
  - `_helpers_speedup.pyx::vector_close`
  - `_helpers_speedup.pyx::in_interval`
  - `_helpers_speedup.pyx::ulps_away`
  - `_surface_speedup.pyx::specialize_surface`
  - `_surface_speedup.pyx::subdivide_nodes`
  - `_surface_speedup.pyx::compute_edge_nodes`
  - `_surface_intersection_speedup.pyx::newton_refine`
  - `_surface_intersection_speedup.pyx::locate_point`
  - `_surface_intersection_speedup.pyx::surface_intersections`
  - `_surface_intersection_speedup.pyx::free_surface_intersections_workspace`

#### Breaking Changes

- [`Curve.intersect()`][0.6.0-89] returns `s-t` parameters rather than
  `x-y` values ([`c309998`][0.6.0-68]).
- [`Surface.intersect()`][0.6.0-90] returns a list with a single `Surface`
  when one of the two surfaces is contained in the
  other ([`05b1fd9`][0.6.0-6]).
- [`Surface.is_valid`][0.6.0-91] will only return `True` if the map `B(s, t)`
  determined by the surface has everywhere positive Jacobian. Previously a
  negative Jacobian was also allowed ([`260fb51`][0.6.0-1]).
- Removed data members from `Curve`:
  - `edge_index` ([`b969488`][0.6.0-63])
  - `next_edge` ([`28619e8`][0.6.0-16])
  - `previous_edge` ([`28619e8`][0.6.0-16])
  - `root` ([`db427f9`][0.6.0-75])
  - `start` ([`39ee98b`][0.6.0-23])
  - `end` ([`39ee98b`][0.6.0-23])
- Removed data members from `Surface`:
  - `base_x` ([`dea75e3`][0.6.0-2])
  - `base_y` ([`dea75e3`][0.6.0-2])
  - `width` ([`dea75e3`][0.6.0-2])
- Remove `dimension` argument in `_curve_speedup.pyx::elevate_nodes` since
  it can be inferred from `nodes` ([`06501c5`][0.6.0-7]).

### ABI Changes

#### New Features

- Fully implemented curve-curve intersection (as
  `curve_intersection.h::curve_intersections`) and surface-surface
  intersection (as `surface_intersection.h::surface_intersections`) at
  the ABI level.
- Added the `surface_intersection.h` header file and implementations for the
  described functions ([`fafd9ff`][0.6.0-84]).
- Newly added functions
  - `curve.h::subdivide_nodes_curve` ([`efb3ce6`][0.6.0-82])
  - `curve.h::newton_refine_curve` ([`2257344`][0.6.0-15])
  - `curve.h::locate_point_curve` ([`2121101`][0.6.0-14],
    [`32b0fa9`][0.6.0-20])
  - `curve.h::elevate_nodes_curve` ([`b03fc28`][0.6.0-60])
  - `curve.h::get_curvature` ([`69cb2f8`][0.6.0-35])
  - `curve.h::reduce_pseudo_inverse` ([`7c3db17`][0.6.0-39])
  - `curve.h::full_reduce` ([`4abd309`][0.6.0-29])
  - `curve.h::compute_length` ([`985a4c0`][0.6.0-53], [`7e71b20`][0.6.0-40])
  - `curve_intersection.h::curve_intersections` ([`c92f98d`][0.6.0-96])
  - `curve_intersection.h::free_curve_intersections_workspace` ([`c92f98d`][0.6.0-96])
  - `helpers.h::contains_nd` ([`36f4b5e`][0.6.0-21])
  - `helpers.h::vector_close` ([`9f3716a`][0.6.0-55])
  - `helpers.h::in_interval` ([`3c0af5d`][0.6.0-24])
  - `helpers.h::ulps_away` ([`0197237`][0.6.0-4])
  - `surface.h::specialize_surface` ([`eb8693e`][0.6.0-81],
    [`fcd5bad`][0.6.0-85])
  - `surface.h::subdivide_nodes_surface` ([`6027210`][0.6.0-32],
    [`4fc5f2a`][0.6.0-30], [`8beb1ac`][0.6.0-47], [`0b2b1f3`][0.6.0-8],
    [`d27b86f`][0.6.0-70], [`88c302b`][0.6.0-46])
  - `surface.h::compute_edge_nodes` ([`2d02590`][0.6.0-17],
    [`f86649a`][0.6.0-83])
  - `surface_intersection.h::newton_refine_surface` ([`93c288d`][0.6.0-50])
  - `surface_intersection.h::locate_point_surface` ([`325ea47`][0.6.0-19],
    [`ca134e6`][0.6.0-69], [`bf69852`][0.6.0-65])
  - `surface_intersection.h::surface_intersections`
  - `surface_intersection.h::free_surface_intersections_workspace`
- Added [`status.h`][0.6.0-97] with an enum for failure states. Each Fortran
  procedure that returns a status documents the possible values and if each
  value is set directly or by a called procedure ([`9fc8575`][0.6.0-56],
  [`c2accf7`][0.6.0-67]).

#### Breaking Changes

- Removed functions
  - `curve.h::specialize_curve_generic` ([`d52453b`][0.6.0-71])
  - `curve.h::specialize_curve_quadratic` ([`d52453b`][0.6.0-71])
  - `curve_intersection.h::from_linearized` ([`d62e462`][0.6.0-72])
  - `curve_intersection.h::bbox_line_intersect` ([`72c0179`][0.6.0-37])
  - `curve_intersection.h::linearization_error` ([`4a3378b`][0.6.0-27])
  - `curve_intersection.h::segment_intersection` ([`4060590`][0.6.0-26])
  - `curve_intersection.h::parallel_different` ([`df3e195`][0.6.0-76])
- Renamed functions
  - `curve.h::newton_refine` to `newton_refine_curve` ([`194ce95`][0.6.0-11])
  - `curve.h::elevate_nodes` to `elevate_nodes_curve` ([`194ce95`][0.6.0-11])
  - `curve_intersection.h::newton_refine_intersect` to
    `newton_refine_curve_intersect` ([`a055525`][0.6.0-57])
- Replaced `degree` with `num_nodes (== degree + 1)` in functions that
  operate on curves:
  - `curve.h::evaluate_curve_barycentric` ([`13eacdd`][0.6.0-10])
  - `curve.h::evaluate_multi` ([`962c288`][0.6.0-52])
  - `curve.h::specialize_curve` ([`ac86233`][0.6.0-59])
  - `curve.h::evaluate_hodograph` ([`9170855`][0.6.0-49])
  - `curve_intersection.h::newton_refine_curve_intersect` ([`80ec491`][0.6.0-42])

### Miscellany

- Added documentation for "native extensions" in
  `DEVELOPMENT` ([`2f9f2c4`][0.6.0-92]).
- Overhauled [`native-libraries` doc][0.6.0-95] with subsections for OS X and
  Windows ([`bfa75ee`][0.6.0-66], [`72005fb`][0.6.0-94], etc.).
- Added Fortran unit tests ([`758bdd1`][0.6.0-38], [`e8afba7`][0.6.0-79],
  [`3164365`][0.6.0-18], etc.).
- Began testing in Mac OS X on Travis ([`9ac5e8e`][0.6.0-54],
  [`85f7619`][0.6.0-44], etc.).
- Added a workaround (`include/bezier/_bool_patch.h`) for the missing support
  for `bool` in old MSVC versions that are required to work with Python
  2.7 ([`5577178`][0.6.0-93]).

[0.6.0-1]: https://github.com/dhermes/bezier/commit/260fb512a67376a7b62b41c37377306c743c8b61
[0.6.0-2]: https://github.com/dhermes/bezier/commit/dea75e3f2999e52f74c3d2603a4e162ae3eb2ef2
[0.6.0-3]: http://www.fortran90.org/src/faq.html
[0.6.0-4]: https://github.com/dhermes/bezier/commit/01972377303afa9a41c79533ee967b8bcc526435
[0.6.0-5]: https://github.com/dhermes/bezier/commit/02d9e115f0c0a7b89bb7b79e3308f1ac00448d71
[0.6.0-6]: https://github.com/dhermes/bezier/commit/05b1fd98c5caea015b87819abdd2d6631ccc9bd4
[0.6.0-7]: https://github.com/dhermes/bezier/commit/06501c5c05e4f646e756f225dc2db0fab98cbbab
[0.6.0-8]: https://github.com/dhermes/bezier/commit/0b2b1f3edeab418bbacc8cd60a419a201cfdd038
[0.6.0-9]: https://github.com/dhermes/bezier/commit/0c9e1be84ce413426066e1abc556a118df759d59
[0.6.0-10]: https://github.com/dhermes/bezier/commit/13eacdd189b81acbbcc3e39bb6643c6edf6a4750
[0.6.0-11]: https://github.com/dhermes/bezier/commit/194ce95a8721e014d3a7d73213d358d89bc81fd8
[0.6.0-12]: https://github.com/dhermes/bezier/commit/1af51cf21866734c3076430d667ac770a3a561b5
[0.6.0-13]: https://github.com/dhermes/bezier/commit/1f55630119980053218d6dcc20cea356e6424705
[0.6.0-14]: https://github.com/dhermes/bezier/commit/21211010bdf5640b6bbfbbd6a270d35dc928efbc
[0.6.0-15]: https://github.com/dhermes/bezier/commit/2257344abc501d4456f9b819969cf8dd9cbefb0b
[0.6.0-16]: https://github.com/dhermes/bezier/commit/28619e8ab4f7b205c7676e45e6339468a5d92460
[0.6.0-17]: https://github.com/dhermes/bezier/commit/2d02590ed972dba902958a07598b95f8099a7295
[0.6.0-18]: https://github.com/dhermes/bezier/commit/3164365b261564c0da158ab46d899357735fbd31
[0.6.0-19]: https://github.com/dhermes/bezier/commit/325ea479665947016844b3ea37cbccf5962f5876
[0.6.0-20]: https://github.com/dhermes/bezier/commit/32b0fa953f7f4dd9ea8f0b3439206f6185f3d863
[0.6.0-21]: https://github.com/dhermes/bezier/commit/36f4b5e9a6f872178781390b352e42ad15a4d9e1
[0.6.0-22]: https://github.com/dhermes/bezier/commit/387798248cf452b27b8e7aa16b83417b1cdcb196
[0.6.0-23]: https://github.com/dhermes/bezier/commit/39ee98b5add3333b90f73279b304dbb1fd0f2c54
[0.6.0-24]: https://github.com/dhermes/bezier/commit/3c0af5d32aa494efc15cb23a3bf1d15c0b5859b1
[0.6.0-25]: https://github.com/dhermes/bezier/commit/3f6280ccfa1b7dbb9415aaf088dcc610ac4bd8ac
[0.6.0-26]: https://github.com/dhermes/bezier/commit/40605901872d679956b384bedf426f6cbf7a43c5
[0.6.0-27]: https://github.com/dhermes/bezier/commit/4a3378b11f54d582f7f223c92939745fe8daaa4c
[0.6.0-28]: https://github.com/dhermes/bezier/commit/4a8f80177b56f1241e96e973639cfc8a1273e080
[0.6.0-29]: https://github.com/dhermes/bezier/commit/4abd309ff0125bf82f91a71946511e82fd7eaf8a
[0.6.0-30]: https://github.com/dhermes/bezier/commit/4fc5f2a552a91ef738f062ff1578d2865672c9f6
[0.6.0-31]: https://github.com/dhermes/bezier/commit/515909531ea864b54f6c6a2b50409bb003e6ccd3
[0.6.0-32]: https://github.com/dhermes/bezier/commit/602721004887cd17e977e7255c7f574cd321d032
[0.6.0-33]: https://github.com/dhermes/bezier/commit/617d9ca6d63d6edeaf3a01d2b7e97135b027f0ec
[0.6.0-34]: https://github.com/dhermes/bezier/commit/699e39b2671daed23a57d37b5e776c6627c72850
[0.6.0-35]: https://github.com/dhermes/bezier/commit/69cb2f852d82076c3e7c98ce68097bc3b8b4a5b6
[0.6.0-36]: https://github.com/dhermes/bezier/commit/71d27d16eda291b4e831ef0caa7cc793b636a8dd
[0.6.0-37]: https://github.com/dhermes/bezier/commit/72c017995f7df9bfa42347dbc0f967e666bbadee
[0.6.0-38]: https://github.com/dhermes/bezier/commit/758bdd15426424c1566bda15c03594bc2e66410a
[0.6.0-39]: https://github.com/dhermes/bezier/commit/7c3db1727c45763ec7a14550764979cb9ceafcb5
[0.6.0-40]: https://github.com/dhermes/bezier/commit/7e71b202cfcb2ea575771c2e1169d0e0c27e481b
[0.6.0-41]: https://github.com/dhermes/bezier/commit/7f6fafb4248ad38ca8b662da1a7ba877a74743ab
[0.6.0-42]: https://github.com/dhermes/bezier/commit/80ec491d5d7094d378de5110b2b254e33ca271a1
[0.6.0-43]: https://github.com/dhermes/bezier/commit/8538af47870849e6f55e0c861e0e2b720aa3ae75
[0.6.0-44]: https://github.com/dhermes/bezier/commit/85f7619929debd3730d6ddafa4ac75789ad8e5f3
[0.6.0-45]: https://github.com/dhermes/bezier/commit/871d23d640e9333f5e7ae02e29ca878c80176682
[0.6.0-46]: https://github.com/dhermes/bezier/commit/88c302b5d070d512339bd1b16960a4e43025005c
[0.6.0-47]: https://github.com/dhermes/bezier/commit/8beb1ace5934d7bce03cb19c483aa1aec57ec06b
[0.6.0-48]: https://github.com/dhermes/bezier/commit/911aeb1ba1ca5b828e53297c277d1c68bbd6755f
[0.6.0-49]: https://github.com/dhermes/bezier/commit/91708552355b71e95adb8454ec69d1f3d1e81c22
[0.6.0-50]: https://github.com/dhermes/bezier/commit/93c288d5e5865986aa2627ea81f12b6370099865
[0.6.0-51]: https://github.com/dhermes/bezier/commit/9578a6c76f6d5382189ba33285026cff16e733b4
[0.6.0-52]: https://github.com/dhermes/bezier/commit/962c288a2eb3c0c7fdeb6f055ebbff57331b7cf5
[0.6.0-53]: https://github.com/dhermes/bezier/commit/985a4c07c25c9bb362c6de615580a4246b3c16b1
[0.6.0-54]: https://github.com/dhermes/bezier/commit/9ac5e8e4a02ce8b64b9e6b5142a6c1fda01ee787
[0.6.0-55]: https://github.com/dhermes/bezier/commit/9f3716ae18ea5e592db01f37242623b718e23a84
[0.6.0-56]: https://github.com/dhermes/bezier/commit/9fc857553fd1b525e801970cfe0eb2e2288f5319
[0.6.0-57]: https://github.com/dhermes/bezier/commit/a055525c1ab81246bc6d040fdce376772cf65703
[0.6.0-58]: https://github.com/dhermes/bezier/commit/aa0f205b552e3436dea104b5068692dc7d68d75f
[0.6.0-59]: https://github.com/dhermes/bezier/commit/ac86233dea35c56f5b5c81fd1020ca480487d87c
[0.6.0-60]: https://github.com/dhermes/bezier/commit/b03fc280053043996138932043de4f6ac69e16ce
[0.6.0-61]: https://github.com/dhermes/bezier/commit/b06bde03025b41986e35039fd312b82d824aff85
[0.6.0-62]: https://github.com/dhermes/bezier/commit/b62460b47faec2666fceb13457bc11558f8079e9
[0.6.0-63]: https://github.com/dhermes/bezier/commit/b969488df308ee78712db0b605a86d3adc0d3da6
[0.6.0-64]: https://github.com/dhermes/bezier/commit/ba99c8d89e060f3b481d484aded306e8a73ca163
[0.6.0-65]: https://github.com/dhermes/bezier/commit/bf698525d99f2424717ba4a8559ea2ab84abe6cb
[0.6.0-66]: https://github.com/dhermes/bezier/commit/bfa75eedda8187bef59ca8e04f9d04ee0fc28b97
[0.6.0-67]: https://github.com/dhermes/bezier/commit/c2accf76047741d7a42327a67c4732b488c56600
[0.6.0-68]: https://github.com/dhermes/bezier/commit/c309998f705de7467e1222e41467190739ff3118
[0.6.0-69]: https://github.com/dhermes/bezier/commit/ca134e63f1061b99404e3dfedefdd5d8cf5956ea
[0.6.0-70]: https://github.com/dhermes/bezier/commit/d27b86f417c4f6227c08bfc0c50da88584d7996b
[0.6.0-71]: https://github.com/dhermes/bezier/commit/d52453ba27c993422da3fbbc78c53aea960fd525
[0.6.0-72]: https://github.com/dhermes/bezier/commit/d62e462507287af66f51043889ae56be21cb8e45
[0.6.0-73]: https://github.com/dhermes/bezier/commit/d798f665fdd0c223a1b0a71919eb9a45ee86951d
[0.6.0-74]: https://github.com/dhermes/bezier/commit/d92357913bdbdfccf43ef7f393f17fde99830b38
[0.6.0-75]: https://github.com/dhermes/bezier/commit/db427f958266103e0a266721e18c66a7025d85ae
[0.6.0-76]: https://github.com/dhermes/bezier/commit/df3e195a77236799ed975d3b9251c45eb4bbf29a
[0.6.0-77]: https://github.com/dhermes/bezier/commit/e1777e767a87ba54689fbd3ae86bd98568212ff8
[0.6.0-78]: https://github.com/dhermes/bezier/commit/e253da2e507b4edf1627ddd753a6e5e73e563684
[0.6.0-79]: https://github.com/dhermes/bezier/commit/e8afba7c64a7ecc2ded84efb5f164513974963cf
[0.6.0-80]: https://github.com/dhermes/bezier/commit/eaec7cfbba84cfb27247647be52286551df9e619
[0.6.0-81]: https://github.com/dhermes/bezier/commit/eb8693e823f7a2af3fb6682f76a667b20d419e5d
[0.6.0-82]: https://github.com/dhermes/bezier/commit/efb3ce65cef671d1745028594bb5f0897e96e053
[0.6.0-83]: https://github.com/dhermes/bezier/commit/f86649aa1f631a11ad314754a236719dd6f0c714
[0.6.0-84]: https://github.com/dhermes/bezier/commit/fafd9ff181755e9b204372f1b94dd10578a16382
[0.6.0-85]: https://github.com/dhermes/bezier/commit/fcd5bad6fe47499c23db32092e7a749e2e866f92
[0.6.0-86]: https://en.wikipedia.org/wiki/QUADPACK
[0.6.0-87]: http://bezier.readthedocs.io/en/0.6.0/reference/bezier.html#bezier.__version__
[0.6.0-88]: http://bezier.readthedocs.io/en/0.6.0/reference/bezier.html#bezier.get_dll
[0.6.0-89]: http://bezier.readthedocs.io/en/0.6.0/reference/bezier.curve.html#bezier.curve.Curve.intersect
[0.6.0-90]: http://bezier.readthedocs.io/en/0.6.0/reference/bezier.surface.html#bezier.surface.Surface.intersect
[0.6.0-91]: http://bezier.readthedocs.io/en/0.6.0/reference/bezier.surface.html#bezier.surface.Surface.is_valid
[0.6.0-92]: https://github.com/dhermes/bezier/commit/2f9f2c49585238024722b5a5b4fb60ea3338b9b3
[0.6.0-93]: https://github.com/dhermes/bezier/commit/5577178a3ef45487667cd72d81146390be2b0c41
[0.6.0-94]: https://github.com/dhermes/bezier/commit/72005fbb0a05715f6832f68dc8c3f04576781047
[0.6.0-95]: http://bezier.readthedocs.io/en/0.6.0/native-libraries.html
[0.6.0-96]: https://github.com/dhermes/bezier/commit/c92f98dde15a2ad3cbdf8db46fd48fbbed105552
[0.6.0-97]: https://github.com/dhermes/bezier/blob/0.6.0/src/bezier/include/bezier/status.h

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
