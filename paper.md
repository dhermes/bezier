---
title: 'Helper for B&#xe9;zier Curves, Triangles, and Higher Order Objects'
tags:
  - algebra
  - bezier
  - curve
  - parametric curve
  - numerical software
  - Python
  - Fortran
authors:
 - name: Danny Hermes
   orcid: 0000-0001-7366-173X
   affiliation: 1
affiliations:
 - name: University of California, Berkeley
   index: 1
date: 22 May 2017
bibliography: paper.bib

---

# Summary

`bezier` is a Python helper for B&#xe9;zier [curves][1] and
[surfaces][2] [@SederbergNotes, @Farin2001]. B&#xe9;zier curves (and
surfaces) are parametric curves with polynomial components, but they are
expressed in terms of the Bernstein basis rather than the traditional power
basis. In addition to being more numerically stable
[@Farouki1987, @Farouki1991, @Farouki1996], this basis allows "intuitive"
manipulation of geometric shapes by controlling a set of points rather
than via algebraic techniques.

B&#xe9;zier curves and surfaces have been widely used for decades
[@Farouki2012] in industrial design (e.g. shape), computer fonts and
graphics, mathematics (e.g. isoparametric elements in finite elements), and
many other fields.

This library provides support for

- 2D plotting
- 2D [intersection][3] (via both geometric
  [@Sederberg1989, @Sederberg1990, @Kim1998, @Sederberg1986] and
  algebraic [@JonssonVavasis, @Manocha:CSD-92-698] algorithms)
- Curve and surface subdivision [@Farouki1990]
- Degree-elevation and reduction
- Evaluation of points on curves / surfaces
- Determining parameters corresponding to a point on a [on a curve][4] or
  [on a surface][5] (i.e. the inverse of evaluation)
- Specialization / reparameterization
- Self-intersection / singularity check for 2D surfaces

-![Surface-surface intersection example](https://cdn.rawgit.com/dhermes/bezier/0.4.0/docs/images/test_surfaces6Q_and_7Q.png)

[1]: https://en.wikipedia.org/wiki/B%C3%A9zier_curve
[2]: https://en.wikipedia.org/wiki/B%C3%A9zier_triangle
[3]: https://bezier.readthedocs.io/en/0.4.0/curve-curve-intersection.html
[4]: https://bezier.readthedocs.io/en/0.4.0/reference/bezier.curve.html#bezier.curve.Curve.locate
[5]: https://bezier.readthedocs.io/en/0.4.0/reference/bezier.surface.html#bezier.surface.Surface.locate

# References
