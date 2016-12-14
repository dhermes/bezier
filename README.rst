``bezier``
==========

    Helper for B |eacute| zier Curves, Triangles, and Higher Order Objects

|build| |coverage| |pypi| |versions| |docs|

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:

Install
-------

.. code-block:: console

   $ pip install --upgrade bezier

Getting Started
---------------

For example, to create a curve:

.. code-block:: python

   >>> nodes1 = np.array([
   ...     [0.0, 0.0],
   ...     [0.5, 1.0],
   ...     [1.0, 0.0],
   ... ])
   >>> curve1 = bezier.Curve(nodes1)

The intersection (points) between two curves can
also be determined:

.. code-block:: python

   >>> nodes2 = np.array([
   ...     [0.0 ,  0.0],
   ...     [0.25,  2.0],
   ...     [0.5 , -2.0],
   ...     [0.75,  2.0],
   ...     [1.0 ,  0.0],
   ... ])
   >>> curve2 = bezier.Curve(nodes2)
   >>> intersections = curve1.intersect(curve2)
   >>> intersections
   array([[ 0.311...,  0.428...],
          [ 0.688...,  0.428...],
          [ 0.      ,  0.      ],
          [ 1.      ,  0.      ]])

and then we can plot these curves (along with their
intersections):

.. code-block:: python

   >>> import matplotlib.pyplot as plt
   >>> import seaborn
   >>>
   >>> ax = curve1.plot(num_pts=256)
   >>> curve2.plot(num_pts=256, ax=ax)
   >>> ax.plot(intersections[:, 0], intersections[:, 1],
   ...         marker='o', linestyle='None', color='black')
   >>> ax.axis('scaled')
   >>> ax.set_xlim(-0.125, 1.125)
   >>> ax.set_ylim(-0.0625, 0.625)
   >>> plt.show()

.. image:: https://cdn.rawgit.com/dhermes/bezier/master/docs/images/test_curves1_and_13.png
   :align: center

Development
-----------

To work on adding a feature or to run the functional tests,
See the `DEVELOPMENT doc`_ for more information on how to get
started.

.. _DEVELOPMENT doc: https://github.com/dhermes/bezier/blob/master/DEVELOPMENT.rst

License
-------

Apache 2.0 - See `the LICENSE`_ for more information.

.. _the LICENSE: https://github.com/dhermes/bezier/blob/master/LICENSE

.. |docs| image:: https://readthedocs.org/projects/bezier/badge/?version=latest
   :target: http://bezier.readthedocs.io/en/latest/
   :alt: Documentation Status
.. |build| image:: https://circleci.com/gh/dhermes/bezier.svg?style=shield
   :target: https://circleci.com/gh/dhermes/bezier
   :alt: CirleCI Build
.. |pypi| image:: https://img.shields.io/pypi/v/bezier.svg
   :target: https://pypi.python.org/pypi/bezier
.. |versions| image:: https://img.shields.io/pypi/pyversions/bezier.svg
   :target: https://pypi.python.org/pypi/bezier
.. |coverage| image:: https://coveralls.io/repos/github/dhermes/bezier/badge.svg?branch=master
   :target: https://coveralls.io/github/dhermes/bezier?branch=master
