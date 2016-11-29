###########
Development
###########


.. contents:: Here are some pointers for hacking on ``bezier``.

***************
Adding Features
***************

In order to add a feature to ``bezier``:

1. **File an issue** to discuss the change

2. **Add tests**: The feature must work fully on the following
   CPython versions: 2.7, 3.4 and 3.5 on both UNIX and Windows.
   In addition, the feature should have 100% line coverage.

3. **Documentation**: The feature must (should) be documented with
   helpful `doctest`_ examples wherever relevant.

.. _doctest: http://www.sphinx-doc.org/en/stable/ext/doctest.html

******************
Running Unit Tests
******************

We recommend using ``tox`` to run unit tests:

.. code-block:: console

   $ tox -e py27
   $ tox -e py34
   $ tox -e py35

However, `pytest`_ can be used directly (though it won't
manage dependencies):

.. _pytest: http://docs.pytest.org/en/stable/

.. code-block:: console

   $ python2.7 -m py.test tests/
   $ python3.4 -m py.test tests/
   $ python3.5 -m py.test tests/

Test Coverage
=============

``bezier`` has 100% `line coverage`_. The coverage is checked
on every build and uploaded to `coveralls.io`_ via the
``COVERALLS_REPO_TOKEN`` environment variable set in
the `CircleCI environment`_.

.. line coverage: https://coveralls.io/github/dhermes/bezier
.. coveralls.io: https://coveralls.io/
.. CircleCI environment: https://circleci.com/gh/dhermes/bezier/edit#env-vars

Slow Tests
==========

To run unit tests without tests that have been (explicitly)
marked slow, use the ``--ignore-slow`` flag:

.. code-block:: console

   $ python -m py.test tests/ --ignore-slow

These slow tests have been identified via:

.. code-block:: console

   $ python -m py.test tests/ --durations=10

and then marked with ``pytest.mark.skipif``.

**********************
Deploying New Versions
**********************

New versions are deployed to `PyPI`_ automatically every time
a new tag is pushed. To allow `twine`_ to authenticate (which
is needed to upload) the ``TWINE_USERNAME`` and ``TWINE_PASSWORD``
environment variables are set in the `CircleCI environment`_.

.. _PyPI: https://pypi.python.org/pypi/bezier
.. twine: https://packaging.python.org/distributing/
