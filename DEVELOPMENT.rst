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

3. The feature must (should) be documented with helpful ``doctest``
   examples wherever relevant.

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

.. _pytest: http://docs.pytest.org/en/latest/

.. code-block:: console

   $ python2.7 -m py.test tests/
   $ python3.4 -m py.test tests/
   $ python3.5 -m py.test tests/

Slow Tests
==========

To run unit tests without tests that have been (explicitly)
marked slow, use the ``--ignore-slow``.

.. code-block:: console

   $ python -m py.test tests/ --ignore-slow

These slow tests have been identified via

.. code-block:: console

   $ python -m py.test tests/ --durations=10

and then marked with ``pytest.mark.skipif``.
