Native Libraries
================

``bezier`` has optional speedups implemented in `Fortran`_.
These are incorporated into the Python interface via
`Cython`_.

.. _Fortran: https://en.wikipedia.org/wiki/Fortran
.. _Cython: https://cython.readthedocs.io/

The subroutines provided there can be called from Fortran,
C, C++, Cython and any other language that can invoke
a foreign C function (e.g. `Go`_ via ``cgo``).

.. _Go: https://golang.org

After ``bezier`` has been installed **with** these speedups,
the library provides helpers to make it easier to build
code that depends on them.

C Headers
---------

The C headers for ``libbezier`` will be included in the source-tree

.. testsetup:: show-headers, show-lib, show-pxd

   import os
   import platform
   import textwrap

   import bezier


   PLATFORM_SYSTEM = platform.system().lower()


   class Path(object):

       def __init__(self, path):
           # Hack to make the doctest work just fine on Windows.
           self.path = path

       def __repr__(self):
           posix_path = self.path.replace(os.path.sep, '/')
           return repr(posix_path)


   def sort_key(name):
       return name.lower().lstrip('_')


   def tree(directory, suffix=None):
       names = sorted(os.listdir(directory), key=sort_key)
       parts = []
       for name in names:
           path = os.path.join(directory, name)
           if os.path.isdir(path):
               sub_part = tree(path, suffix=suffix)
               if sub_part is not None:
                   # NOTE: We **always** use posix separator.
                   parts.append(name + '/')
                   parts.append(textwrap.indent(sub_part, '  '))
           else:
               if suffix is None or name.endswith(suffix):
                   parts.append(name)

       if parts:
           return '\n'.join(parts)
       else:
           return None


   def print_tree(directory, suffix=None):
       assert isinstance(directory, Path)
       directory = directory.path
       # NOTE: We **always** use posix separator.
       print(os.path.basename(directory) + '/')

       if directory.endswith('lib') and PLATFORM_SYSTEM == 'windows':
           # NOTE: This is a (temporary) hack so that doctests can pass
           #       on Windows even though the directory contents differ.
           assert not os.path.exists(directory)
           full_tree = 'libbezier.a'
       else:
           full_tree = tree(directory, suffix=suffix)

       print(textwrap.indent(full_tree, '  '))


   def parent_directory(directory):
       assert isinstance(directory, Path)
       return Path(os.path.dirname(directory.path))


   # Monkey-patch functions to return a ``Path``.
   original_get_include = bezier.get_include
   original_get_lib = bezier.get_lib

   def get_include():
       return Path(original_get_include())

   def get_lib():
       return Path(original_get_lib())

   bezier.get_include = get_include
   bezier.get_lib = get_lib

   # Allow this value to be re-used.
   include_directory = get_include()

.. doctest:: show-headers

   >>> include_directory = bezier.get_include()
   >>> include_directory
   '.../site-packages/bezier/include'
   >>> print_tree(include_directory)
   include/
     bezier/
       curve.h
       curve_intersection.h
       helpers.h
       surface.h
     bezier.h

.. testcleanup:: show-headers, show-lib, show-pxd

   # Restore the monkey-patched functions.
   bezier.get_include = original_get_include
   bezier.get_lib = original_get_lib

Note that this includes a catch-all ``bezier.h`` that just includes all of
the headers.

Cython ``.pxd`` Declarations
----------------------------

In addition to the header files, several ``cimport``-able ``.pxd``
Cython declaration files are provided:

.. doctest:: show-pxd

   >>> bezier_directory = parent_directory(include_directory)
   >>> bezier_directory
   '.../site-packages/bezier'
   >>> print_tree(bezier_directory, suffix='.pxd')
   bezier/
     _curve.pxd
     _curve_intersection.pxd
     _helpers.pxd
     _surface.pxd

For example, ``cimport bezier._curve`` will provide all the functions
in ``bezier/curve.h``.

Static Library
--------------

The actual library ``libbezier`` is included as a single static library
(a ``.lib`` file on Windows and a ``.a`` file elsewhere):

.. doctest:: show-lib

   >>> lib_directory = bezier.get_lib()
   >>> lib_directory
   '.../site-packages/bezier/lib'
   >>> print_tree(lib_directory)
   lib/
     libbezier.a

.. note::

   A static library is used (rather than a shared or dynamic library)
   because the "final" install location of the Python package is not
   dependable. Even on the same machine with the same operating system,
   ``bezier`` can be installed in virtual environments, in different
   Python versions, as an egg or wheel, and so on.

.. warning::

   When ``bezier`` is installed via `pip`_, it will likely be installed
   from a `Python wheel`_. These wheels will be pre-built and the Fortran
   extensions will be compiled with `GNU Fortran`_ (``gfortran``). As a
   result, ``libbezier`` will depend on ``libgfortran``.

   This can be problematic due to version conflicts, ABI incompatibility,
   a desire to use a different Fortran compiler (e.g. ``ifort``) and a host
   of other reasons. Some of the standard `tooling`_ for building wheels
   will try to address this by adding a ``bezier/.libs`` directory with a
   version of ``libgfortran`` that is compatible with ``libbezier``, e.g.

   .. code-block:: rest

      .../site-packages/bezier/.libs/libgfortran-ed201abd.so.3.0.0

   If present, this directory can be used when linking. If that is not
   feasible, then ``bezier`` can be built from source via:

   .. code-block:: console

       python setup.py build
       python setup.py build --fcompiler=${FC}

.. _pip: https://pip.pypa.io
.. _Python wheel: https://wheel.readthedocs.io
.. _GNU Fortran: https://gcc.gnu.org/fortran/
.. _tooling: https://github.com/pypa/auditwheel

Building a Python Extension
---------------------------

To incorporate ``libbezier`` into a Python extension, either via
Cython, C, C++ or some other means, simply include the header
and library directories:

.. testsetup:: setup-extension

   import bezier

.. doctest:: setup-extension

   >>> import setuptools
   >>>
   >>> extension = setuptools.Extension(
   ...     'wrapper',
   ...     ['wrapper.c'],
   ...     include_dirs=[
   ...         bezier.get_include(),
   ...     ],
   ...     libraries=['bezier'],
   ...     library_dirs=[
   ...         bezier.get_lib(),
   ...     ],
   ... )
   >>> extension
   <setuptools.extension.Extension('wrapper') at 0x...>
