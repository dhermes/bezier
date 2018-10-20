############################
Cython ``.pxd`` Declarations
############################

In addition to the header files, several ``cimport``-able ``.pxd``
Cython declaration files are provided:

.. testsetup:: show-pxd

   import os
   import textwrap

   import bezier


   class Path(object):
       """This class is a hack for Windows.

       It wraps a simple string but prints / repr-s it with Windows
       path separator converted to the standard *nix separator.

       This way doctests will succeed on Windows without modification.
       """

       def __init__(self, path):
           self.path = path

       def __repr__(self):
           posix_path = self.path.replace(os.path.sep, "/")
           return repr(posix_path)


   def sort_key(name):
       return name.lower().lstrip("_")


   def tree(directory, suffix=None):
       names = sorted(os.listdir(directory), key=sort_key)
       parts = []
       for name in names:
           path = os.path.join(directory, name)
           if os.path.isdir(path):
               sub_part = tree(path, suffix=suffix)
               if sub_part is not None:
                   # NOTE: We **always** use posix separator.
                   parts.append(name + "/")
                   parts.append(textwrap.indent(sub_part, "  "))
           else:
               if suffix is None or name.endswith(suffix):
                   parts.append(name)

       if parts:
           return "\n".join(parts)
       else:
           return None


   def print_tree(directory, suffix=None):
       if isinstance(directory, Path):
           # Make Windows act like posix.
           directory = directory.path
           separator = "/"
       else:
           separator = os.path.sep
       print(os.path.basename(directory) + separator)
       full_tree = tree(directory, suffix=suffix)
       print(textwrap.indent(full_tree, "  "))


   include_directory = bezier.get_include()
   bezier_directory = Path(os.path.dirname(include_directory))

.. doctest:: show-pxd

   >>> bezier_directory
   '.../site-packages/bezier'
   >>> print_tree(bezier_directory, suffix=".pxd")
   bezier/
     _curve.pxd
     _curve_intersection.pxd
     _helpers.pxd
     _status.pxd
     _surface.pxd
     _surface_intersection.pxd

For example, ``cimport bezier._curve`` will provide all the functions
in ``bezier/curve.h``.

.. toctree::
   :titlesonly:

   curve
   curve_intersection
   helpers
   status
   surface
   surface_intersection
