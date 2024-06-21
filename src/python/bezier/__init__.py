# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

r"""Helper for B |eacute| zier Curves, Triangles, and Higher Order Objects.

Intended to perform basic operations on B |eacute| zier objects such
as intersections, length / area / etc. computations, subdivision,
implicitization and other relevant information.

Plotting utilities are also provided.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

# NOTE: ``__config__`` **must** be the first import because it (may)
#       modify the search path used to locate shared libraries.
from bezier import __config__
from bezier.curve import Curve
from bezier.curved_polygon import CurvedPolygon
from bezier.hazmat.helpers import UnsupportedDegree
from bezier.triangle import Triangle

try:
    import bezier._speedup  # noqa: F401

    _HAS_SPEEDUP = True
except ImportError as exc:  # pragma: NO COVER
    __config__.handle_import_error(exc, "_speedup")
    _HAS_SPEEDUP = False
# NOTE: The ``__version__`` and ``__author__`` are hard-coded here, rather
#       than using ``importlib.metadata.distribution("bezier").version``
#       and related. This is **entirely** to accommodate builds where
#       ``bezier`` is imported from source (and not installed).
__author__ = "Danny Hermes"
__version__ = "2024.6.20"
"""str: The current version of :mod:`bezier`."""
__all__ = [
    "__author__",
    "__version__",
    "Curve",
    "CurvedPolygon",
    "Triangle",
    "UnsupportedDegree",
]
