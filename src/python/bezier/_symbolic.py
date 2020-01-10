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

"""Helper for doing symbolic algebra with data.

Includes functions for:

* Converting floating point arrays to rational SymPy matrices
* Computing B |eacute| zier curve and surface polynomial representations
* Implicitizing parametric equations

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

try:
    import sympy
except ImportError:  # pragma: NO COVER
    sympy = None


def to_symbolic(nodes):
    """Convert a 2D NumPy array to a SymPy matrix of rational numbers.

    Args:
        nodes (numpy.ndarray): Nodes defining an object, a 2D NumPy array.

    Returns:
        sympy.Matrix: The nodes as a SymPy matrix of rational numbers.

    Raises:
        OSError: If SymPy is not installed.
        ValueError: If ``nodes`` is not 2D.
    """
    if sympy is None:
        raise OSError("This function requires SymPy.")

    if nodes.ndim != 2:
        raise ValueError("Nodes must be 2-dimensional, not", nodes.ndim)

    return sympy.Matrix(
        [[sympy.Rational(value) for value in row] for row in nodes]
    )
