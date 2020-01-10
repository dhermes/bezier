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

import functools

try:
    import sympy
except ImportError:  # pragma: NO COVER
    sympy = None


def require_sympy(wrapped):
    """Function decorator to require :mod:`sympy` to exist.

    Args:
        wrapped (Callable): A function to be wrapped.

    Returns:
        Callable: The wrapped function.
    """

    @functools.wraps(wrapped)
    def ensure_sympy(*args, **kwargs):
        if sympy is None:
            raise OSError("This function requires SymPy.")
        return wrapped(*args, **kwargs)

    return ensure_sympy


@require_sympy
def to_symbolic(nodes):
    """Convert a 2D NumPy array to a SymPy matrix of rational numbers.

    Args:
        nodes (numpy.ndarray): Nodes defining an object, a 2D NumPy array.

    Returns:
        sympy.Matrix: The nodes as a SymPy matrix of rational numbers.

    Raises:
        ValueError: If ``nodes`` is not 2D.
    """
    if nodes.ndim != 2:
        raise ValueError("Nodes must be 2-dimensional, not", nodes.ndim)

    return sympy.Matrix(
        [[sympy.Rational(value) for value in row] for row in nodes]
    )


@require_sympy
def curve_weights(degree, s):
    """Compute de Casteljau weights for a curve.

    .. note::

       This function could be optimized by introducing a cache. However, any
       code using SymPy is (for now) assumed not to be performance critical.

    Args:
        degree (int): The degree of a curve.
        s (sympy.Symbol): The sympy to be used in the weights.

    Returns:
        sympy.Matrix: The de Casteljau weights for a curve as a
        ``(degree + 1) x 1`` matrix.
    """
    return sympy.Matrix(
        [
            [sympy.binomial(degree, k) * s ** k * (1 - s) ** (degree - k)]
            for k in range(degree + 1)
        ]
    )


@require_sympy
def curve_as_polynomial(nodes):
    """Convert ``nodes`` into a SymPy polynomial array :math:`B(s)`.

    Args:
        nodes (numpy.ndarray): Nodes defining a B |eacute| zier curve.

    Returns:
        sympy.Matrix: The curve :math:`B(s)`.
    """
    nodes_sym = to_symbolic(nodes)
    degree = nodes_sym.cols - 1

    s = sympy.Symbol("s")
    b_polynomial = nodes_sym * curve_weights(degree, s)
    b_polynomial.simplify()

    factored = [value.factor() for value in b_polynomial]
    return sympy.Matrix(factored).reshape(*b_polynomial.shape)
