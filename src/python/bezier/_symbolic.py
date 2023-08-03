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
* Computing B |eacute| zier curve and triangle polynomial representations
* Implicitizing parametric equations

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""


def to_symbolic(nodes):
    """Convert a 2D NumPy array to a SymPy matrix of rational numbers.

    Args:
        nodes (numpy.ndarray): Nodes defining an object, a 2D NumPy array.

    Returns:
        sympy.Matrix: The nodes as a SymPy matrix of rational numbers.

    Raises:
        ValueError: If ``nodes`` is not 2D.
    """
    # NOTE: We import SymPy at runtime to avoid the import-time cost for users
    #       that don't want to do symbolic computation. The ``sympy`` import is
    #       a tad expensive.
    import sympy  # pylint: disable=import-outside-toplevel

    if nodes.ndim != 2:
        raise ValueError("Nodes must be 2-dimensional, not", nodes.ndim)

    return sympy.Matrix(
        [[sympy.Rational(value) for value in row] for row in nodes]
    )


def curve_weights(degree, s):
    """Compute de Casteljau weights for a curve.

    .. note::

       This function could be optimized by introducing a cache. However, any
       code using SymPy is (for now) assumed not to be performance critical.

    Args:
        degree (int): The degree of a curve.
        s (sympy.Symbol): The symbol to be used in the weights.

    Returns:
        sympy.Matrix: The de Casteljau weights for the curve as a
        ``(degree + 1) x 1`` matrix.
    """
    # NOTE: We import SymPy at runtime to avoid the import-time cost for users
    #       that don't want to do symbolic computation. The ``sympy`` import is
    #       a tad expensive.
    import sympy  # pylint: disable=import-outside-toplevel

    return sympy.Matrix(
        [
            [sympy.binomial(degree, k) * s**k * (1 - s) ** (degree - k)]
            for k in range(degree + 1)
        ]
    )


def curve_as_polynomial(nodes, degree):
    """Convert ``nodes`` into a SymPy polynomial array :math:`B(s)`.

    Args:
        nodes (numpy.ndarray): Nodes defining a B |eacute| zier curve.
        degree (int): The degree of the curve. This is assumed to
            correctly correspond to the number of ``nodes``.

    Returns:
        Tuple[sympy.Symbol, sympy.Matrix]: Pair of
        * The symbol ``s`` used in the polynomial
        * The curve :math:`B(s)`.
    """
    # NOTE: We import SymPy at runtime to avoid the import-time cost for users
    #       that don't want to do symbolic computation. The ``sympy`` import is
    #       a tad expensive.
    import sympy  # pylint: disable=import-outside-toplevel

    nodes_sym = to_symbolic(nodes)

    s = sympy.Symbol("s")
    b_polynomial = nodes_sym * curve_weights(degree, s)
    b_polynomial.simplify()

    factored = [value.factor() for value in b_polynomial]
    return s, sympy.Matrix(factored).reshape(*b_polynomial.shape)


def implicitize_2d(x_fn, y_fn, s):
    """Implicitize a 2D parametric curve.

    Args:
        x_fn (sympy.Expr): Function :math:`x(s)` in the curve.
        y_fn (sympy.Expr): Function :math:`y(s)` in the curve.
        s (sympy.Symbol): The symbol used to define ``x_fn`` and ``y_fn``.

    Returns:
        sympy.Expr: The implicitized function :math:`f(x, y)` such that the
        curve satisfies :math:`f(x(s), y(s)) = 0`.
    """
    # NOTE: We import SymPy at runtime to avoid the import-time cost for users
    #       that don't want to do symbolic computation. The ``sympy`` import is
    #       a tad expensive.
    import sympy  # pylint: disable=import-outside-toplevel

    x_sym, y_sym = sympy.symbols("x, y")
    return sympy.resultant(x_fn - x_sym, y_fn - y_sym, s).factor()


def implicitize_curve(nodes, degree):
    """Implicitize a 2D parametric curve, given the nodes.

    .. note::

       This function assumes (but does not check) that the caller has verified
       ``nodes`` represents a 2D curve. If it **does not**, this function will
       error out when tuple-unpacking the polynomial.

    Args:
        nodes (numpy.ndarray): Nodes defining a B |eacute| zier curve.
        degree (int): The degree of the curve. This is assumed to
            correctly correspond to the number of ``nodes``.

    Returns:
        sympy.Expr: The implicitized function :math:`f(x, y)` such that the
        curve satisfies :math:`f(x(s), y(s)) = 0`.
    """
    s, b_polynomial = curve_as_polynomial(nodes, degree)
    x_fn, y_fn = b_polynomial
    return implicitize_2d(x_fn, y_fn, s)


def triangle_weights(degree, s, t):
    """Compute de Casteljau weights for a triangle.

    .. note::

       This function could be optimized by introducing a cache. However, any
       code using SymPy is (for now) assumed not to be performance critical.

    Args:
        degree (int): The degree of a triangle.
        s (sympy.Symbol): The first symbol to be used in the weights.
        t (sympy.Symbol): The second symbol to be used in the weights.

    Returns:
        sympy.Matrix: The de Casteljau weights for the triangle as an ``N x 1``
        matrix, where ``N == (degree + 1)(degree + 2) / 2``.
    """
    # NOTE: We import SymPy at runtime to avoid the import-time cost for users
    #       that don't want to do symbolic computation. The ``sympy`` import is
    #       a tad expensive.
    import sympy  # pylint: disable=import-outside-toplevel

    lambda1 = 1 - s - t
    lambda2 = s
    lambda3 = t

    values = []
    for k in range(degree + 1):
        coeff_k = sympy.binomial(degree, k)
        # N! / i! j! k! = N! / [k! (N - k)!] (i + j)! / [i! j!]
        for j in range(degree - k + 1):
            i = degree - j - k
            coeff = coeff_k * sympy.binomial(degree - k, j)
            values.append(coeff * lambda1**i * lambda2**j * lambda3**k)

    return sympy.Matrix(values).reshape(len(values), 1)


def triangle_as_polynomial(nodes, degree):
    """Convert ``nodes`` into a SymPy polynomial array :math:`B(s, t)`.

    Args:
        nodes (numpy.ndarray): Nodes defining a B |eacute| zier triangle.
        degree (int): The degree of the triangle. This is assumed to
            correctly correspond to the number of ``nodes``.

    Returns:
        Tuple[sympy.Symbol, sympy.Symbol, sympy.Matrix]: Triple of
        * The symbol ``s`` used in the polynomial
        * The symbol ``t`` used in the polynomial
        * The triangle :math:`B(s, t)`.
    """
    # NOTE: We import SymPy at runtime to avoid the import-time cost for users
    #       that don't want to do symbolic computation. The ``sympy`` import is
    #       a tad expensive.
    import sympy  # pylint: disable=import-outside-toplevel

    nodes_sym = to_symbolic(nodes)

    s, t = sympy.symbols("s, t")
    b_polynomial = nodes_sym * triangle_weights(degree, s, t)
    b_polynomial.simplify()

    factored = [value.factor() for value in b_polynomial]
    return s, t, sympy.Matrix(factored).reshape(*b_polynomial.shape)


def implicitize_3d(x_fn, y_fn, z_fn, s, t):
    """Implicitize a 3D parametric triangle.

    Args:
        x_fn (sympy.Expr): Function :math:`x(s)` in the triangle.
        y_fn (sympy.Expr): Function :math:`y(s)` in the triangle.
        z_fn (sympy.Expr): Function :math:`z(s)` in the triangle.
        s (sympy.Symbol): The first symbol used to define ``x_fn``, ``y_fn``
            and ``z_fn``.
        t (sympy.Symbol): The second symbol used to define ``x_fn``, ``y_fn``
            and ``z_fn``.

    Returns:
        sympy.Expr: The implicitized function :math:`f(x, y, z)` such that the
        triangle satisfies :math:`f(x(s, t), y(s, t), z(s, t)) = 0`.
    """
    # NOTE: We import SymPy at runtime to avoid the import-time cost for users
    #       that don't want to do symbolic computation. The ``sympy`` import is
    #       a tad expensive.
    import sympy  # pylint: disable=import-outside-toplevel

    x_sym, y_sym, z_sym = sympy.symbols("x, y, z")

    f_xy = sympy.resultant(x_fn - x_sym, y_fn - y_sym, s)
    f_yz = sympy.resultant(y_fn - x_sym, z_fn - z_sym, s)
    return sympy.resultant(f_xy, f_yz, t).factor()


def implicitize_triangle(nodes, degree):
    """Implicitize a 3D parametric triangle, given the nodes.

    .. note::

       This function assumes (but does not check) that the caller has verified
       ``nodes`` represents a 3D triangle. If it **does not**, this function
       will error out when tuple-unpacking the polynomial.

    Args:
        nodes (numpy.ndarray): Nodes defining a B |eacute| zier triangle.
        degree (int): The degree of the triangle. This is assumed to
            correctly correspond to the number of ``nodes``.

    Returns:
        sympy.Expr: The implicitized function :math:`f(x, y, z)` such that the
        triangle satisfies :math:`f(x(s, t), y(s, t), z(s, t)) = 0`.
    """
    s, t, b_polynomial = triangle_as_polynomial(nodes, degree)
    x_fn, y_fn, z_fn = b_polynomial
    return implicitize_3d(x_fn, y_fn, z_fn, s, t)
