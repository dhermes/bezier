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

# NOTE: This module is written in the ``pytest`` style, i.e. it does not use
#       ``unittest``. The other test modules in this project use ``unittest``,
#       they were written before ``pytest`` was used in this project (the
#       original test runner was ``nose``).

import unittest
import unittest.mock

import numpy as np
import pytest

try:
    import sympy
except ImportError:  # pragma: NO COVER
    sympy = None


def sympy_equal(value1, value2):
    difference = (value1 - value2).expand()
    return difference == 0


def sympy_matrix_equal(value1, value2):
    if value1.shape != value2.shape:  # pragma: NO COVER
        return False

    difference = (value1 - value2).expand()
    return difference == sympy.zeros(*value1.shape)


class Test_to_symbolic:
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _symbolic

        return _symbolic.to_symbolic(nodes)

    @unittest.mock.patch("bezier._symbolic.sympy", new=None)
    def test_sympy_missing(self):
        with pytest.raises(OSError) as exc_info:
            self._call_function_under_test(None)

        exc_args = exc_info.value.args
        assert exc_args == ("This function requires SymPy.",)

    def test_bad_shape(self):
        nodes = np.empty((1, 1, 1), order="F")
        with pytest.raises(ValueError) as exc_info:
            self._call_function_under_test(nodes)

        exc_args = exc_info.value.args
        assert exc_args == ("Nodes must be 2-dimensional, not", 3)

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_success(self):
        nodes = np.asfortranarray([[1.0, 5.5], [2.0, 2.0]])
        nodes_sym = self._call_function_under_test(nodes)
        expected = sympy.Matrix([[1, sympy.Rational(11, 2)], [2, 2],])
        assert nodes_sym == expected


class Test_curve_weights:
    @staticmethod
    def _call_function_under_test(degree, s):
        from bezier import _symbolic

        return _symbolic.curve_weights(degree, s)

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_it(self):
        t = sympy.Symbol("t")
        weights = self._call_function_under_test(3, t)
        expected = sympy.Matrix(
            [
                [
                    (1 - t) ** 3,
                    3 * t * (1 - t) ** 2,
                    3 * t ** 2 * (1 - t),
                    t ** 3,
                ]
            ]
        ).T
        assert sympy_matrix_equal(weights, expected)


class Test_curve_as_polynomial:
    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _symbolic

        return _symbolic.curve_as_polynomial(nodes, degree)

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_it(self):
        nodes = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, 1.0, 0.0]])
        s, b_polynomial = self._call_function_under_test(nodes, 2)
        expected = sympy.Matrix([s, 2 * s * (1 - s)])
        assert sympy_matrix_equal(b_polynomial, expected)


class Test_implicitize_2d:
    @staticmethod
    def _call_function_under_test(x_fn, y_fn, s):
        from bezier import _symbolic

        return _symbolic.implicitize_2d(x_fn, y_fn, s)

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_it(self):
        s, x_sym, y_sym = sympy.symbols("s, x, y")
        x_fn = 1 - s ** 2
        y_fn = 1 + s + s ** 2
        f_polynomial = self._call_function_under_test(x_fn, y_fn, s)

        expected = (
            x_sym ** 2
            + 2 * x_sym * y_sym
            - 3 * x_sym
            + y_sym ** 2
            - 4 * y_sym
            + 3
        )
        assert sympy_equal(f_polynomial, expected)


class Test_surface_weights:
    @staticmethod
    def _call_function_under_test(degree, s, t):
        from bezier import _symbolic

        return _symbolic.surface_weights(degree, s, t)

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_it(self):
        s, t = sympy.symbols("s, t")
        weights = self._call_function_under_test(2, s, t)
        expected = sympy.Matrix(
            [
                [
                    (1 - s - t) ** 2,
                    2 * (1 - s - t) * s,
                    s ** 2,
                    2 * (1 - s - t) * t,
                    2 * s * t,
                    t ** 2,
                ]
            ]
        ).T
        assert sympy_matrix_equal(weights, expected)


class Test_surface_as_polynomial:
    @staticmethod
    def _call_function_under_test(nodes, degree):
        from bezier import _symbolic

        return _symbolic.surface_as_polynomial(nodes, degree)

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_it(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0, 1.0], [0.0, 0.0, 1.0], [1.0, 1.0, 1.0]]
        )
        s, t, b_polynomial = self._call_function_under_test(nodes, 1)
        expected = sympy.Matrix([s + t, t, 1])
        assert sympy_matrix_equal(b_polynomial, expected)


class Test_implicitize_3d:
    @staticmethod
    def _call_function_under_test(x_fn, y_fn, z_fn, s, t):
        from bezier import _symbolic

        return _symbolic.implicitize_3d(x_fn, y_fn, z_fn, s, t)

    @unittest.skipIf(sympy is None, "SymPy not installed")
    def test_it(self):
        s, t, x_sym, y_sym, z_sym = sympy.symbols("s, t, x, y, z")
        x_fn = s * t
        y_fn = s + t
        z_fn = s ** 2 + t ** 2
        f_polynomial = self._call_function_under_test(x_fn, y_fn, z_fn, s, t)

        expected = (
            x_sym ** 4
            - 2 * x_sym ** 3 * y_sym
            + 2 * x_sym ** 2 * y_sym ** 2
            - 4 * x_sym ** 2 * y_sym
            - 2 * x_sym ** 2 * z_sym
            + 4 * x_sym ** 2
            + 2 * x_sym * y_sym * z_sym
            + 4 * x_sym * z_sym
            - 2 * y_sym ** 2 * z_sym
            + z_sym ** 2
        )
        assert sympy_equal(f_polynomial, expected)
