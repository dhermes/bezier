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

import unittest

import numpy as np


SPACING = np.spacing  # pylint: disable=no-member


class Test_compute_implicit_line(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _clipping

        return _clipping.compute_implicit_line(nodes)

    def test_no_rounding(self):
        nodes = np.asfortranarray([
            [1.0, 5.0],
            [2.0, 2.0],
        ])
        coeff_a, coeff_b, coeff_c = self._call_function_under_test(nodes)
        self.assertEqual(coeff_a, 0.0)
        self.assertEqual(coeff_b, 1.0)
        self.assertEqual(coeff_c, -2.0)

    def test_rational_length(self):
        nodes = np.asfortranarray([
            [3.0, 7.0],
            [2.0, 5.0],
        ])
        coeff_a, coeff_b, coeff_c = self._call_function_under_test(nodes)
        self.assertEqual(coeff_a, -3.0 / 5.0)
        self.assertEqual(coeff_b, 4.0 / 5.0)
        self.assertEqual(coeff_c, 1.0 / 5.0)

    def test_irrational_length(self):
        nodes = np.asfortranarray([
            [4.0, 5.0],
            [7.0, 8.0],
        ])
        coeff_a, coeff_b, coeff_c = self._call_function_under_test(nodes)
        sqrt_half = np.sqrt(0.5)
        delta = SPACING(sqrt_half)
        self.assertAlmostEqual(coeff_a, -sqrt_half, delta=delta)
        self.assertAlmostEqual(coeff_b, sqrt_half, delta=delta)
        self.assertAlmostEqual(coeff_c, -3 * sqrt_half, delta=4 * delta)


class Test_compute_fat_line(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _clipping

        return _clipping.compute_fat_line(nodes)

    def test_line(self):
        nodes = np.asfortranarray([
            [1.0, 5.0],
            [2.0, 2.0],
        ])
        result = self._call_function_under_test(nodes)
        coeff_a, coeff_b, coeff_c, d_min, d_max = result
        self.assertEqual(coeff_a, 0.0)
        self.assertEqual(coeff_b, 1.0)
        self.assertEqual(coeff_c, -2.0)
        self.assertEqual(d_min, 0.0)
        self.assertEqual(d_max, 0.0)

    def test_quadratic(self):
        nodes = np.asfortranarray([
            [0.0, 1.0, 0.0],
            [0.0, 1.0, 2.0],
        ])
        result = self._call_function_under_test(nodes)
        coeff_a, coeff_b, coeff_c, d_min, d_max = result
        self.assertEqual(coeff_a, -1.0)
        self.assertEqual(coeff_b, 0.0)
        self.assertEqual(coeff_c, 0.0)
        self.assertEqual(d_min, -1.0)
        self.assertEqual(d_max, 0.0)
