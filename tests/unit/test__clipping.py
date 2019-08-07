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
        nodes = np.asfortranarray([[1.0, 5.0], [2.0, 2.0]])
        coeff_a, coeff_b, coeff_c = self._call_function_under_test(nodes)
        self.assertEqual(coeff_a, 0.0)
        self.assertEqual(coeff_b, 1.0)
        self.assertEqual(coeff_c, -2.0)

    def test_rational_length(self):
        nodes = np.asfortranarray([[3.0, 7.0], [2.0, 5.0]])
        coeff_a, coeff_b, coeff_c = self._call_function_under_test(nodes)
        self.assertEqual(coeff_a, -3.0 / 5.0)
        self.assertEqual(coeff_b, 4.0 / 5.0)
        self.assertEqual(coeff_c, 1.0 / 5.0)

    def test_irrational_length(self):
        nodes = np.asfortranarray([[4.0, 5.0], [7.0, 8.0]])
        coeff_a, coeff_b, coeff_c = self._call_function_under_test(nodes)
        sqrt_half = np.sqrt(0.5)
        delta = SPACING(sqrt_half)  # pylint: disable=assignment-from-no-return
        self.assertAlmostEqual(coeff_a, -sqrt_half, delta=delta)
        self.assertAlmostEqual(coeff_b, sqrt_half, delta=delta)
        self.assertAlmostEqual(coeff_c, -3 * sqrt_half, delta=4 * delta)


class Test_compute_fat_line(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes):
        from bezier import _clipping

        return _clipping.compute_fat_line(nodes)

    def test_line(self):
        nodes = np.asfortranarray([[1.0, 5.0], [2.0, 2.0]])
        result = self._call_function_under_test(nodes)
        coeff_a, coeff_b, coeff_c, d_min, d_max = result
        self.assertEqual(coeff_a, 0.0)
        self.assertEqual(coeff_b, 1.0)
        self.assertEqual(coeff_c, -2.0)
        self.assertEqual(d_min, 0.0)
        self.assertEqual(d_max, 0.0)

    def test_quadratic(self):
        nodes = np.asfortranarray([[0.0, 1.0, 0.0], [0.0, 1.0, 2.0]])
        result = self._call_function_under_test(nodes)
        coeff_a, coeff_b, coeff_c, d_min, d_max = result
        self.assertEqual(coeff_a, -1.0)
        self.assertEqual(coeff_b, 0.0)
        self.assertEqual(coeff_c, 0.0)
        self.assertEqual(d_min, -1.0)
        self.assertEqual(d_max, 0.0)

    def test_many_interior(self):
        nodes = np.asfortranarray(
            [[0.0, 1.0, 2.0, 3.0, 4.0], [0.0, 4.0, -4.0, 2.0, 0.0]]
        )
        result = self._call_function_under_test(nodes)
        coeff_a, coeff_b, coeff_c, d_min, d_max = result
        self.assertEqual(coeff_a, 0.0)
        self.assertEqual(coeff_b, 1.0)
        self.assertEqual(coeff_c, 0.0)
        self.assertEqual(d_min, -4.0)
        self.assertEqual(d_max, 4.0)


class Test__update_parameters(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(s_min, s_max, start0, end0, start1, end1):
        from bezier import _clipping

        return _clipping._update_parameters(
            s_min, s_max, start0, end0, start1, end1
        )

    def test_parallel(self):
        from bezier import _clipping

        start0 = np.asfortranarray([0.0, 0.0])
        end0 = np.asfortranarray([1.0, 0.0])
        start1 = np.asfortranarray([0.0, -1.0])
        end1 = np.asfortranarray([1.0, -1.0])
        with self.assertRaises(NotImplementedError) as exc_info:
            self._call_function_under_test(
                None, None, start0, end0, start1, end1
            )
        expected_args = (_clipping.NO_PARALLEL,)
        self.assertEqual(exc_info.exception.args, expected_args)

    def test_t_outside(self):
        start0 = np.asfortranarray([0.0, -1.0])
        end0 = np.asfortranarray([1.0, -1.0])
        start1 = np.asfortranarray([0.5, 0.0])
        end1 = np.asfortranarray([1.0, 2.0])
        s_min, s_max = self._call_function_under_test(
            2.0, -1.0, start0, end0, start1, end1
        )
        self.assertEqual(s_min, 2.0)
        self.assertEqual(s_max, -1.0)

    def _update_helper(self, s_min, s_max):
        start0 = np.asfortranarray([0.0, 2.0])
        end0 = np.asfortranarray([1.0, 2.0])
        start1 = np.asfortranarray([0.0, 1.0])
        end1 = np.asfortranarray([0.5, 3.0])
        return self._call_function_under_test(
            s_min, s_max, start0, end0, start1, end1
        )

    def test_update_s_min(self):
        s_min, s_max = self._update_helper(1.0, -1.0)
        self.assertEqual(s_min, 0.25)
        self.assertEqual(s_max, -1.0)

    def test_update_s_max(self):
        s_min, s_max = self._update_helper(0.125, -1.0)
        self.assertEqual(s_min, 0.125)
        self.assertEqual(s_max, 0.25)

    def test_s_not_updated(self):
        s_min, s_max = self._update_helper(0.125, 0.5)
        self.assertEqual(s_min, 0.125)
        self.assertEqual(s_max, 0.5)


class Test__check_parameter_range(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(s_min, s_max):
        from bezier import _clipping

        return _clipping._check_parameter_range(s_min, s_max)

    def test_default_both(self):
        from bezier import _clipping

        s_min, s_max = self._call_function_under_test(
            _clipping.DEFAULT_S_MIN, _clipping.DEFAULT_S_MAX
        )
        self.assertEqual(s_min, 0.0)
        self.assertEqual(s_max, 1.0)

    def test_default_max(self):
        from bezier import _clipping

        s_min, s_max = self._call_function_under_test(
            0.25, _clipping.DEFAULT_S_MAX
        )
        self.assertEqual(s_min, 0.25)
        self.assertEqual(s_max, 0.25)

    def test_both_set(self):
        s_min, s_max = self._call_function_under_test(0.25, 0.75)
        self.assertEqual(s_min, 0.25)
        self.assertEqual(s_max, 0.75)


class Test_clip_range(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _clipping

        return _clipping.clip_range(nodes1, nodes2)

    def test_it(self):
        nodes1 = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 2.0, 0.0]])
        nodes2 = np.asfortranarray([[0.0, 2.0, 0.0], [-1.0, 1.0, 3.0]])
        start_s, end_s = self._call_function_under_test(nodes1, nodes2)
        self.assertEqual(start_s, 0.25)
        self.assertEqual(end_s, 0.75)

    def test_parallel(self):
        from bezier import _clipping

        nodes1 = np.asfortranarray([[0.0, 1.0, 2.0], [1.0, 3.0, 1.0]])
        nodes2 = np.asfortranarray(
            [[0.0, 0.5, 1.0, 1.5, 2.0], [0.0, 4.0, 4.0, 4.0, 0.0]]
        )
        with self.assertRaises(NotImplementedError) as exc_info:
            self._call_function_under_test(nodes1, nodes2)
        expected_args = (_clipping.NO_PARALLEL,)
        self.assertEqual(exc_info.exception.args, expected_args)
