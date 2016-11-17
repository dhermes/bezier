# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import unittest

import numpy as np


class Test_bbox_intersect(unittest.TestCase):

    UNIT_SQUARE = np.array([
        [0.0, 0.0],
        [1.0, 0.0],
        [1.0, 1.0],
        [0.0, 1.0],
    ])

    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _intersection_helpers

        return _intersection_helpers.bbox_intersect(nodes1, nodes2)

    def test_intersect(self):
        nodes = self.UNIT_SQUARE + np.array([[0.5, 0.5]])
        self.assertTrue(self._call_function_under_test(
            self.UNIT_SQUARE, nodes))

    def test_far_apart(self):
        nodes = self.UNIT_SQUARE + np.array([[100.0, 100.0]])
        self.assertFalse(self._call_function_under_test(
            self.UNIT_SQUARE, nodes))

    def test_tangent(self):
        nodes = self.UNIT_SQUARE + np.array([[1.0, 0.0]])
        self.assertFalse(self._call_function_under_test(
            self.UNIT_SQUARE, nodes))

    def test_almost_tangent(self):
        x_val = 1.0 + np.spacing(1.0)
        nodes = self.UNIT_SQUARE + np.array([[x_val, 0.0]])
        self.assertFalse(self._call_function_under_test(
            self.UNIT_SQUARE, nodes))
