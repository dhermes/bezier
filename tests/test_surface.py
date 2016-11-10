# See the License for the specific language governing permissions and
# limitations under the License.
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


class Test__get_degree(unittest.TestCase):

    @staticmethod
    def _call_function_under_test(degree):
        from bezier import surface

        return surface._get_degree(degree)

    def test_valid(self):
        self.assertEqual(0, self._call_function_under_test(1))
        self.assertEqual(1, self._call_function_under_test(3))
        self.assertEqual(2, self._call_function_under_test(6))
        self.assertEqual(3, self._call_function_under_test(10))
        self.assertEqual(11, self._call_function_under_test(78))

    def test_invalid(self):
        with self.assertRaises(ValueError):
            self._call_function_under_test(2)

        with self.assertRaises(ValueError):
            self._call_function_under_test(9)


class TestSurface(unittest.TestCase):

    @staticmethod
    def _get_target_class():
        from bezier import surface

        return surface.Surface

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        import numpy as np

        nodes = np.array([
            [0.0, 0.0],
            [0.625, 0.5],
            [1.0, 0.75],
        ])
        surface = self._make_one(nodes)
        self.assertEqual(surface._degree, 1)
        self.assertEqual(surface._dimension, 2)
        self.assertIs(surface._nodes, nodes)

    def test_constructor_wrong_dimension(self):
        import numpy as np

        nodes = np.array([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

        nodes = np.zeros((2, 2, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test_constructor_bad_degree(self):
        import numpy as np

        nodes = np.array([
            [0.0, 0.0],
        ])
        with self.assertRaises(ValueError):
            self._make_one(nodes)
