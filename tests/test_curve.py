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


class TestCurve(unittest.TestCase):

    @staticmethod
    def _get_target_class():
        from bezier import curve
        return curve.Curve

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
        curve = self._make_one(nodes)
        self.assertEqual(curve._degree, 2)
        self.assertEqual(curve._dimension, 2)
        self.assertIs(curve._nodes, nodes)

    def test_constructor_bad_input(self):
        import numpy as np

        nodes = np.array([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

        nodes = np.zeros((2, 2, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test_degree_property(self):
        import numpy as np

        degree = 6
        num_nodes = degree + 1
        nodes = np.zeros((num_nodes, 2))
        curve = self._make_one(nodes)
        self.assertEqual(curve.degree, degree)

    def test_dimension_property(self):
        import numpy as np

        dimension = 4
        nodes = np.zeros((3, dimension))
        curve = self._make_one(nodes)
        self.assertEqual(curve.dimension, dimension)

    def test_evaluate(self):
        import numpy as np

        s = 0.25
        nodes = np.array([
            [0.0, 0.0],
            [0.5, 0.5],
            [1.0, 1.25],
        ])
        curve = self._make_one(nodes)
        expected = np.array([0.25,  0.265625])
        result = curve.evaluate(s)
        self.assertTrue(np.all(expected == result))
