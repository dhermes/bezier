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


class TestBase(unittest.TestCase):

    @staticmethod
    def _get_target_class():
        from bezier import _base

        return _base.Base

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        import numpy as np

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
        ])
        shape = self._make_one(nodes)
        self.assertEqual(shape._degree, 3)
        self.assertEqual(shape._dimension, 2)
        self.assertIs(shape._nodes, nodes)

    def test_virtual__get_degree(self):
        klass = self._get_target_class()
        num_nodes = 41
        self.assertEqual(num_nodes, klass._get_degree(num_nodes))

    def test_constructor_wrong_dimension(self):
        import numpy as np

        nodes = np.array([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

        nodes = np.zeros((2, 2, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test_constructor_insufficient_nodes(self):
        import numpy as np

        nodes = np.zeros((0, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test___repr__(self):
        import numpy as np

        nodes = np.zeros((4, 3))
        shape = self._make_one(nodes)
        expected = '<Base (degree=4, dimension=3)>'
        self.assertEqual(repr(shape), expected)

    def test_degree_property(self):
        import numpy as np

        degree = 6
        nodes = np.zeros((degree, 2))
        shape = self._make_one(nodes)
        self.assertEqual(shape.degree, degree)
        self.assertEqual(shape._degree, degree)

    def test_dimension_property(self):
        import numpy as np

        dimension = 4
        nodes = np.zeros((3, dimension))
        shape = self._make_one(nodes)
        self.assertEqual(shape.dimension, dimension)
        self.assertEqual(shape._dimension, dimension)

    def test_nodes_property(self):
        import numpy as np

        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        shape = self._make_one(nodes)
        self.assertTrue(np.all(shape.nodes == nodes))
        self.assertIsNot(shape.nodes, nodes)
