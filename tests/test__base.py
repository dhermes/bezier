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

import mock
import numpy as np

from tests import utils


class TestBase(utils.NumPyTestCase):

    @staticmethod
    def _get_target_class():
        from bezier import _base

        return _base.Base

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
        ])
        shape = self._make_one(nodes)
        self.assertEqual(shape._degree, 3)
        self.assertEqual(shape._dimension, 2)
        self.assertIsNot(shape._nodes, nodes)
        self.assertEqual(shape._nodes, nodes)

    def test_constructor_without_copy(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 1.0],
            [2.0, 3.0],
        ])
        shape = self._make_one(nodes, _copy=False)
        self.assertEqual(shape._degree, 3)
        self.assertEqual(shape._dimension, 2)
        self.assertIs(shape._nodes, nodes)

    def test_constructor_wrong_dimension(self):
        nodes = np.array([1.0, 2.0])
        with self.assertRaises(ValueError):
            self._make_one(nodes)

        nodes = np.zeros((2, 2, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test_constructor_insufficient_nodes(self):
        nodes = np.zeros((0, 2))
        with self.assertRaises(ValueError):
            self._make_one(nodes)

    def test__get_degree(self):
        klass = self._get_target_class()
        num_nodes = 41
        self.assertEqual(num_nodes, klass._get_degree(num_nodes))

    def test_degree_property(self):
        degree = 6
        nodes = np.zeros((degree, 2))
        shape = self._make_one(nodes)
        self.assertEqual(shape.degree, degree)
        self.assertEqual(shape._degree, degree)

    def test_dimension_property(self):
        dimension = 4
        nodes = np.zeros((3, dimension))
        shape = self._make_one(nodes)
        self.assertEqual(shape.dimension, dimension)
        self.assertEqual(shape._dimension, dimension)

    def test_nodes_property(self):
        nodes = np.array([
            [0.0, 0.0],
            [1.0, 2.0],
        ])
        shape = self._make_one(nodes)
        self.assertEqual(shape.nodes, nodes)
        self.assertIsNot(shape.nodes, nodes)

    def test_copy(self):
        np_shape = (2, 2)
        shape = self._make_one(np.zeros(np_shape))
        fake_nodes = mock.Mock(ndim=2, shape=np_shape)
        shape._nodes = fake_nodes

        copied_nodes = np.zeros(np_shape)
        fake_nodes.copy.return_value = copied_nodes

        new_shape = shape.copy()
        self.assertIsInstance(new_shape, self._get_target_class())
        self.assertIs(new_shape._nodes, copied_nodes)

        fake_nodes.copy.assert_called_once_with()

    def test___eq__(self):
        nodes1 = np.zeros((4, 3))
        shape1 = self._make_one(nodes1)
        nodes2 = np.zeros((4, 3))
        shape2 = self._make_one(nodes2)
        self.assertEqual(shape1, shape2)

    def test___ne__different_degree(self):
        nodes1 = np.zeros((1, 3))
        shape1 = self._make_one(nodes1)
        nodes2 = np.zeros((4, 3))
        shape2 = self._make_one(nodes2)
        self.assertNotEqual(shape1, shape2)

    def test___ne__different_dimension(self):
        nodes1 = np.zeros((4, 2))
        shape1 = self._make_one(nodes1)
        nodes2 = np.zeros((4, 3))
        shape2 = self._make_one(nodes2)
        self.assertNotEqual(shape1, shape2)

    def test___ne__different_class(self):
        nodes = np.zeros((4, 2))

        base_class = self._get_target_class()
        new_class1 = type('NewOne', (base_class,), {})
        new_class2 = type('NewTwo', (base_class,), {})

        shape1 = new_class1(nodes)
        shape2 = new_class2(nodes)
        self.assertNotEqual(shape1, shape2)

    def test___repr__(self):
        nodes = np.zeros((4, 3))
        shape = self._make_one(nodes)
        expected = '<Base (degree=4, dimension=3)>'
        self.assertEqual(repr(shape), expected)
