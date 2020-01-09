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

import numpy as np

from tests.unit import utils


class TestBase(utils.NumPyTestCase):
    @staticmethod
    def _get_target_class():
        from bezier import _base

        return _base.Base

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_constructor(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 3.0]])
        shape = self._make_one(nodes)
        self.assertEqual(shape._degree, -1)
        self.assertEqual(shape._dimension, 2)
        self.assertIsNot(shape._nodes, nodes)
        self.assertEqual(shape._nodes, nodes)

    def test_constructor_without_copy(self):
        nodes = np.asfortranarray([[0.0, 1.0, 2.0], [0.0, 1.0, 3.0]])
        shape = self._make_one(nodes, copy=False)
        self.assertEqual(shape._degree, -1)
        self.assertEqual(shape._dimension, 2)
        self.assertIs(shape._nodes, nodes)

    def test_constructor_wrong_dimension(self):
        nodes = np.asfortranarray([1.0, 2.0])
        with self.assertRaises(ValueError) as exc_info:
            self._make_one(nodes)
        expected = ("Nodes must be 2-dimensional, not", 1)
        self.assertEqual(exc_info.exception.args, expected)

        nodes = np.zeros((2, 2, 2), order="F")
        with self.assertRaises(ValueError) as exc_info:
            self._make_one(nodes)
        expected = ("Nodes must be 2-dimensional, not", 3)
        self.assertEqual(exc_info.exception.args, expected)

    def test_constructor_change_order(self):
        nodes = np.array([[10.0, 1.0], [3.5, 2.0]], order="C")
        shape = self._make_one(nodes, copy=False)

        self.assertFalse(nodes.flags.f_contiguous)
        self.assertTrue(shape._nodes.flags.f_contiguous)
        self.assertTrue(np.all(nodes == shape._nodes))

    def test_constructor_non_array(self):
        nodes = [[10.25, 20.0], [30.5, 4.0]]
        shape = self._make_one(nodes, copy=False)

        self.assertIsInstance(shape._nodes, np.ndarray)
        self.assertTrue(np.all(nodes == shape._nodes))

    def test_constructor_convert_dtype(self):
        nodes = np.asfortranarray([[10, 20], [30, 4]])
        shape = self._make_one(nodes, copy=False)

        self.assertEqual(nodes.dtype, np.dtype(int))
        self.assertEqual(shape._nodes.dtype, np.float64)
        self.assertTrue(np.all(nodes == shape._nodes))

    def test_constructor_rounding_failure(self):
        nodes = np.asfortranarray([[0], [73786976294838206463]])
        with self.assertRaises(ValueError) as exc_info:
            self._make_one(nodes)

        expected = ("Array cannot be converted to floating point",)
        self.assertEqual(exc_info.exception.args, expected)

    def test_degree_property(self):
        shape = self._make_one(np.zeros((1, 0), order="F"))
        self.assertEqual(shape.degree, -1)
        self.assertEqual(shape._degree, -1)

    def test_dimension_property(self):
        dimension = 4
        nodes = np.zeros((dimension, 3), order="F")
        shape = self._make_one(nodes)
        self.assertEqual(shape.dimension, dimension)
        self.assertEqual(shape._dimension, dimension)

    def test_nodes_property(self):
        nodes = np.asfortranarray([[0.0, 1.0], [0.0, 2.0]])
        shape = self._make_one(nodes)
        self.assertEqual(shape.nodes, nodes)
        self.assertIsNot(shape.nodes, nodes)

    def test___repr__(self):
        nodes = np.zeros((3, 4), order="F")
        shape = self._make_one(nodes)
        expected = "<Base (degree=-1, dimension=3)>"
        self.assertEqual(repr(shape), expected)
