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

import numpy as np


class TestCurvedPolygon(unittest.TestCase):

    NODES0 = np.array([
        [0.0, 0.0],
        [0.5, -1.0],
        [1.0, 0.0],
    ])
    NODES1 = np.array([
        [1.0, 0.0],
        [0.5, 1.0],
        [0.0, 0.0],
    ])

    @staticmethod
    def _get_target_class():
        from bezier import curved_polygon

        return curved_polygon.CurvedPolygon

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def _make_default(self):
        import bezier

        edge0 = bezier.Curve(self.NODES0)
        edge1 = bezier.Curve(self.NODES1)
        return self._make_one(edge0, edge1)

    def test_constructor(self):
        import bezier

        edge0 = bezier.Curve(self.NODES0)
        edge1 = bezier.Curve(self.NODES1)
        curved_poly = self._make_one(edge0, edge1)
        self.assertEqual(curved_poly._edges, (edge0, edge1))
        self.assertEqual(curved_poly._num_sides, 2)

    def test__verify(self):
        curved_poly = self._make_one()
        self.assertIsNone(curved_poly._verify())

    def test_num_sides_property(self):
        curved_poly = self._make_default()
        self.assertIs(curved_poly.num_sides, 2)

    def test___repr__(self):
        curved_poly = self._make_default()
        self.assertEqual(repr(curved_poly),
                         '<CurvedPolygon (num_sides=2)>')
