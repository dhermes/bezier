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
import unittest.mock

import numpy as np

from tests.unit import utils


class TestCurvedPolygon(utils.NumPyTestCase):
    NODES0 = np.asfortranarray([[0.0, 0.5, 1.0], [0.0, -1.0, 0.0]])
    NODES1 = np.asfortranarray([[1.0, 0.5, 0.0], [0.0, 1.0, 0.0]])
    COLOR = (0.125, 0.125, 0.0)

    @staticmethod
    def _get_target_class():
        from bezier import curved_polygon

        return curved_polygon.CurvedPolygon

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def _make_default(self):
        import bezier

        edge0 = bezier.Curve(self.NODES0, 2)
        edge1 = bezier.Curve(self.NODES1, 2)
        return self._make_one(edge0, edge1)

    def test_constructor(self):
        import bezier

        edge0 = bezier.Curve(self.NODES0, 2)
        edge1 = bezier.Curve(self.NODES1, 2)
        curved_poly = self._make_one(edge0, edge1)
        self.assertEqual(curved_poly._edges, (edge0, edge1))
        self.assertEqual(curved_poly._num_sides, 2)
        self.assertIsNone(curved_poly._metadata)

    def test_constructor_without_verify(self):
        import bezier

        edge0 = bezier.Curve(self.NODES0, 2)
        with self.assertRaises(ValueError):
            self._make_one(edge0)
        curved_poly = self._make_one(edge0, _verify=False)
        self.assertEqual(curved_poly._edges, (edge0,))
        self.assertEqual(curved_poly._num_sides, 1)
        self.assertIsNone(curved_poly._metadata)

    def test_constructor_with_metadata(self):
        import bezier

        edge0 = bezier.Curve(self.NODES0, 2)
        edge1 = bezier.Curve(self.NODES1, 2)
        metadata = ((0, 0.0, 0.5), (4, 0.5, 1.0))
        curved_poly = self._make_one(edge0, edge1, metadata=metadata)
        self.assertEqual(curved_poly._edges, (edge0, edge1))
        self.assertEqual(curved_poly._num_sides, 2)
        self.assertEqual(curved_poly._metadata, metadata)

    def test__verify_too_few(self):
        with self.assertRaises(ValueError):
            self._make_one()
        with self.assertRaises(ValueError):
            self._make_one(None)

    def test__verify_bad_dimension(self):
        import bezier

        nodes0 = np.asfortranarray([[1.0, 2.0], [1.0, 2.0]])
        edge0 = bezier.Curve(nodes0, 1)
        edge1 = bezier.Curve(self.NODES1, 2)
        with self.assertRaises(ValueError):
            self._make_one(edge0, edge1)

    def test__verify_not_aligned(self):
        import bezier

        edge0 = bezier.Curve(np.asfortranarray([[0.0, 0.0]]), 1)
        edge1 = bezier.Curve(self.NODES1, 2)
        with self.assertRaises(ValueError):
            self._make_one(edge0, edge1)

    def test_num_sides_property(self):
        curved_poly = self._make_default()
        self.assertIs(curved_poly.num_sides, 2)

    def test___dict___property(self):
        curved_poly = self._make_default()
        props_dict = curved_poly.__dict__
        expected = {
            "_edges": curved_poly._edges,
            "_num_sides": curved_poly._num_sides,
        }
        self.assertEqual(props_dict, expected)
        # Check that modifying ``props_dict`` won't modify ``curved_poly``.
        expected["_num_sides"] = 5
        self.assertNotEqual(curved_poly._num_sides, expected["_num_sides"])

    def test_area(self):
        curved_poly = self._make_default()
        self.assertEqual(curved_poly.area, 2.0 / 3.0)

    def test___repr__(self):
        curved_poly = self._make_default()
        self.assertEqual(repr(curved_poly), "<CurvedPolygon (num_sides=2)>")

    @unittest.mock.patch("bezier._plot_helpers.new_axis")
    @unittest.mock.patch("bezier._plot_helpers.add_patch")
    def test_plot_defaults(self, add_patch_mock, new_axis_mock):
        ax = unittest.mock.Mock(spec=[])
        new_axis_mock.return_value = ax
        curved_poly = self._make_default()
        pts_per_edge = 16
        result = curved_poly.plot(pts_per_edge)
        self.assertIs(result, ax)
        # Verify mocks.
        new_axis_mock.assert_called_once_with()
        add_patch_mock.assert_called_once_with(
            ax, None, pts_per_edge, *curved_poly._edges
        )

    @unittest.mock.patch("bezier._plot_helpers.new_axis")
    @unittest.mock.patch("bezier._plot_helpers.add_patch")
    def test_plot_explicit(self, add_patch_mock, new_axis_mock):
        ax = unittest.mock.Mock(spec=[])
        color = (0.5, 0.5, 0.5)
        curved_poly = self._make_default()
        pts_per_edge = 16
        result = curved_poly.plot(pts_per_edge, color=color, ax=ax)
        self.assertIs(result, ax)
        # Verify mocks.
        new_axis_mock.assert_not_called()
        add_patch_mock.assert_called_once_with(
            ax, color, pts_per_edge, *curved_poly._edges
        )
