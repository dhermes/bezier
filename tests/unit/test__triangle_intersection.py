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

import threading
import unittest

from tests import utils as base_utils
from tests.unit import test__py_triangle_intersection
from tests.unit import utils


@utils.needs_speedup
class Test_speedup_newton_refine(
    test__py_triangle_intersection.Test_newton_refine
):
    @staticmethod
    def _call_function_under_test(nodes, degree, x_val, y_val, s, t):
        from bezier import _speedup

        return _speedup.newton_refine_triangle(
            nodes, degree, x_val, y_val, s, t
        )


@utils.needs_speedup
class Test_speedup_locate_point(
    test__py_triangle_intersection.Test_locate_point
):
    @staticmethod
    def _call_function_under_test(nodes, degree, x_val, y_val):
        from bezier import _speedup

        return _speedup.locate_point_triangle(nodes, degree, x_val, y_val)


@utils.needs_speedup
class Test_speedup_geometric_intersect(
    test__py_triangle_intersection.Test_geometric_intersect
):
    BAD_BOUNDARY_ARGS = ("Unexpected number of edges",)
    BAD_BOUNDARY_TYPE = RuntimeError
    BAD_BOUNDARY_INCREASE_ULPS = True

    @staticmethod
    def _call_function_under_test(nodes1, degree1, nodes2, degree2, **kwargs):
        from bezier import _speedup

        return _speedup.triangle_intersections(
            nodes1, degree1, nodes2, degree2, **kwargs
        )

    def test_two_curved_polygons(self):
        # Make sure there is enough space so that no resize is needed.
        sizes = triangle_workspace_sizes()
        segment_ends_size, segments_size = sizes
        self.assertGreaterEqual(segment_ends_size, 2)
        self.assertGreaterEqual(segments_size, 6)
        super_ = super(Test_speedup_geometric_intersect, self)
        super_.test_two_curved_polygons()
        # Make sure the workspace was **not** resized.
        self.assertEqual(triangle_workspace_sizes(), sizes)

    def test_resize_both(self):
        reset_triangle_workspaces(segment_ends_size=1, segments_size=1)
        super_ = super(Test_speedup_geometric_intersect, self)
        super_.test_two_curved_polygons()
        # Make sure the sizes were resized from (1, 1).
        self.assertEqual(triangle_workspace_sizes(), (2, 6))

    def test_insufficient_segment_ends(self):
        from bezier import _speedup

        reset_triangle_workspaces(segment_ends_size=1)
        sizes = triangle_workspace_sizes()
        with self.assertRaises(ValueError) as exc_info:
            self._two_curved_polygons(resizes_allowed=0)
        exc_args = exc_info.exception.args
        template = _speedup.SEGMENT_ENDS_TOO_SMALL
        self.assertEqual(exc_args, (template.format(2, 1),))
        # Make sure the workspace was **not** resized.
        self.assertEqual(triangle_workspace_sizes(), sizes)

    def test_insufficient_segments(self):
        from bezier import _speedup

        reset_triangle_workspaces(segment_ends_size=2, segments_size=2)
        sizes = triangle_workspace_sizes()
        with self.assertRaises(ValueError) as exc_info:
            self._two_curved_polygons(resizes_allowed=0)
        exc_args = exc_info.exception.args
        template = _speedup.SEGMENTS_TOO_SMALL
        self.assertEqual(exc_args, (template.format(6, 2),))
        # Make sure the workspace was **not** resized.
        self.assertEqual(triangle_workspace_sizes(), sizes)


@utils.needs_speedup
class Test_speedup__type_info(unittest.TestCase):
    @staticmethod
    def _call_function_under_test():
        from bezier import _speedup

        return _speedup._type_info()

    def test_it(self):
        result = self._call_function_under_test()
        is_native, item_size, dtype_num, size_of_struct = result
        self.assertTrue(is_native)
        self.assertEqual(dtype_num, 20)
        if base_utils.IS_64_BIT or base_utils.IS_WINDOWS:
            self.assertEqual(item_size, 24)
            self.assertEqual(size_of_struct, 24)
        else:  # pragma: NO COVER
            self.assertEqual(item_size, 20)
            self.assertEqual(size_of_struct, 20)


@utils.needs_speedup
class Test_reset_triangle_workspaces(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(**kwargs):
        return reset_triangle_workspaces(**kwargs)

    def test_it(self):
        return_value = self._call_function_under_test(
            segment_ends_size=1, segments_size=2
        )
        self.assertIsNone(return_value)
        self.assertEqual(triangle_workspace_sizes(), (1, 2))

    @unittest.expectedFailure
    def test_threadsafe(self):
        sizes_main = (4, 3)
        self._call_function_under_test(
            segment_ends_size=sizes_main[0], segments_size=sizes_main[1]
        )
        worker = WorkspaceThreadedAccess()
        self.assertIsNone(worker.sizes1)
        self.assertIsNone(worker.sizes2)
        sizes1 = (1, 3)
        sizes2 = (2, 2)
        thread1 = threading.Thread(target=worker.task1, args=(sizes1,))
        thread2 = threading.Thread(target=worker.task2, args=(sizes2,))
        thread1.start()
        thread2.start()
        thread1.join()
        thread2.join()
        # This check demonstrates the **broken-ness** of the implementation.
        # The sizes for each thread should be the sizes actually **set** in
        # the given thread and the workspace in the main thread should be
        # unchanged (i.e. should have ``sizes_main``). What we'll actually
        # observe is ``(sizes2, sizes1, sizes2)``.
        expected = (sizes1, sizes2, sizes_main)
        actual = (worker.sizes1, worker.sizes2, triangle_workspace_sizes())
        self.assertEqual(actual, expected)


@utils.needs_speedup
class Test_triangle_workspace_sizes(unittest.TestCase):
    @staticmethod
    def _call_function_under_test():
        return triangle_workspace_sizes()

    def test_it(self):
        reset_triangle_workspaces(segment_ends_size=3, segments_size=5)
        self.assertEqual(self._call_function_under_test(), (3, 5))
        reset_triangle_workspaces(segment_ends_size=1)
        self.assertEqual(self._call_function_under_test(), (1, 5))
        reset_triangle_workspaces(segments_size=2)
        self.assertEqual(self._call_function_under_test(), (1, 2))


def reset_triangle_workspaces(**kwargs):
    from bezier import _speedup

    return _speedup.reset_triangle_workspaces(**kwargs)


def triangle_workspace_sizes():
    from bezier import _speedup

    return _speedup.triangle_workspace_sizes()


class WorkspaceThreadedAccess:
    def __init__(self):
        self.barrier1 = threading.Event()
        self.barrier2 = threading.Event()
        self.barrier3 = threading.Event()
        self.sizes1 = None
        self.sizes2 = None

    def event1(self, sizes):
        # NOTE: There is no need to ``wait`` since this is the first event.
        reset_triangle_workspaces(
            segment_ends_size=sizes[0], segments_size=sizes[1]
        )
        self.barrier1.set()

    def event2(self):
        self.barrier1.wait()
        result = triangle_workspace_sizes()
        self.barrier2.set()
        return result

    def event3(self, sizes):
        self.barrier2.wait()
        reset_triangle_workspaces(
            segment_ends_size=sizes[0], segments_size=sizes[1]
        )
        self.barrier3.set()

    def event4(self):
        self.barrier3.wait()
        # NOTE: There is no barrier to ``set`` since this is the last event.
        return triangle_workspace_sizes()

    def task1(self, sizes):
        self.event1(sizes)
        self.sizes1 = self.event4()

    def task2(self, sizes):
        self.sizes2 = self.event2()
        self.event3(sizes)
