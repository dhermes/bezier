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

import numpy as np

from tests.unit import utils
from tests.unit.hazmat import test_geometric_intersection


@utils.needs_speedup
class Test_speedup_bbox_intersect(
    test_geometric_intersection.Test_bbox_intersect
):
    @staticmethod
    def _call_function_under_test(nodes1, nodes2):
        from bezier import _speedup

        return _speedup.bbox_intersect(nodes1, nodes2)


@utils.needs_speedup
class Test_speedup_all_intersections(
    test_geometric_intersection.Test_all_intersections
):
    @staticmethod
    def _call_function_under_test(nodes_first, nodes_second, **kwargs):
        from bezier import _speedup

        return _speedup.curve_intersections(
            nodes_first, nodes_second, **kwargs
        )

    @staticmethod
    def reset_curves_workspace(workspace_size):
        from bezier import _speedup

        return _speedup.reset_curves_workspace(workspace_size)

    @staticmethod
    def curves_workspace_size():
        from bezier import _speedup

        return _speedup.curves_workspace_size()

    def test_workspace_resize(self):
        nodes1 = np.asfortranarray([[-3.0, 5.0], [0.0, 0.0]])
        nodes2 = np.asfortranarray(
            [[-7.0, 9.0, -7.0, 9.0], [-9.0, 13.0, -13.0, 9.0]]
        )
        # NOTE: These curves intersect 3 times, so a workspace of
        #       2 is not large enough.
        self.reset_curves_workspace(2)
        intersections, coincident = self._call_function_under_test(
            nodes1, nodes2
        )
        expected = np.asfortranarray([[0.5, 0.375, 0.625], [0.5, 0.25, 0.75]])
        self.assertEqual(intersections, expected)
        self.assertFalse(coincident)
        # Make sure the workspace was resized.
        self.assertEqual(self.curves_workspace_size(), 3)

    def test_workspace_too_small(self):
        from bezier import _speedup

        nodes1 = np.asfortranarray([[-3.0, 5.0], [0.0, 0.0]])
        nodes2 = np.asfortranarray(
            [[-7.0, 9.0, -7.0, 9.0], [-9.0, 13.0, -13.0, 9.0]]
        )
        # NOTE: These curves intersect 3 times, so a workspace of
        #       2 is not large enough.
        self.reset_curves_workspace(2)
        with self.assertRaises(ValueError) as exc_info:
            self._call_function_under_test(nodes1, nodes2, allow_resize=False)
        exc_args = exc_info.exception.args
        expected = _speedup.TOO_SMALL_TEMPLATE.format(3, 2)
        self.assertEqual(exc_args, (expected,))
        # Make sure the workspace was **not** resized.
        self.assertEqual(self.curves_workspace_size(), 2)


@utils.needs_speedup
class Test_reset_curves_workspace(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(workspace_size):
        from bezier import _speedup

        return _speedup.reset_curves_workspace(workspace_size)

    def test_it(self):
        from bezier import _speedup

        size = 5
        return_value = self._call_function_under_test(size)
        self.assertIsNone(return_value)
        self.assertEqual(_speedup.curves_workspace_size(), size)

    @unittest.expectedFailure
    def test_threadsafe(self):
        from bezier import _speedup

        size_main = 3
        self._call_function_under_test(size_main)
        worker = WorkspaceThreadedAccess()
        self.assertIsNone(worker.size1)
        self.assertIsNone(worker.size2)
        size1 = 7
        size2 = 8
        thread1 = threading.Thread(target=worker.task1, args=(size1,))
        thread2 = threading.Thread(target=worker.task2, args=(size2,))
        thread1.start()
        thread2.start()
        thread1.join()
        thread2.join()
        # This check demonstrates the **broken-ness** of the implementation.
        # The sizes for each thread should be the sizes actually **set** in
        # the given thread and the workspace in the main thread should be
        # unchanged (i.e. should have ``size_main``). What we'll actually
        # observe is ``(size2, size1, size2)``.
        expected = (size1, size2, size_main)
        actual = (worker.size1, worker.size2, _speedup.curves_workspace_size())
        self.assertEqual(actual, expected)


@utils.needs_speedup
class Test_curves_workspace_size(unittest.TestCase):
    @staticmethod
    def _call_function_under_test():
        from bezier import _speedup

        return _speedup.curves_workspace_size()

    def test_it(self):
        from bezier import _speedup

        size = 5
        _speedup.reset_curves_workspace(size)
        self.assertEqual(self._call_function_under_test(), size)


class WorkspaceThreadedAccess:
    def __init__(self):
        self.barrier1 = threading.Event()
        self.barrier2 = threading.Event()
        self.barrier3 = threading.Event()
        self.size1 = None
        self.size2 = None

    def event1(self, size):
        from bezier import _speedup

        # NOTE: There is no need to ``wait`` since this is the first event.
        _speedup.reset_curves_workspace(size)
        self.barrier1.set()

    def event2(self):
        from bezier import _speedup

        self.barrier1.wait()
        result = _speedup.curves_workspace_size()
        self.barrier2.set()
        return result

    def event3(self, size):
        from bezier import _speedup

        self.barrier2.wait()
        _speedup.reset_curves_workspace(size)
        self.barrier3.set()

    def event4(self):
        from bezier import _speedup

        self.barrier3.wait()
        # NOTE: There is no barrier to ``set`` since this is the last event.
        return _speedup.curves_workspace_size()

    def task1(self, size):
        self.event1(size)
        self.size1 = self.event4()

    def task2(self, size):
        self.size2 = self.event2()
        self.event3(size)
