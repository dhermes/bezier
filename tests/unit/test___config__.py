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

import os
import unittest
import unittest.mock


class Test_add_dll_directory(unittest.TestCase):
    @staticmethod
    def _call_function_under_test():
        from bezier import __config__

        return __config__.add_dll_directory()

    @unittest.mock.patch.multiple(os, name="not-nt", environ={})
    @unittest.mock.patch("os.path.isdir", return_value=True)
    @unittest.mock.patch("os.add_dll_directory", create=True)
    def test_non_windows(self, os_add_dll_directory, isdir):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mocks.
        isdir.assert_not_called()
        os_add_dll_directory.assert_not_called()

    @unittest.mock.patch.multiple(os, name="nt", environ={})
    @unittest.mock.patch("os.path.isdir", return_value=True)
    @unittest.mock.patch("os.add_dll_directory", create=True)
    def test_windows_without_dll_env(self, os_add_dll_directory, isdir):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mocks.
        isdir.assert_not_called()
        os_add_dll_directory.assert_not_called()

    @unittest.mock.patch.multiple(
        os, name="nt", environ={"BEZIER_EXTRA_DLL": "builtdir\\raw"}
    )
    @unittest.mock.patch("os.path.isdir", return_value=True)
    @unittest.mock.patch("os.add_dll_directory", create=True)
    def test_windows_with_dll_env(self, os_add_dll_directory, isdir):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {"BEZIER_EXTRA_DLL": "builtdir\\raw"})
        # Check mocks.
        isdir.assert_called_once_with("builtdir\\raw")
        os_add_dll_directory.assert_called_once_with("builtdir\\raw")

    @unittest.mock.patch.multiple(
        os, name="nt", environ={"BEZIER_EXTRA_DLL": "builtdir\\raw"}
    )
    @unittest.mock.patch("os.path.isdir", return_value=False)
    @unittest.mock.patch("os.add_dll_directory", create=True)
    def test_windows_with_dll_env_but_not_a_dir(
        self, os_add_dll_directory, isdir
    ):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {"BEZIER_EXTRA_DLL": "builtdir\\raw"})
        # Check mocks.
        isdir.assert_called_once_with("builtdir\\raw")
        os_add_dll_directory.assert_not_called()


class Test_handle_import_error(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(caught_exc, name):
        from bezier import __config__

        return __config__.handle_import_error(caught_exc, name)

    def test_valid_exception(self):
        from bezier import __config__

        name = "this_module"
        single_arg = __config__.TEMPLATE.format(name)
        caught_exc = ImportError(single_arg)
        return_value = self._call_function_under_test(caught_exc, name)
        self.assertIsNone(return_value)

    def test_invalid_exception(self):
        caught_exc = ImportError("two", "args")
        with self.assertRaises(ImportError) as exc_info:
            self._call_function_under_test(caught_exc, "name")
        self.assertIs(exc_info.exception, caught_exc)
