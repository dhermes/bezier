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


class Test_modify_path(unittest.TestCase):
    @staticmethod
    def _call_function_under_test():
        from bezier import __config__

        return __config__.modify_path()

    @unittest.mock.patch.multiple(os, name="not-nt", environ={})
    def test_non_windows(self):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})

    @unittest.mock.patch.multiple(os, name="nt", environ={})
    @unittest.mock.patch(
        "pkg_resources.resource_filename", side_effect=ImportError
    )
    def test_windows_without_dll(self, resource_filename):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mock.
        resource_filename.assert_called_once_with("bezier", "extra-dll")

    @unittest.mock.patch.multiple(
        os, name="nt", environ={"PATH": "before" + os.pathsep}
    )
    @unittest.mock.patch("os.path.isdir", return_value=True)
    @unittest.mock.patch(
        "pkg_resources.resource_filename", return_value="not-a-path"
    )
    @unittest.mock.patch("bezier.__config__.OS_ADD_DLL_DIRECTORY", new=None)
    def test_windows_with_dll(self, resource_filename, isdir):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        expected_path = "before" + os.pathsep + resource_filename.return_value
        self.assertEqual(os.environ, {"PATH": expected_path})
        # Check mocks.
        resource_filename.assert_called_once_with("bezier", "extra-dll")
        isdir.assert_called_once_with(resource_filename.return_value)

    @unittest.mock.patch.multiple(os, name="nt", environ={})
    @unittest.mock.patch("os.path.isdir", return_value=True)
    @unittest.mock.patch(
        "pkg_resources.resource_filename", return_value="not-a-path"
    )
    @unittest.mock.patch("bezier.__config__.OS_ADD_DLL_DIRECTORY", new=None)
    def test_windows_with_dll_no_path(self, resource_filename, isdir):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {"PATH": resource_filename.return_value})
        # Check mocks.
        resource_filename.assert_called_once_with("bezier", "extra-dll")
        isdir.assert_called_once_with(resource_filename.return_value)

    @unittest.mock.patch.multiple(os, name="nt", environ={})
    @unittest.mock.patch("os.path.isdir", return_value=True)
    @unittest.mock.patch(
        "pkg_resources.resource_filename", return_value="not-a-path"
    )
    @unittest.mock.patch("bezier.__config__.OS_ADD_DLL_DIRECTORY")
    def test_windows_with_dll_at_least_38(
        self, os_add_dll_directory, resource_filename, isdir
    ):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mocks.
        resource_filename.assert_called_once_with("bezier", "extra-dll")
        isdir.assert_called_once_with(resource_filename.return_value)
        os_add_dll_directory.assert_called_once_with(
            resource_filename.return_value
        )

    @unittest.mock.patch.multiple(os, name="nt", environ={})
    @unittest.mock.patch("os.path.isdir", return_value=False)
    @unittest.mock.patch(
        "pkg_resources.resource_filename", return_value="not-a-path"
    )
    def test_windows_with_dll_but_not_a_dir(self, resource_filename, isdir):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mocks.
        resource_filename.assert_called_once_with("bezier", "extra-dll")
        isdir.assert_called_once_with(resource_filename.return_value)


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
