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

try:
    import importlib.metadata as importlib_metadata
except ImportError:  # pragma: NO COVER
    import importlib_metadata
import os
import pathlib
import unittest
import unittest.mock


MOCK_PACKAGE_PATH = importlib_metadata.PackagePath("extra-dll", "bezier.dll")
MOCK_PACKAGE_PATH.dist = importlib_metadata.PathDistribution(pathlib.Path(""))


class Test__get_extra_dll_dir(unittest.TestCase):
    @staticmethod
    def _call_function_under_test(bezier_files):
        from bezier import __config__

        return __config__._get_extra_dll_dir(bezier_files)

    def test_no_matches(self):
        extra_dll_dir = self._call_function_under_test(())
        self.assertIsNone(extra_dll_dir)

    def test_multiple_choices_with_match(self):
        mock_path1 = importlib_metadata.PackagePath("bezier", "__config__.py")
        mock_path1.dist = MOCK_PACKAGE_PATH.dist
        mock_path2 = importlib_metadata.PackagePath(
            "bezier", "extra-dll", "bezier.dll"
        )
        mock_path2.dist = MOCK_PACKAGE_PATH.dist
        bezier_files = (mock_path1, mock_path2)

        extra_dll_dir = self._call_function_under_test(bezier_files)
        expected = os.path.sep.join(mock_path2.parent.parts)
        self.assertEqual(extra_dll_dir, expected)


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
    @unittest.mock.patch.object(
        importlib_metadata,
        "files",
        side_effect=importlib_metadata.PackageNotFoundError,
    )
    def test_windows_without_package(self, metadata_files):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mock.
        metadata_files.assert_called_once_with("bezier")

    @unittest.mock.patch.multiple(os, name="nt", environ={})
    @unittest.mock.patch.object(
        importlib_metadata,
        "files",
        return_value=(),
    )
    def test_windows_without_dll(self, metadata_files):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mock.
        metadata_files.assert_called_once_with("bezier")

    @unittest.mock.patch.multiple(os, name="nt", environ={})
    @unittest.mock.patch("os.path.isdir", return_value=True)
    @unittest.mock.patch.object(
        importlib_metadata, "files", return_value=(MOCK_PACKAGE_PATH,)
    )
    @unittest.mock.patch("os.add_dll_directory", create=True)
    def test_windows_with_dll(
        self, os_add_dll_directory, metadata_files, isdir
    ):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mocks.
        metadata_files.assert_called_once_with("bezier")
        isdir.assert_called_once_with("extra-dll")
        os_add_dll_directory.assert_called_once_with("extra-dll")

    @unittest.mock.patch.multiple(os, name="nt", environ={})
    @unittest.mock.patch("os.path.isdir", return_value=False)
    @unittest.mock.patch.object(
        importlib_metadata, "files", return_value=(MOCK_PACKAGE_PATH,)
    )
    def test_windows_with_dll_but_not_a_dir(self, metadata_files, isdir):
        return_value = self._call_function_under_test()
        self.assertIsNone(return_value)
        self.assertEqual(os.environ, {})
        # Check mocks.
        metadata_files.assert_called_once_with("bezier")
        isdir.assert_called_once_with("extra-dll")


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
