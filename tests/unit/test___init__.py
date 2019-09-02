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

import email
import os
import unittest
import unittest.mock

import pkg_resources

CHECK_PKG_MSG = """\
path     = {!r}
suffix   = {!r}
site_pkg = {!r}
from_egg = {!r}"""


class Test___version__(unittest.TestCase):

    # NOTE: The ``__version__`` is hard-coded in ``__init__.py`` to
    #       accomodate builds where ``bezier`` is imported from source
    #       but not installed.

    def test_it(self):
        import bezier

        hardcoded_version = bezier.__version__
        installed_version = pkg_resources.get_distribution("bezier").version
        self.assertEqual(hardcoded_version, installed_version)


class Test___author__(unittest.TestCase):

    # NOTE: The ``__author__`` is hard-coded in ``__init__.py`` to
    #       accomodate builds where ``bezier`` is imported from source
    #       but not installed.

    def test_it(self):
        import bezier

        hardcoded_author = bezier.__author__
        distrib = pkg_resources.get_distribution("bezier")
        metadata = distrib.get_metadata(distrib.PKG_INFO)
        installed_author = email.message_from_string(metadata).get("Author")
        self.assertEqual(hardcoded_author, installed_author)


class Test_get_include(unittest.TestCase):
    @staticmethod
    def _call_function_under_test():
        import bezier

        return bezier.get_include()

    def test_it(self):
        include_directory = self._call_function_under_test()
        _check_pkg_filename(self, include_directory, "include")


class Test_get_lib(unittest.TestCase):
    @staticmethod
    def _call_function_under_test():
        import bezier

        return bezier.get_lib()

    def test_it(self):
        lib_directory = self._call_function_under_test()
        _check_pkg_filename(self, lib_directory, "lib")


class Test_get_dll(unittest.TestCase):
    @staticmethod
    def _call_function_under_test():
        import bezier

        return bezier.get_dll()

    @unittest.mock.patch("os.name", new="nt")
    def test_windows(self):
        dll_directory = self._call_function_under_test()
        _check_pkg_filename(self, dll_directory, "extra-dll")

    @unittest.mock.patch("os.name", new="posix")
    def test_non_windows(self):
        with self.assertRaises(OSError):
            self._call_function_under_test()


def _check_pkg_filename(test_case, path, last_segment):
    short = os.path.join("bezier", last_segment)
    from_egg = path.endswith(short) and ".egg" in path
    src_path = os.path.join("src", "python", short)
    from_source = path.endswith(src_path)
    verbose = os.path.join("site-packages", short)
    site_pkg = path.endswith(verbose)
    msg = CHECK_PKG_MSG.format(path, verbose, site_pkg, from_egg)
    test_case.assertTrue(from_egg or from_source or site_pkg, msg=msg)
