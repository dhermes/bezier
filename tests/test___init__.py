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

import distutils.sysconfig
import os
import platform
import unittest


PLATFORM_SYSTEM = platform.system().lower()
CHECK_PKG_MSG = """\
path     = {!r}
expected = {!r}
site_pkg = {!r}
from_egg = {!r}"""


class Test_get_include(unittest.TestCase):

    @staticmethod
    def _call_function_under_test():
        import bezier

        return bezier.get_include()

    def test_it(self):
        include_directory = self._call_function_under_test()
        _check_pkg_filename(self, include_directory, 'include')


class Test_get_lib(unittest.TestCase):

    @staticmethod
    def _call_function_under_test():
        import bezier

        return bezier.get_lib()

    @unittest.skipIf(
        PLATFORM_SYSTEM == 'windows',
        'Static library not yet built on Windows')
    def test_it(self):
        lib_directory = self._call_function_under_test()
        _check_pkg_filename(self, lib_directory, 'lib')


def _check_pkg_filename(test_case, path, last_segment):
    suffix = os.path.join('bezier', last_segment)
    from_egg = path.endswith(suffix) and '.egg' in path

    site_packages = distutils.sysconfig.get_python_lib()
    expected = os.path.join(site_packages, suffix)
    site_pkg = path == expected

    msg = CHECK_PKG_MSG.format(
        path, expected, site_pkg, from_egg)
    test_case.assertTrue(site_pkg or from_egg, msg=msg)
