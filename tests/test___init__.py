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
import sys
import unittest


class Test_get_include(unittest.TestCase):

    @staticmethod
    def _call_function_under_test():
        import bezier

        return bezier.get_include()

    def test_it(self):
        include_directory = self._call_function_under_test()
        postfix = _postfix('include')
        self.assertTrue(include_directory.endswith(postfix))


class Test_get_lib(unittest.TestCase):

    @staticmethod
    def _call_function_under_test():
        import bezier

        return bezier.get_lib()

    def test_it(self):
        lib_directory = self._call_function_under_test()
        postfix = _postfix('lib')
        self.assertTrue(lib_directory.endswith(postfix))


def _postfix(final_directory):
    python_version = 'python{}.{}'.format(*sys.version_info[:2])
    return os.path.join(
        'lib', python_version, 'site-packages', 'bezier', final_directory)
