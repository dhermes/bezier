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
import unittest

import pkg_resources


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
