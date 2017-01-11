# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Check that the current README.rst is built from the template."""


from __future__ import print_function

import os


_SCRIPTS_DIR = os.path.dirname(__file__)
_ROOT_DIR = os.path.abspath(os.path.join(_SCRIPTS_DIR, '..'))
TEMPLATE_FILE = os.path.join(_ROOT_DIR, 'README.rst.template')
ACTUAL_FILE = os.path.join(_ROOT_DIR, 'README.rst')

PYPI_IMG = """
.. |pypi| image:: https://img.shields.io/pypi/v/bezier.svg
   :target: https://pypi.python.org/pypi/bezier
   :alt: PyPI Latest"""
VERSIONS_IMG = """
.. |versions| image:: https://img.shields.io/pypi/pyversions/bezier.svg
   :target: https://pypi.python.org/pypi/bezier
   :alt: Package Versions"""


def main():
    """Populate the template and compare values.

    Raises:
        ValueError: If the current README doesn't agree with the expected
            value computed from the template.
    """
    with open(TEMPLATE_FILE, 'r') as file_obj:
        template = file_obj.read()

    expected = template.format(
        pypi='|pypi| ',
        pypi_img=PYPI_IMG,
        versions='|versions| ',
        versions_img=VERSIONS_IMG,
        rtd_version='latest',
        coveralls_branch='master',
        revision='master',
    )

    with open(ACTUAL_FILE, 'r') as file_obj:
        contents = file_obj.read()

    if contents != expected:
        raise ValueError('README.rst is not up to date with template',
                         'Expected', expected,
                         'Actual', contents)
    else:
        print('REAMDE contents are as expected.')


if __name__ == '__main__':
    main()
