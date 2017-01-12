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

"""Setup file for bezier."""

from __future__ import print_function

import os
import pkg_resources
import sys


VERSION = '0.3.0'
PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(PACKAGE_ROOT, 'README.rst.template')) as file_obj:
    TEMPLATE = file_obj.read()

NUMPY_MESSAGE = """\
Error: NumPy needs to be installed first. It can be installed via:

$ pip install numpy
$ # OR
$ conda install numpy
"""
README = TEMPLATE.format(
    pypi='',
    pypi_img='',
    versions='',
    versions_img='',
    rtd_version=VERSION,
    coveralls_branch=VERSION,
    revision=VERSION,
)

REQUIREMENTS = (
    'enum34',
    'matplotlib >= 1.5.3',
    'numpy >= 1.11.2',
    'six >= 1.9.0',
)
DESCRIPTION = (
    u'Helper for B\u00e9zier Curves, Triangles, and Higher Order Objects')


def is_installed(requirement):
    try:
        pkg_resources.require(requirement)
    except pkg_resources.ResolutionError:
        return False
    else:
        return True


def extension_modules():
    # NOTE: This assumes is_installed('numpy') has already passed.
    #       H/T to http://stackoverflow.com/a/41575848/1068170
    from numpy.distutils import core

    bezier_path = os.path.join(PACKAGE_ROOT, 'bezier')
    extension = core.Extension(
        name='bezier._speedup',
        sources=[
            os.path.join(bezier_path, 'speedup.f90'),
        ],
    )
    return [extension]


def setup():
    from numpy.distutils import core

    core.setup(
        name='bezier',
        version=VERSION,
        description=DESCRIPTION,
        author='Danny Hermes',
        author_email='daniel.j.hermes@gmail.com',
        long_description=README,
        scripts=(),
        url='https://github.com/dhermes/bezier',
        packages=['bezier'],
        license='Apache 2.0',
        platforms='Posix; MacOS X; Windows',
        include_package_data=True,
        zip_safe=True,
        install_requires=REQUIREMENTS,
        ext_modules=extension_modules(),
        classifiers=(
            'Development Status :: 3 - Alpha',
            'Intended Audience :: Developers',
            'License :: OSI Approved :: Apache Software License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Topic :: Internet',
        ),
    )


def main():
    if not is_installed('numpy>=1.9.0'):
        print(NUMPY_MESSAGE, file=sys.stderr)
        sys.exit(1)

    setup()


if __name__ == '__main__':
    main()
