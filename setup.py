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

"""Setup file for bezier."""

from __future__ import print_function

import os
import pkg_resources
import sys

import setuptools

import setup_helpers
import setup_helpers_osx


VERSION = '0.5.0.dev1'  # Also in codemeta.json
README_FILENAME = os.path.join(os.path.dirname(__file__), 'README.rst')
NUMPY_MESSAGE = """\
Error: NumPy needs to be installed first. It can be installed via:

$ pip              install numpy
$ python    -m pip install numpy --user
$ python2.7 -m pip install numpy --user
$ python3.6 -m pip install numpy --user
$ # OR
$ conda install numpy
"""
MISSING_F90_MESSAGE = """\
No Fortran 90 compiler found.

Skipping Fortran extension speedups.
"""
WINDOWS_MESSAGE = """\
Skipping Fortran extension speedups on Windows.

Sorry for this inconvenience. For more information or to help, visit:

    https://github.com/dhermes/bezier/issues/26
"""
NO_EXTENSIONS_ENV = 'BEZIER_NO_EXTENSIONS'
NO_SPEEDUPS_MESSAGE = """\
The {} environment variable has been used to explicitly disable the
building of extension modules.
""".format(NO_EXTENSIONS_ENV)
REQUIREMENTS = (
    'numpy >= 1.13.3',
    'six >= 1.11.0',
)
EXTRAS_REQUIRE = {
    ':python_version<"3.4"': ['enum34'],
}
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
    if NO_EXTENSIONS_ENV in os.environ:
        print(NO_SPEEDUPS_MESSAGE)
        return []
    elif os.name == 'nt':
        print(WINDOWS_MESSAGE, file=sys.stderr)
        return []
    elif setup_helpers.BuildFortranThenExt.has_f90_compiler():
        return setup_helpers.extension_modules()
    else:
        print(MISSING_F90_MESSAGE, file=sys.stderr)
        return []


def make_readme():
    with open(README_FILENAME, 'r') as file_obj:
        return file_obj.read()


def setup():
    setuptools.setup(
        name='bezier',
        version=VERSION,
        description=DESCRIPTION,
        author='Danny Hermes',
        author_email='daniel.j.hermes@gmail.com',
        long_description=make_readme(),
        scripts=(),
        url='https://github.com/dhermes/bezier',
        packages=['bezier'],
        package_dir={'': 'src'},
        license='Apache 2.0',
        platforms='Posix; MacOS X; Windows',
        package_data={
            'bezier': [
                '*.pxd',
                os.path.join('include', '*.h'),
                os.path.join('include', 'bezier', '*.h'),
                os.path.join('lib', '*.a'),
                os.path.join('lib', '*.lib'),
            ],
        },
        zip_safe=True,
        install_requires=REQUIREMENTS,
        extras_require=EXTRAS_REQUIRE,
        ext_modules=extension_modules(),
        classifiers=(
            'Development Status :: 4 - Beta',
            'Intended Audience :: Developers',
            'Intended Audience :: Science/Research',
            'Topic :: Scientific/Engineering :: Mathematics',
            'License :: OSI Approved :: Apache Software License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 2',
            'Programming Language :: Python :: 2.7',
            'Programming Language :: Python :: 3',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ),
        cmdclass={'build_ext': setup_helpers.BuildFortranThenExt},
    )


def main():
    if not is_installed('numpy>=1.9.0'):
        print(NUMPY_MESSAGE, file=sys.stderr)
        sys.exit(1)

    # Add any "patches" needed for the Fortran compiler.
    setup_helpers.BuildFortranThenExt.PATCH_FUNCTIONS[:] = [
        setup_helpers.patch_f90_compiler,
        setup_helpers_osx.patch_f90_compiler,
    ]

    setup()


if __name__ == '__main__':
    main()
