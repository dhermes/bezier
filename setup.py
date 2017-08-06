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
import platform
import sys

import setuptools


VERSION = '0.4.0.dev1'  # Also in codemeta.json
PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))
TEMPLATE_FILENAME = os.path.join(PACKAGE_ROOT, 'RELEASE_README.rst.template')
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
REQUIREMENTS = (
    'numpy >= 1.11.2',
    'six >= 1.9.0',
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


def has_f90_compiler():
    from distutils.ccompiler import new_compiler
    from numpy.distutils.fcompiler import new_fcompiler

    c_compiler = new_compiler()
    f90_compiler = new_fcompiler(requiref90=True, c_compiler=c_compiler)
    return f90_compiler is not None


def _extension_modules():
    # NOTE: This assumes is_installed('numpy') has already passed.
    #       H/T to https://stackoverflow.com/a/41575848/1068170
    from numpy.distutils import core

    bezier_path = os.path.join(PACKAGE_ROOT, 'src', 'bezier')
    extension = core.Extension(
        name='bezier._speedup',
        sources=[
            os.path.join(bezier_path, '_speedup.pyf'),
            os.path.join(bezier_path, 'speedup.f90'),
        ],
        language='f90',
    )
    if 'config_fc' not in sys.argv:
        sys.argv.extend(['config_fc', '--opt=-O3'])
    return [extension]


def extension_modules():
    if platform.system().lower() == 'windows':
        print(WINDOWS_MESSAGE, file=sys.stderr)
        return []
    elif has_f90_compiler():
        return _extension_modules()
    else:
        print(MISSING_F90_MESSAGE, file=sys.stderr)
        return []


def make_readme():
    with open(TEMPLATE_FILENAME, 'r') as file_obj:
        template = file_obj.read()

    return template.format(version=VERSION)


def setup():
    from numpy.distutils import core

    core.setup(
        name='bezier',
        version=VERSION,
        description=DESCRIPTION,
        author='Danny Hermes',
        author_email='daniel.j.hermes@gmail.com',
        long_description=make_readme(),
        scripts=(),
        url='https://github.com/dhermes/bezier',
        packages=setuptools.find_packages('src'),
        package_dir={'': 'src'},
        license='Apache 2.0',
        platforms='Posix; MacOS X; Windows',
        include_package_data=True,
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
    )


def main():
    if not is_installed('numpy>=1.9.0'):
        print(NUMPY_MESSAGE, file=sys.stderr)
        sys.exit(1)

    setup()


if __name__ == '__main__':
    main()
