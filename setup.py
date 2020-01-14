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

import os
import sys

import pkg_resources
import setuptools

import setup_helpers
import setup_helpers_macos
import setup_helpers_windows


VERSION = "2020.1.14"  # Also in ``codemeta.json`` and ``__init__.py``.
AUTHOR = "Danny Hermes"  # Also in ``__init__.py``.
README_FILENAME = os.path.join(os.path.dirname(__file__), "README.rst")
NUMPY_MESSAGE = """\
Error: NumPy needs to be installed first. It can be installed via:

$ python    -m pip install numpy
$ python3.8 -m pip install numpy
$ # OR
$ conda install numpy
"""
MISSING_F90_MESSAGE = """\
No Fortran 90 compiler found.

Skipping Fortran speedups via binary extension module.
"""
NO_EXTENSION_ENV = "BEZIER_NO_EXTENSION"
NO_SPEEDUPS_MESSAGE = """\
The {} environment variable has been used to explicitly disable the
building of the binary extension module.
""".format(
    NO_EXTENSION_ENV
)
REQUIREMENTS = ("numpy >= 1.18.1",)
EXTRAS_REQUIRE = {"full": ["scipy >= 1.4.1", "sympy >= 1.5.1"]}
DESCRIPTION = (
    u"Helper for B\u00e9zier Curves, Triangles, and Higher Order Objects"
)


def is_installed(requirement):
    try:
        pkg_resources.require(requirement)
    except pkg_resources.ResolutionError:
        return False

    else:
        return True


def require_numpy():
    if not is_installed("numpy>=1.9.0"):
        print(NUMPY_MESSAGE, file=sys.stderr)
        sys.exit(1)


def extension_modules():
    if NO_EXTENSION_ENV in os.environ:
        print(NO_SPEEDUPS_MESSAGE, file=sys.stderr)
        return []

    require_numpy()
    if setup_helpers.BuildFortranThenExt.has_f90_compiler():
        return setup_helpers.extension_modules()

    else:
        print(MISSING_F90_MESSAGE, file=sys.stderr)
        return []


def make_readme():
    with open(README_FILENAME, "r") as file_obj:
        return file_obj.read()


def setup():
    setuptools.setup(
        name="bezier",
        version=VERSION,
        description=DESCRIPTION,
        author=AUTHOR,
        author_email="daniel.j.hermes@gmail.com",
        long_description=make_readme(),
        scripts=(),
        url="https://github.com/dhermes/bezier",
        project_urls={
            "Documentation": "https://bezier.readthedocs.io/",
            "Changelog": (
                "https://bezier.readthedocs.io/en/2020.1.14/releases/index.html"
            ),
            "Issue Tracker": "https://github.com/dhermes/bezier/issues",
        },
        keywords=["Geometry", "Curve", "Bezier", "Intersection", "Python"],
        packages=["bezier"],
        package_dir={"": os.path.join("src", "python")},
        license="Apache 2.0",
        platforms="Posix; macOS; Windows",
        package_data={
            "bezier": [
                "*.pxd",
                os.path.join("include", "*.h"),
                os.path.join("include", "bezier", "*.h"),
                os.path.join("lib", "*.a"),
                os.path.join("lib", "*.lib"),
                os.path.join("extra-dll", "*.dll"),
            ]
        },
        zip_safe=True,
        install_requires=REQUIREMENTS,
        extras_require=EXTRAS_REQUIRE,
        ext_modules=extension_modules(),
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Mathematics",
            "License :: OSI Approved :: Apache Software License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.6",
            "Programming Language :: Python :: 3.7",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: Implementation :: CPython",
            "Programming Language :: Python :: Implementation :: PyPy",
        ],
        cmdclass={"build_ext": setup_helpers.BuildFortranThenExt},
    )


def main():
    # Add any "patches" needed for the Fortran compiler.
    setup_helpers.BuildFortranThenExt.PATCH_FUNCTIONS[:] = [
        setup_helpers.patch_f90_compiler,
        setup_helpers_macos.patch_f90_compiler,
        setup_helpers_windows.patch_f90_compiler,
    ]
    setup_helpers_windows.patch_cmd(setup_helpers.BuildFortranThenExt)
    setup()


if __name__ == "__main__":
    main()
