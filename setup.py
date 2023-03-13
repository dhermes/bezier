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

"""Setup file for ``bezier``."""

import hashlib
import os
import pathlib
import shutil
import sys

import pkg_resources
import setuptools
import setuptools.command.build_ext
import setuptools.dist


VERSION = "2021.2.13.dev1"  # Also in ``codemeta.json`` and ``__init__.py``.
AUTHOR = "Danny Hermes"  # Also in ``__init__.py``.
README_FILENAME = os.path.join(os.path.dirname(__file__), "README.rst")
NUMPY_MESSAGE = """\
Error: NumPy needs to be installed first. It can be installed via:

$ python     -m pip install numpy
$ python3.11 -m pip install numpy
$ # OR
$ conda install numpy
"""
NO_EXTENSION_ENV = "BEZIER_NO_EXTENSION"
NO_SPEEDUPS_MESSAGE = """\
The `{}` environment variable has been used to explicitly disable the
building of the binary extension module.
""".format(
    NO_EXTENSION_ENV
)
IGNORE_VERSION_CHECK_ENV = "BEZIER_IGNORE_VERSION_CHECK"
INVALID_VERSION_MESSAGE = """\
The current Python version ({major}.{minor}) is not supported.

The supported versions are: {versions}

Using `bezier` on an unsupported version of Python is not known to work. You
may be seeing this message as part of a source distribution (`sdist`) install
because no wheels exist on PyPI for your current Python environment. The
Python environment is uniquely identified by Python version, operating
system (ABI) and architecture (platform). You are likely seeing this message
because a new version of Python has been released. To disable this check, set
the `BEZIER_IGNORE_VERSION_CHECK` environment variable.
"""
READTHEDOCS_ENV = "READTHEDOCS"
ON_READTHEDOCS_MESSAGE = """\
The `{}` environment variable has been detected, the binary extension module
will not be built.
""".format(
    READTHEDOCS_ENV
)
INSTALL_PREFIX_ENV = "BEZIER_INSTALL_PREFIX"
NO_INSTALL_PREFIX_MESSAGE = """\
The `{install_prefix}` environment variable must be set when installing
`bezier` from source with the binary extension module. If you are not
intending to install from source, check on PyPI to see if a prebuilt wheel
exists for your current Python environment. If no wheel exists, installing
from source is the only option (even via `pip`, which will use the `sdist`
source distribution).

For a pure Python install (i.e. without the binary extension module), set
the `{no_extension}` environment variable.

For a source install using the binary extension module, see installation
instructions for the C ABI (`libbezier`).
""".format(
    install_prefix=INSTALL_PREFIX_ENV, no_extension=NO_EXTENSION_ENV
)
WHEEL_ENV = "BEZIER_WHEEL"
"""Environment variable used to indicate a wheel is being built.

If this is present (e.g. ``BEZIER_WHEEL="True"``) then copied DLL on Windows
will be modified.
"""
# NOTE: This is a workaround, put in place for "deterministic" hashing of the
#       DLL in cases where that matters (i.e. ``doctest``.)
DLL_HASH_ENV = "BEZIER_DLL_HASH"
REQUIREMENTS = ("numpy >= 1.24.2",)
# See: https://www.python.org/dev/peps/pep-0508/
#      Dependency specification for Python Software Packages
EXTRAS_REQUIRE = {
    "full": ["matplotlib >= 3.4.3", "scipy >= 1.10.1", "sympy >= 1.11.1"],
}
DESCRIPTION = (
    "Helper for B\u00e9zier Curves, Triangles, and Higher Order Objects"
)
_IS_WINDOWS = os.name == "nt"
_IS_PYPY = sys.implementation.name == "pypy"
_EXTRA_DLL = "extra-dll"
_DLL_FILENAME = "bezier.dll"


def is_installed(requirement):
    try:
        pkg_resources.require(requirement)
    except pkg_resources.ResolutionError:
        return False

    else:
        return True


def numpy_include_dir():
    if not is_installed("numpy >= 1.9.0"):
        print(NUMPY_MESSAGE, file=sys.stderr, end="")
        sys.exit(1)

    import numpy as np

    return np.get_include()


def _sha256_hash(filename, blocksize=65536):
    """Hash the contents of an open file handle with SHA256"""
    hash_obj = hashlib.sha256()

    with open(filename, "rb") as file_obj:
        block = file_obj.read(blocksize)
        while block:
            hash_obj.update(block)
            block = file_obj.read(blocksize)

    return hash_obj.hexdigest()


def _sha256_short_hash(filename):
    full_hash = _sha256_hash(filename)
    return full_hash[:8]


def extension_modules():
    if os.environ.get(READTHEDOCS_ENV) == "True":
        print(ON_READTHEDOCS_MESSAGE, file=sys.stderr, end="")
        return []

    if NO_EXTENSION_ENV in os.environ:
        print(NO_SPEEDUPS_MESSAGE, file=sys.stderr, end="")
        return []

    install_prefix = os.environ.get(INSTALL_PREFIX_ENV)
    if install_prefix is None:
        print(NO_INSTALL_PREFIX_MESSAGE, file=sys.stderr, end="")
        sys.exit(1)

    rpath = os.path.join(install_prefix, "lib")
    if not os.path.isdir(rpath):
        rpath = os.path.join(install_prefix, "lib64")
    extra_link_args = []
    if not _IS_WINDOWS:
        extra_link_args.append("-Wl,-rpath,{}".format(rpath))

    extension = setuptools.Extension(
        "bezier._speedup",
        [os.path.join("src", "python", "bezier", "_speedup.c")],
        include_dirs=[
            numpy_include_dir(),
            os.path.join(install_prefix, "include"),
        ],
        libraries=["bezier"],
        library_dirs=[rpath],
        extra_link_args=extra_link_args,
    )
    return [extension]


def make_readme():
    with open(README_FILENAME, "r") as file_obj:
        return file_obj.read()


def copy_dll(build_lib):
    if not _IS_WINDOWS:
        return None

    install_prefix = os.environ.get(INSTALL_PREFIX_ENV)
    if install_prefix is None:
        return None

    # NOTE: ``bin`` is hardcoded here, expected to correspond to
    #       ``CMAKE_INSTALL_BINDIR`` on Windows.
    installed_dll = os.path.join(install_prefix, "bin", _DLL_FILENAME)
    build_lib_extra_dll = os.path.join(build_lib, "bezier", _EXTRA_DLL)
    os.makedirs(build_lib_extra_dll, exist_ok=True)

    if WHEEL_ENV in os.environ:
        short_hash = os.environ.get(DLL_HASH_ENV)
        if short_hash is None:
            short_hash = _sha256_short_hash(installed_dll)
        dll_name = f"bezier-{short_hash}.dll"
        return_value = dll_name
    else:
        dll_name = _DLL_FILENAME
        return_value = None

    relocated_dll = os.path.join(build_lib_extra_dll, dll_name)
    shutil.copyfile(installed_dll, relocated_dll)

    return return_value


def _find_speedup_pyd(build_lib):
    bezier_dir = pathlib.Path(build_lib) / "bezier"
    speedup_glob = "_speedup*.pyd"
    matches = list(bezier_dir.glob(speedup_glob))
    if len(matches) != 1:
        raise ValueError(f"Could not find unique ``{speedup_glob}``", matches)

    return str(matches[0])


def rewrite_pyd_reference(dll_name, build_lib):
    if dll_name is None:
        return

    import machomachomangler.pe

    speedup_filename = _find_speedup_pyd(build_lib)
    with open(speedup_filename, "rb") as file_obj:
        pyd_bytes = file_obj.read()

    new_pyd_bytes = machomachomangler.pe.redll(
        pyd_bytes, {_DLL_FILENAME.encode("ascii"): dll_name.encode("ascii")}
    )
    with open(speedup_filename, "wb") as file_obj:
        file_obj.write(new_pyd_bytes)


class BuildExtWithDLL(setuptools.command.build_ext.build_ext):
    def run(self):
        dll_name = copy_dll(self.build_lib)
        setuptools.command.build_ext.build_ext.run(self)
        rewrite_pyd_reference(dll_name, self.build_lib)


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
                "https://bezier.readthedocs.io/en/latest/releases/index.html"
            ),
            "Issue Tracker": "https://github.com/dhermes/bezier/issues",
        },
        keywords=["Geometry", "Curve", "Bezier", "Intersection", "Python"],
        packages=["bezier", "bezier.hazmat"],
        package_dir={"": os.path.join("src", "python")},
        license="Apache 2.0",
        platforms="Posix; macOS; Windows",
        package_data={
            "bezier": [
                "*.pxd",
                os.path.join("extra-dll", "*.dll"),
            ]
        },
        zip_safe=True,
        install_requires=REQUIREMENTS,
        extras_require=EXTRAS_REQUIRE,
        ext_modules=extension_modules(),
        python_requires=">=3.8",
        classifiers=[
            "Development Status :: 4 - Beta",
            "Intended Audience :: Developers",
            "Intended Audience :: Science/Research",
            "Topic :: Scientific/Engineering :: Mathematics",
            "License :: OSI Approved :: Apache Software License",
            "Operating System :: OS Independent",
            "Programming Language :: Python :: 3",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: 3.11",
            "Programming Language :: Python :: Implementation :: CPython",
            "Programming Language :: Python :: Implementation :: PyPy",
        ],
        cmdclass={"build_ext": BuildExtWithDLL},
    )


def _check_python_version():
    """Check that this is being installed in a valid version of Python.

    If the current version of Python is unsupported, this will exit with an
    error.

    The ``BEZIER_IGNORE_VERSION_CHECK`` environment variable can be set to
    opt out of this check.
    """
    if IGNORE_VERSION_CHECK_ENV in os.environ:
        return

    major = sys.version_info.major
    minor = sys.version_info.minor
    if (major, minor) in ((3, 8), (3, 9), (3, 10), (3, 11)):
        return

    message = INVALID_VERSION_MESSAGE.format(
        major=major, minor=minor, versions="3.8, 3.9, 3.10 and 3.11"
    )
    print(message, file=sys.stderr, end="")
    sys.exit(1)


def _patch_setuptools():
    """Patch ``setuptools`` to address known issues.

    Known issues:
    * In some PyPy installs, the ``setuptools.build_py.build_package_data()``
      method depends on the ``convert_2to3_doctests`` being set on an instance
      of ``setuptools.dist.Distribution``, but it is unset. We handle this
      by setting it as a **class attribute** (vs. monkey-patching ``__init__``
      to set it on instances).
    """
    if not _IS_PYPY:
        return

    setuptools.dist.Distribution.convert_2to3_doctests = []


def main():
    _check_python_version()
    _patch_setuptools()
    setup()


if __name__ == "__main__":
    main()
