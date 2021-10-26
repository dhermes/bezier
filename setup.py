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


VERSION = "2021.2.13.dev1"  # Also in ``codemeta.json`` and ``__init__.py``.
AUTHOR = "Danny Hermes"  # Also in ``__init__.py``.
README_FILENAME = os.path.join(os.path.dirname(__file__), "README.rst")
NUMPY_MESSAGE = """\
Error: NumPy needs to be installed first. It can be installed via:

$ python    -m pip install numpy
$ python3.9 -m pip install numpy
$ # OR
$ conda install numpy
"""
NO_EXTENSION_ENV = "BEZIER_NO_EXTENSION"
NO_SPEEDUPS_MESSAGE = """\
The {} environment variable has been used to explicitly disable the
building of the binary extension module.
""".format(
    NO_EXTENSION_ENV
)
READTHEDOCS_ENV = "READTHEDOCS"
ON_READTHEDOCS_MESSAGE = """\
The {} environment variable has been detected, the binary extension module
will not be built.
""".format(
    READTHEDOCS_ENV
)
INSTALL_PREFIX_ENV = "BEZIER_INSTALL_PREFIX"
NO_INSTALL_PREFIX_MESSAGE = (
    "The {} environment variable must be set."
).format(INSTALL_PREFIX_ENV)
WHEEL_ENV = "BEZIER_WHEEL"
"""Environment variable used to indicate a wheel is being built.

If this is present (e.g. ``BEZIER_WHEEL="True"``) then copied DLL on Windows
will be modified.
"""
# NOTE: This is a workaround, put in place for "deterministic" hashing of the
#       DLL in cases where that matters (i.e. ``doctest``.)
DLL_HASH_ENV = "BEZIER_DLL_HASH"
REQUIREMENTS = ("numpy >= 1.21.3",)
# See: https://www.python.org/dev/peps/pep-0508/
#      Dependency specification for Python Software Packages
EXTRAS_REQUIRE = {
    "full": ["matplotlib >= 3.4.3", "scipy >= 1.7.1", "sympy >= 1.9"],
}
DESCRIPTION = (
    "Helper for B\u00e9zier Curves, Triangles, and Higher Order Objects"
)
_IS_WINDOWS = os.name == "nt"
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
        print(NUMPY_MESSAGE, file=sys.stderr)
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
        print(ON_READTHEDOCS_MESSAGE, file=sys.stderr)
        return []

    if NO_EXTENSION_ENV in os.environ:
        print(NO_SPEEDUPS_MESSAGE, file=sys.stderr)
        return []

    install_prefix = os.environ.get(INSTALL_PREFIX_ENV)
    if install_prefix is None:
        print(NO_INSTALL_PREFIX_MESSAGE, file=sys.stderr)
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
            "Programming Language :: Python :: Implementation :: CPython",
            "Programming Language :: Python :: Implementation :: PyPy",
        ],
        cmdclass={"build_ext": BuildExtWithDLL},
    )


def main():
    setup()


if __name__ == "__main__":
    main()
