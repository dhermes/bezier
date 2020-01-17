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

import distutils.ccompiler
import os
import shutil
import sys

import pkg_resources
import setuptools
import setuptools.command.build_ext


VERSION = "2020.1.15.dev1"  # Also in ``codemeta.json`` and ``__init__.py``.
AUTHOR = "Danny Hermes"  # Also in ``__init__.py``.
README_FILENAME = os.path.join(os.path.dirname(__file__), "README.rst")
NUMPY_MESSAGE = """\
Error: NumPy needs to be installed first. It can be installed via:

$ python    -m pip install numpy
$ python3.8 -m pip install numpy
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
INSTALL_PREFIX_ENV = "BEZIER_INSTALL_PREFIX"
NO_INSTALL_PREFIX_MESSAGE = (
    "The {} environment variable must be set."
).format(INSTALL_PREFIX_ENV)
REQUIREMENTS = ("numpy >= 1.18.1",)
EXTRAS_REQUIRE = {"full": ["scipy >= 1.4.1", "sympy >= 1.5.1"]}
DESCRIPTION = (
    u"Helper for B\u00e9zier Curves, Triangles, and Higher Order Objects"
)
_IS_WINDOWS = os.name == "nt"
_EXTRA_DLL = "extra-dll"
_DLL_FILENAME = "bezier.dll"
BAD_JOURNAL = "Saving journal failed with {!r}."
JOURNAL_ENV = "BEZIER_JOURNAL"
"""Environment variable to specify a text file for saving compiler commands.

Can be used to determine how the binary extension was compiled. This can be
useful, for example, to track changes across different systems or simply
to make sure the build is occurring as expected.
"""


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


def extension_modules():
    if NO_EXTENSION_ENV in os.environ:
        print(NO_SPEEDUPS_MESSAGE, file=sys.stderr)
        return []

    install_prefix = os.environ.get(INSTALL_PREFIX_ENV)
    if install_prefix is None:
        print(NO_INSTALL_PREFIX_MESSAGE, file=sys.stderr)
        sys.exit(1)

    rpath = os.path.join(install_prefix, "lib")
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
        return

    install_prefix = os.environ.get(INSTALL_PREFIX_ENV)
    if install_prefix is None:
        return

    # NOTE: ``bin`` is hardcoded here, expected to correspond to
    #       ``CMAKE_INSTALL_BINDIR`` on Windows.
    installed_dll = os.path.join(install_prefix, "bin", _DLL_FILENAME)
    build_lib_extra_dll = os.path.join(build_lib, "bezier", _EXTRA_DLL)
    os.makedirs(build_lib_extra_dll, exist_ok=True)
    relocated_dll = os.path.join(build_lib_extra_dll, _DLL_FILENAME)
    shutil.copyfile(installed_dll, relocated_dll)


class BuildExtWithDLL(setuptools.command.build_ext.build_ext):
    def __init__(self, *args, **kwargs):
        setuptools.command.build_ext.build_ext.__init__(self, *args, **kwargs)
        self.journal_file = os.environ.get(JOURNAL_ENV)
        self.commands = []

    def start_journaling(self):
        """Capture calls to the system by compilers.

        See: https://github.com/numpy/numpy/blob/v1.18.1/\
        numpy/distutils/ccompiler.py#L178

        Intercepts all calls to ``CCompiler.spawn`` and keeps the
        arguments around to be stored in the local ``commands``
        instance attribute.
        """
        import numpy.distutils.ccompiler

        if self.journal_file is None:
            return

        def journaled_spawn(patched_self, cmd, display=None):
            self.commands.append(cmd)
            return numpy.distutils.ccompiler.CCompiler_spawn(
                patched_self, cmd, display=None
            )

        numpy.distutils.ccompiler.replace_method(
            distutils.ccompiler.CCompiler, "spawn", journaled_spawn
        )

    @staticmethod
    def _command_to_text(command):
        # NOTE: This assumes, but doesn't check that the command has 3
        #       or more arguments.
        first_line = "$ {} \\"
        middle_line = ">   {} \\"
        last_line = ">   {}"
        parts = [first_line.format(command[0])]
        for argument in command[1:-1]:
            parts.append(middle_line.format(argument))
        parts.append(last_line.format(command[-1]))
        return "\n".join(parts)

    def _commands_to_text(self):
        separator = "-" * 40
        parts = [separator]
        for command in self.commands:
            command_text = self._command_to_text(command)
            parts.extend([command_text, separator])
        parts.append("")  # Trailing newline in file.
        return "\n".join(parts)

    def save_journal(self):
        """Save journaled commands to file.

        If there is no active journal, does nothing.

        If saving the commands to a file fails, a message will be printed to
        STDERR but the failure will be swallowed so that the extension can
        be built successfully.
        """
        if self.journal_file is None:
            return

        try:
            as_text = self._commands_to_text()
            with open(self.journal_file, "w") as file_obj:
                file_obj.write(as_text)
        except Exception as exc:
            msg = BAD_JOURNAL.format(exc)
            print(msg, file=sys.stderr)

    def run(self):
        self.start_journaling()
        copy_dll(self.build_lib)
        result = setuptools.command.build_ext.build_ext.run(self)
        self.save_journal()
        return result


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
        packages=["bezier"],
        package_dir={"": os.path.join("src", "python")},
        license="Apache 2.0",
        platforms="Posix; macOS; Windows",
        package_data={
            "bezier": [
                "*.pxd",
                os.path.join("include", "*.h"),
                os.path.join("include", "bezier", "*.h"),
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
        cmdclass={"build_ext": BuildExtWithDLL},
    )


def main():
    setup()


if __name__ == "__main__":
    main()
