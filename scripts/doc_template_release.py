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

"""Populate documentation for a release.

This will introduce one-off changes in

* ``README.rst``
* ``docs/index.rst``
* ``docs/python/binary-extension.rst``
* ``DEVELOPMENT.rst``

that are not intended to be checked into ``main`` (except maybe
to be reverted after a release).

This changes will cause ``nox -s lint`` to fail because it will make
those documents look incorrect to the ``check_doc_templates.py``
script.
"""

import importlib.machinery
import os


_SCRIPTS_DIR = os.path.dirname(__file__)
_ROOT_DIR = os.path.dirname(_SCRIPTS_DIR)
README_FILE = os.path.join(_ROOT_DIR, "README.rst")
RELEASE_README_FILE = os.path.join(_ROOT_DIR, "README.rst.release.template")
INDEX_FILE = os.path.join(_ROOT_DIR, "docs", "index.rst")
RELEASE_INDEX_FILE = os.path.join(
    _ROOT_DIR, "docs", "index.rst.release.template"
)
DEVELOPMENT_TEMPLATE = os.path.join(_ROOT_DIR, "DEVELOPMENT.rst.template")
DEVELOPMENT_FILE = os.path.join(_ROOT_DIR, "DEVELOPMENT.rst")


def get_version():
    """Get the current version from ``setup.py``.

    Assumes that importing ``setup.py`` will have no side-effects (i.e.
    assumes the behavior is guarded by ``if __name__ == "__main__"``).

    Returns:
        str: The current version in ``setup.py``.
    """
    filename = os.path.join(_ROOT_DIR, "setup.py")
    loader = importlib.machinery.SourceFileLoader("setup", filename)
    setup_mod = loader.load_module()
    return setup_mod.VERSION


def populate_readme(
    version,
    linux_run,
    windows_run,
    coveralls_build,
    macos_run,
):
    """Populates ``README.rst`` with release-specific data.

    This is because ``README.rst`` is used on PyPI.

    Args:
        version (str): The current version.
        linux_run (int): The GitHub Actions run ID (for Linux) corresponding to
            the release.
        windows_run (int): The GitHub Actions run ID (for Windows)
            corresponding to the release.
        coveralls_build (Union[str, int]): The Coveralls.io build ID
            corresponding to the release.
        macos_run (int): The GitHub Actions run ID (for macOS) corresponding to
            the release.
    """
    with open(RELEASE_README_FILE, "r") as file_obj:
        template = file_obj.read()
    contents = template.format(
        version=version,
        linux_run=linux_run,
        windows_run=windows_run,
        coveralls_build=coveralls_build,
        macos_run=macos_run,
    )
    with open(README_FILE, "w") as file_obj:
        file_obj.write(contents)


def populate_index(
    version,
    linux_run,
    windows_run,
    coveralls_build,
    macos_run,
):
    """Populates ``docs/index.rst`` with release-specific data.

    Args:
        version (str): The current version.
        linux_run (int): The GitHub Actions run ID (for Linux) corresponding to
            the release.
        windows_run (int): The GitHub Actions run ID (for Windows)
            corresponding to the release.
        coveralls_build (Union[str, int]): The Coveralls.io build ID
            corresponding to the release.
        macos_run (int): The GitHub Actions run ID (for macOS) corresponding to
            the release.
    """
    with open(RELEASE_INDEX_FILE, "r") as file_obj:
        template = file_obj.read()
    contents = template.format(
        version=version,
        linux_run=linux_run,
        windows_run=windows_run,
        coveralls_build=coveralls_build,
        macos_run=macos_run,
    )
    with open(INDEX_FILE, "w") as file_obj:
        file_obj.write(contents)


def populate_development(version):
    """Populates ``DEVELOPMENT.rst`` with release-specific data.

    This is because ``DEVELOPMENT.rst`` is used in the Sphinx documentation.

    Args:
        version (str): The current version.
    """
    with open(DEVELOPMENT_TEMPLATE, "r") as file_obj:
        template = file_obj.read()
    contents = template.format(revision=version, rtd_version=version)
    with open(DEVELOPMENT_FILE, "w") as file_obj:
        file_obj.write(contents)


def main():
    """Populate the templates with release-specific fields.

    Requires user input for the GitHub Actions (Linux, macOS and Windows)
    and Coveralls.io build IDs.
    """
    version = get_version()
    linux_run = input("Linux GitHub Actions Run ID: ")
    windows_run = input("Windows GitHub Actions Run ID: ")
    coveralls_build = input("Coveralls Build ID: ")
    macos_run = input("macOS GitHub Actions Run ID: ")
    populate_readme(
        version,
        linux_run,
        windows_run,
        coveralls_build,
        macos_run,
    )
    populate_index(
        version,
        linux_run,
        windows_run,
        coveralls_build,
        macos_run,
    )
    populate_development(version)


if __name__ == "__main__":
    main()
