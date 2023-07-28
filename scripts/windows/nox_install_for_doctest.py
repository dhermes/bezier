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

"""Install and repair wheel with ``delvewheel`` on Windows.

Used in:

.. code-block:: console

   $ nox --session doctest
"""

import os
import shutil
import subprocess
import tempfile


_SCRIPTS_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.dirname(_SCRIPTS_DIR)
INSTALL_PREFIX_ENV = "BEZIER_INSTALL_PREFIX"


def main():
    # 0. Make sure we are running in the project root and install prefix
    #    environment variable is set.
    os.chdir(ROOT_DIR)
    install_prefix = os.environ.get(INSTALL_PREFIX_ENV)
    if install_prefix is None:
        raise RuntimeError("Install prefix must be set", INSTALL_PREFIX_ENV)

    # 1. Install the ``delocate`` tool.
    subprocess.call(
        ("python", "-m", "pip", "install", "--upgrade", "delvewheel")
    )

    # 2. Build the wheel from source.
    basic_dir = tempfile.TemporaryDirectory()
    # NOTE: ``pip wheel`` requires ``BEZIER_INSTALL_PREFIX`` to be set.
    subprocess.call(
        ("python", "-m", "pip", "wheel", ".", "--wheel-dir", basic_dir)
    )

    # 3. repair the built wheel.
    repaired_dir = tempfile.TemporaryDirectory()
    subprocess.call(
        (
            "python",
            "-m",
            "delvewheel",
            "repair",
            "--wheel-dir",
            repaired_dir,
            "--add-path",
            os.path.join(install_prefix, "bin"),
            os.path.join(basic_dir, "bezier*.whl"),
        )
    )

    # 4. Install from the repaired wheel.
    subprocess.call(
        (
            "python",
            "-m",
            "pip",
            "install",
            "bezier",
            "--no-index",
            "--find-links",
            repaired_dir,
        )
    )

    # 5. Clean up temporary directories.
    shutil.rmtree(basic_dir, ignore_errors=True)
    shutil.rmtree(repaired_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
