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

"""Remove all Cython auto-generated ``*.c`` files.

Used in:

.. code-block:: console

   $ nox --session "update_generated(check=True)"
"""

import os
import subprocess

_SCRIPTS_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.dirname(_SCRIPTS_DIR)


def main():
    # Make sure we are running in the project root.
    os.chdir(ROOT_DIR)
    c_glob = os.path.join("src", "python", "bezier", "*.c")
    all_files = subprocess.check_output(["git", "ls-files", c_glob])
    all_files = all_files.decode("utf-8").strip()
    if not all_files:
        return

    for filename in all_files.split("\n"):
        os.remove(filename)


if __name__ == "__main__":
    main()
