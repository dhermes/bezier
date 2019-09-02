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

"""Check if any auto-generated ``*.c`` files have changed.

Used in:

.. code-block:: console

   $ nox -s "update_generated(check=True)"
"""

from __future__ import print_function

import os
import subprocess
import sys

_SCRIPTS_DIR = os.path.abspath(os.path.dirname(__file__))
ROOT_DIR = os.path.dirname(_SCRIPTS_DIR)


def main():
    # Make sure we are running in the project root.
    os.chdir(ROOT_DIR)
    c_glob = os.path.join("src", "python", "bezier", "*.c")
    subprocess.call(["git", "add", c_glob])
    updated = subprocess.check_output(
        ["git", "diff", "HEAD", "--name-status", "--", c_glob]
    )
    if updated != b"":
        updated = updated.decode("utf-8")
        msg = "Some generated files have changed:\n{}".format(updated)
        print(msg, file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
