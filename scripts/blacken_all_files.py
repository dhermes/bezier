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

"""Run `black` on all Python files."""

import pathlib
import subprocess
import sys


_SCRIPT_FILE = pathlib.Path(__file__).resolve()
ROOT_DIR = _SCRIPT_FILE.parent.parent


def main():
    all_files = subprocess.check_output(["git", "ls-files", "*.py"])
    all_files = all_files.decode("utf-8").strip()
    if not all_files:
        return

    cmd = ["black", "--line-length", "79"] + all_files.split("\n")
    status_code = subprocess.call(cmd)
    sys.exit(status_code)


if __name__ == "__main__":
    main()
