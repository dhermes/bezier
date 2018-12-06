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

"""Post-process a generated journal file.

This is for "quality" checks that the correct compiler flags are used
on target platforms.
"""

import argparse
import os


def post_process_travis_macos(journal_filename):
    """Post-process a generated journal file on Travis macOS.

    Args:
        journal_filename (str): The name of the journal file.
    """
    travis_build_dir = os.environ.get("TRAVIS_BUILD_DIR", "")
    with open(journal_filename, "r") as file_obj:
        content = file_obj.read()
    processed = content.replace(travis_build_dir, "${TRAVIS_BUILD_DIR}")
    with open(journal_filename, "w") as file_obj:
        file_obj.write(processed)


def post_process_journal(journal_filename, machine):
    """Post-process a generated journal file.

    Args:
        journal_filename (str): The name of the journal file.
        machine (str): The machine type where the journal was generated.
    """
    if machine == "travis-macos":
        post_process_travis_macos(journal_filename)


def main():
    parser = argparse.ArgumentParser(
        description="Post-process generated journal file."
    )
    parser.add_argument(
        "--journal-filename",
        required=True,
        help="Filename for generated journal.",
    )
    parser.add_argument(
        "--machine",
        required=True,
        help="Machine type where journal was generated.",
    )
    args = parser.parse_args()
    post_process_journal(args.journal_filename, args.machine)


if __name__ == "__main__":
    main()
