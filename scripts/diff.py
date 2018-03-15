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
"""Compare two files.

Intended as a platform-independent version of ``diff``.
"""

from __future__ import print_function

import argparse
import difflib
import os
import sys

NOT_EXISTS = '{} does not exist.'


def diff(filename1, filename2):
    if not os.path.exists(filename1):
        msg = NOT_EXISTS.format(filename1)
        print(msg, file=sys.stderr)
        sys.exit(1)
    if not os.path.exists(filename2):
        msg = NOT_EXISTS.format(filename2)
        print(msg, file=sys.stderr)
        sys.exit(1)
    with open(filename1, 'r') as file_obj:
        lines1 = file_obj.readlines()
    with open(filename2, 'r') as file_obj:
        lines2 = file_obj.readlines()
    if lines1 == lines2:
        msg = '{} and {} are identical'.format(filename1, filename2)
        print(msg)
    else:
        diff_lines = difflib.context_diff(
            lines1, lines2, fromfile=filename1, tofile=filename2
        )
        print(''.join(diff_lines))
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(
        description='Print diff between two files.'
    )
    parser.add_argument('filename1', help='First file to compare.')
    parser.add_argument('filename2', help='Second file to compare.')
    args = parser.parse_args()
    diff(args.filename1, args.filename2)


if __name__ == '__main__':
    main()
