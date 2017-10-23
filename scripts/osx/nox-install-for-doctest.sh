#!/bin/bash
#
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

set -e -x

# 0. Install the ``delocate`` tool.
pip install --upgrade delocate

# 1. Build the wheel from source.
BASIC_DIR=$(mktemp -d)
pip wheel . --wheel-dir ${BASIC_DIR}

# 2. "delocate" the built wheel.
DELOCATED_DIR=$(mktemp -d)
# NOTE: This intentionally does not use ``--check-archs``.
delocate-wheel \
    --wheel-dir ${DELOCATED_DIR} \
    --verbose \
    ${BASIC_DIR}/bezier*.whl

# 3. Install from the "delocated" wheel.
pip install ${DELOCATED_DIR}/bezier*.whl

# Clean up temporary directories.
rm -fr ${BASIC_DIR}
rm -fr ${DELOCATED_DIR}
