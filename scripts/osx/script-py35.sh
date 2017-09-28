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

set -e

# NOTE: This is the Python.org version of Python.
export BIN_DIR="/Library/Frameworks/Python.framework/Versions/3.5/bin"
export PY_BIN="${BIN_DIR}/python3"
export PY_TAG="cp35-cp35m"

# ``readlink -f`` is not our friend on OS X.
SCRIPT_FI=$(${PY_BIN} -c "import os; print(os.path.realpath('${0}'))")
CURR_DIR=$(dirname ${SCRIPT_FI})
${CURR_DIR}/build-wheels.sh
