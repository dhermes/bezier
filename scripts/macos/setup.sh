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

if [[ -z "${Python_ROOT_DIR}" ]]; then
    echo "Python_ROOT_DIR environment variable should be set by the caller." >&2
    exit 1
fi

PIP_CMD="${Python_ROOT_DIR}/bin/python -m pip"
echo "PIP_CMD=${PIP_CMD}"

# Make sure our installed CPython is set up for testing.
${PIP_CMD} install --ignore-installed virtualenv pip
${PIP_CMD} install --upgrade "nox >= 2022.11.21" numpy
# Make sure there is a working ``cmake`` on the ``${PATH}``.
${PIP_CMD} install --upgrade "cmake >= 3.25.2"
command -v cmake
