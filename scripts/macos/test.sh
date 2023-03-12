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

if [[ -z "${PY_VERSION}" ]]; then
    echo "PY_VERSION environment variable should be set by the caller."
    exit 1
fi

if [[ -z "${PY_BIN_DIR}" ]]; then
    echo "PY_BIN_DIR environment variable should be set by the caller."
    exit 1
fi

if [[ "${PY_VERSION}" == "3.8" ]]; then
    ${PY_BIN_DIR}/nox --session "unit-3.8"
elif [[ "${PY_VERSION}" == "3.9" ]]; then
    ${PY_BIN_DIR}/nox --session "unit-3.9"
elif [[ "${PY_VERSION}" == "3.10" ]]; then
    ${PY_BIN_DIR}/nox --session "unit-3.10"
elif [[ "${PY_VERSION}" == "3.11" ]]; then
    ${PY_BIN_DIR}/nox --session cover
    ${PY_BIN_DIR}/nox --session "functional-3.11"
    ${PY_BIN_DIR}/nox --session doctest
else
    echo "Unexpected version: ${PY_VERSION}"
    exit 1
fi
