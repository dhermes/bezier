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

if [[ "${PY_VERSION}" == "2.7" ]]; then
    ${PY_BIN_DIR}/nox -s "unit(py='2.7')"
    arch -32 .nox/unit-py-2-7/bin/python -m pytest tests/unit/
elif [[ "${PY_VERSION}" == "3.5" ]]; then
    ${PY_BIN_DIR}/nox -s "unit(py='3.5')"
    arch -32 .nox/unit-py-3-5/bin/python -m pytest tests/unit/
elif [[ "${PY_VERSION}" == "3.6" ]]; then
    ${PY_BIN_DIR}/nox -s cover
    arch -32 .nox/cover/bin/python -m pytest tests/unit/
    ${PY_BIN_DIR}/nox -s "functional(py='3.6')"
    arch -32 .nox/functional-py-3-6/bin/python -m pytest tests/functional/
    # Apply patch to get ``nox -s doctest`` to work.
    # See: https://github.com/dhermes/bezier/issues/73
    git apply scripts/osx/issue-73.patch
    ${PY_BIN_DIR}/nox -s doctest
    ${PY_BIN_DIR}/nox -s "check_journal(machine='travis-osx')"
else
    echo "Unexpected version: ${PY_VERSION}"
    exit 1
fi
