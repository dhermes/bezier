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

if [[ -z "${PY_ROOT}" ]]; then
    echo "PY_ROOT environment variable should be set by the caller." >&2
    exit 1
fi

if [[ -z "${WHEELHOUSE}" ]]; then
    echo "WHEELHOUSE environment variable should be set by the caller." >&2
    exit 1
fi

if [[ -z "${BEZIER_ROOT}" ]]; then
    echo "BEZIER_ROOT environment variable should be set by the caller." >&2
    exit 1
fi

# 0. Install the Python dependencies
"${PY_ROOT}/bin/python" -m pip install --upgrade pip
"${PY_ROOT}/bin/python" -m pip install --upgrade auditwheel cmake nox 'numpy >= 1.26.4, < 2'

# 1. Make sure no previous build artifacts are still around
cd "${BEZIER_ROOT}"
"${PY_ROOT}/bin/python" -m nox --session clean
rm -fr .nox/

# 2. Build and install ``libbezier``. (This script assumes it's running in an
#    ephemeral directory in an ephemeral container.)
export TARGET_NATIVE_ARCH=OFF
"${PY_ROOT}/bin/python" -m nox --session libbezier-release

# 3. Build the wheel
INSTALL_PREFIX="${BEZIER_ROOT}/.nox/.cache/libbezier-release/usr"
DIST_WHEELS=$(mktemp -d)
BEZIER_INSTALL_PREFIX="${INSTALL_PREFIX}" "${PY_ROOT}/bin/python" -m pip wheel \
    "${BEZIER_ROOT}" \
    --wheel-dir "${DIST_WHEELS}"

# 4. "repair" the built wheel.
# NOTE: This will **fail** if not run on a ``manylinux`` compatible Linux.
auditwheel repair \
    --wheel-dir "${WHEELHOUSE}" \
    "${DIST_WHEELS}"/bezier*.whl

# NOTE: We don't clean up temporary directories because this script assumes
#       it's running in an ephemeral directory in an ephemeral container.
