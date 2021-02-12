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
    echo "PY_ROOT environment variable should be set by the caller."
    exit 1
fi

if [[ -z "${WHEELHOUSE}" ]]; then
    echo "WHEELHOUSE environment variable should be set by the caller."
    exit 1
fi

if [[ -z "${BEZIER_ROOT}" ]]; then
    echo "BEZIER_ROOT environment variable should be set by the caller."
    exit 1
fi

# 0. Install the Python dependencies
"${PY_ROOT}/bin/python" -m pip install --upgrade pip
"${PY_ROOT}/bin/python" -m pip install --upgrade auditwheel "cmake >= 3.18.4.post1" numpy

# 1. Build and install ``libbezier`` into a custom location.
BUILD_DIR=$(mktemp -d)
INSTALL_PREFIX=$(mktemp -d)
mkdir -p "${BUILD_DIR}"
"${PY_ROOT}/bin/cmake" \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX:PATH="${INSTALL_PREFIX}" \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DTARGET_NATIVE_ARCH:BOOL=OFF \
    -S "${BEZIER_ROOT}/src/fortran/" \
    -B "${BUILD_DIR}"
"${PY_ROOT}/bin/cmake" \
    --build "${BUILD_DIR}" \
    --config Release \
    --target install
"${PY_ROOT}/bin/cmake" -L "${BUILD_DIR}"

# 2. Build the wheel
DIST_WHEELS=$(mktemp -d)
BEZIER_INSTALL_PREFIX="${INSTALL_PREFIX}" "${PY_ROOT}/bin/python" -m pip wheel \
    "${BEZIER_ROOT}" \
    --wheel-dir "${DIST_WHEELS}"

# 3. "repair" the built wheel.
# NOTE: This will **fail** if not run on a ``manylinux`` compatible Linux.
auditwheel repair \
    --wheel-dir "${WHEELHOUSE}" \
    "${DIST_WHEELS}"/bezier*.whl

# 4. Clean up temporary directories.
rm -fr "${BUILD_DIR}"
rm -fr "${INSTALL_PREFIX}"
rm -fr "${DIST_WHEELS}"
