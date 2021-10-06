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

# Some helpful links:
# - https://docs.docker.com/engine/installation/linux/ubuntu/
# - https://github.com/pypa/python-manylinux-demo/blob/master/.travis.yml
# - https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh

# NOTE: On the `manylinux2010` images (which `rpm -q centos-release` indicates
#       are CentOS 6, release 10.el6.centos.12.3) there is a pre-installed
#       `gfortran` on the path under `/opt/rh`. Running
#       `yum install gcc-gfortran` installs an **older** `/usr/bin/gfortran`
#       so we don't do that.

# Install (new) CMake into Python 3.10 environment.
/opt/python/cp310-cp310/bin/python -m pip install --upgrade pip
/opt/python/cp310-cp310/bin/python -m pip install "cmake >= 3.21.3"

# Build and install ``libbezier`` into a custom location.
SRC_DIR="/io/src/fortran/"
BUILD_DIR="${HOME}/libbezier-release/build"
INSTALL_PREFIX="${HOME}/libbezier-release/usr"
mkdir -p "${BUILD_DIR}"
/opt/python/cp310-cp310/bin/cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX:PATH="${INSTALL_PREFIX}" \
    -DCMAKE_VERBOSE_MAKEFILE:BOOL=ON \
    -DTARGET_NATIVE_ARCH:BOOL=OFF \
    -S "${SRC_DIR}" \
    -B "${BUILD_DIR}"
/opt/python/cp310-cp310/bin/cmake \
    --build "${BUILD_DIR}" \
    --config Release \
    --target install
/opt/python/cp310-cp310/bin/cmake -L "${BUILD_DIR}"

VERSION_WHITELIST=""
for PYBIN in /opt/python/*/bin; do
    # H/T: https://stackoverflow.com/a/229606/1068170
    if [[ "${PYBIN}" == *"37"* ]]; then
        VERSION_WHITELIST="${VERSION_WHITELIST} ${PYBIN}"
        continue
    elif [[ "${PYBIN}" == *"38"* ]]; then
        VERSION_WHITELIST="${VERSION_WHITELIST} ${PYBIN}"
        continue
    elif [[ "${PYBIN}" == *"39"* ]]; then
        VERSION_WHITELIST="${VERSION_WHITELIST} ${PYBIN}"
        continue
    elif [[ "${PYBIN}" == *"310"* ]]; then
        VERSION_WHITELIST="${VERSION_WHITELIST} ${PYBIN}"
        continue
    else
        echo "Ignoring unsupported version: ${PYBIN}"
        echo "====================================="
    fi
done

# Compile wheels
for PYBIN in ${VERSION_WHITELIST}; do
    "${PYBIN}/python" -m pip install --upgrade pip
    "${PYBIN}/python" -m pip install --requirement /io/scripts/requirements.txt
    BEZIER_INSTALL_PREFIX="${INSTALL_PREFIX}" \
        "${PYBIN}/python" -m pip \
        wheel \
        /io/ \
        --wheel-dir "${HOME}/wheelhouse/"
done

# Bundle external shared libraries into the wheels
for whl in "${HOME}/wheelhouse/bezier*.whl"; do
    auditwheel repair "${whl}" --wheel-dir /io/wheelhouse/
    rm -f "${whl}"
done

# Install packages and test
for PYBIN in ${VERSION_WHITELIST}; do
    "${PYBIN}/python" -m pip install bezier \
        --no-index \
        --find-links /io/wheelhouse
    (cd "${HOME}" && "${PYBIN}/python" -m pytest /io/tests/unit/)
    (cd "${HOME}" && "${PYBIN}/python" -m pytest /io/tests/functional/)
done
