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

PKG_NAME="bezier"
# Install a system package required by our library
yum install -y gcc-gfortran

VERSION_WHITELIST=""
for PYBIN in /opt/python/*/bin; do
    # H/T: https://stackoverflow.com/a/229606/1068170
    if [[ "${PYBIN}" == *"36"* ]]; then
        VERSION_WHITELIST="${VERSION_WHITELIST} ${PYBIN}"
        continue
    elif [[ "${PYBIN}" == *"37"* ]]; then
        VERSION_WHITELIST="${VERSION_WHITELIST} ${PYBIN}"
        continue
    elif [[ "${PYBIN}" == *"38"* ]]; then
        VERSION_WHITELIST="${VERSION_WHITELIST} ${PYBIN}"
        continue
    else
        echo "Ignoring unsupported version: ${PYBIN}"
        echo "====================================="
    fi
done

# Compile wheels
for PYBIN in ${VERSION_WHITELIST}; do
    "${PYBIN}/pip" install -r /io/scripts/requirements.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/${PKG_NAME}*.whl; do
    auditwheel repair "${whl}" -w /io/wheelhouse/
    rm -f "${whl}"
done

# Install packages and test
for PYBIN in ${VERSION_WHITELIST}; do
    "${PYBIN}/pip" install bezier --no-index \
        --find-links /io/wheelhouse --find-links wheelhouse
    (cd "${HOME}" && "${PYBIN}/pytest" /io/tests/unit/)
    (cd "${HOME}" && "${PYBIN}/pytest" /io/tests/functional/)
done
