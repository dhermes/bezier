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

# Avoid adding ``multibuild`` as a ``git`` submodule.
MB_COMMON_UTILS="https://raw.githubusercontent.com/matthew-brett/multibuild/devel/common_utils.sh"
MB_OSX_UTILS="https://raw.githubusercontent.com/matthew-brett/multibuild/devel/osx_utils.sh"

if [[ -z "${PY_VERSION}" ]]; then
    echo "PY_VERSION environment variable should be set by the caller."
    exit 1
fi

# Get the "bare minimum" from ``matthew-brett/multibuild``
curl -O ${MB_COMMON_UTILS}
source common_utils.sh
curl -O ${MB_OSX_UTILS}
source osx_utils.sh

# Use ``multibuild`` to install the given python.org version of CPython.
get_macpython_environment ${PY_VERSION}

echo "PYTHON_EXE=${PYTHON_EXE}"
# NOTE: This assumes that ``multibuild`` uses ``sudo .../bin/python/pipX.Y``
#       as the command (i.e. it's missing the ``-H`` flag).
PIP_CMD=${PIP_CMD/sudo/sudo -H}
echo "PIP_CMD=${PIP_CMD}"

# Make sure our installed CPython is set up for testing.
${PIP_CMD} install --ignore-installed virtualenv pip
${PIP_CMD} install --upgrade "nox >= 2021.10.1" numpy
# Make sure there is a working ``cmake`` on the ``${PATH}``.
${PIP_CMD} install --upgrade "cmake >= 3.21.4"
command -v cmake

PY_BIN_DIR=$(dirname "${PYTHON_EXE}")
echo "PY_BIN_DIR=${PY_BIN_DIR}"
