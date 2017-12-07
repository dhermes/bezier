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

# Avoid adding ``multibuild`` as a ``git`` submodule.
MB_COMMON_UTILS="https://raw.githubusercontent.com/matthew-brett/multibuild/master/common_utils.sh"
MB_OSX_UTILS="https://raw.githubusercontent.com/matthew-brett/multibuild/master/osx_utils.sh"

if [[ -z "${PY_VERSION}" ]]; then
    echo "PY_VERSION environment variable should be set by the caller."
    exit 1
fi

# For ``gfortran``
brew cask uninstall --force oclint
brew install gcc

# Get the "bare minimum" from ``matthew-brett/multibuild``
curl -O ${MB_COMMON_UTILS}
source common_utils.sh
curl -O ${MB_OSX_UTILS}
source osx_utils.sh

# Use ``multibuild`` to install the given python.org version of CPython.
get_macpython_environment ${PY_VERSION}
echo "PYTHON_EXE}=${PYTHON_EXE}"
echo "PIP_CMD=${PIP_CMD}"

# Make sure our installed CPython is set up for testing.
${PIP_CMD} install --ignore-installed virtualenv pip
${PIP_CMD} install --upgrade 'nox-automation >= 0.18.1' numpy

export PY_BIN_DIR=$(dirname "${PYTHON_EXE}")
echo "PY_BIN_DIR=${PY_BIN_DIR}"

# Make sure there is a universal ``libgfortran``.
${PYTHON_EXE} ${TRAVIS_BUILD_DIR}/scripts/osx/make_universal_libgfortran.py
export GFORTRAN_LIB="${TRAVIS_BUILD_DIR}/scripts/osx/frankenstein"

# Set up the tempfile directories for universal builds.
export TEMPDIR_I386=${TRAVIS_BUILD_DIR}/scripts/osx/tempdir-i386
mkdir -p ${TEMPDIR_I386}
export TEMPDIR_X86_64=${TRAVIS_BUILD_DIR}/scripts/osx/tempdir-x86_64
mkdir -p ${TEMPDIR_X86_64}
