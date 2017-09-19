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

PKG_NAME="bezier";

if [[ -z "${BIN_DIR}" ]]; then
  echo "BIN_DIR environment variable should be set by the caller."
  exit 1
fi
DELOCATE_WHEEL="${BIN_DIR}/delocate-wheel"

if [[ -z "${PY_BIN}" ]]; then
  echo "PY_BIN environment variable should be set by the caller."
  exit 1
fi

if [[ -z "${PY_TAG}" ]]; then
  echo "PY_TAG environment variable should be set by the caller."
  exit 1
fi

if [[ -f ${PY_BIN} ]]; then
    ${PY_BIN} --version
else
    echo "Python.org executable not found:"
    echo ${PY_BIN}
    exit 1
fi

# ``readlink -f`` is not our friend on OS X.
SCRIPT_FI=$(${PY_BIN} -c "import os; print(os.path.realpath('${0}'))");
OSX_DIR=$(dirname ${SCRIPT_FI});
SCRIPTS_DIR=$(dirname ${OSX_DIR});
REPO_ROOT=$(dirname ${SCRIPTS_DIR});

# Make sure packaging tools are installed / up-to-date.
${PY_BIN} -m pip install --upgrade delocate wheel virtualenv
# Make sure install requirement are installed / up-to-date.
${PY_BIN} -m pip install --upgrade setuptools numpy

# Make sure the custom ``libgfortran`` is in place.
${PY_BIN} ${OSX_DIR}/make_universal_libgfortran.py

# Create the wheel (make sure to use the custom ``libgfortran``).
DIST_WHEELS="${OSX_DIR}/dist_wheels";
mkdir -p ${DIST_WHEELS}
FRANKENSTEIN="${OSX_DIR}/frankenstein";
GFORTRAN_LIB="${FRANKENSTEIN}" \
    ${PY_BIN} -m pip wheel ${REPO_ROOT} \
    --wheel-dir ${DIST_WHEELS}

# Delocate the wheel.
FIXED_WHEELS="${OSX_DIR}/fixed_wheels";
mkdir -p ${FIXED_WHEELS}
${DELOCATE_WHEEL} \
    --check-archs \
    --wheel-dir ${FIXED_WHEELS} \
    --verbose \
    ${DIST_WHEELS}/${PKG_NAME}*${PY_TAG}*.whl

# Test out the newly created wheel in a virtual environment.
VENV="${OSX_DIR}/test-venv";
${PY_BIN} -m virtualenv ${VENV}
${VENV}/bin/pip install \
    --upgrade \
    --requirement ${OSX_DIR}/requirements.txt
${VENV}/bin/pip install \
    ${PKG_NAME} \
    --no-index \
    --find-links ${FIXED_WHEELS} \
    --find-links ${DIST_WHEELS}

set +e  # Allow tests to fail
# Run unit tests (in both 32 and 64 bit).
arch -i386 ${VENV}/bin/py.test ${REPO_ROOT}/tests/unit
arch -x86_64 ${VENV}/bin/py.test ${REPO_ROOT}/tests/unit
# Run functional tests (in both 32 and 64 bit).
PYTHONPATH=${REPO_ROOT}/functional_tests \
    arch -i386 ${VENV}/bin/py.test ${REPO_ROOT}/functional_tests
PYTHONPATH=${REPO_ROOT}/functional_tests \
    arch -x86_64 ${VENV}/bin/py.test ${REPO_ROOT}/functional_tests

# Clean-up.
rm -fr ${VENV}
rm -fr ${FRANKENSTEIN}
