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

GIT_ROOT=$(git rev-parse --show-toplevel)
cd "${GIT_ROOT}"

nox --session libbezier-release --reuse-existing-virtualenvs

export CIBW_BEFORE_BUILD='pip install numpy'
export CIBW_BUILD='cp39-*arm64* cp310-*arm64* cp311-*arm64*'
export CIBW_TEST_REQUIRES=pytest
export CIBW_TEST_COMMAND='pytest {project}/tests/unit'
BEZIER_INSTALL_PREFIX="$(pwd)/.nox/.cache/libbezier-release/usr"
export BEZIER_INSTALL_PREFIX

cibuildwheel --platform macos
