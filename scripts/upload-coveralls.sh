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
#
# Upload coverage report to coveralls.io.

set -e -x

if [[ -z "${GITHUB_REF}" ]]; then
    echo "GITHUB_REF environment variable should be set by the caller." >&2
    exit 1
fi

if [[ -z "${COVERALLS_REPO_TOKEN}" ]]; then
    echo "COVERALLS_REPO_TOKEN environment variable should be set by the caller." >&2
    exit 1
fi

if [[ "${GITHUB_REF}" != "refs/heads/main" ]]; then
    echo "Coverage upload only happens on main"
    echo "Currently on ${GITHUB_REF}, doing nothing"
    exit
fi

.nox/cover/bin/python -m pip install coveralls
.nox/cover/bin/coveralls
