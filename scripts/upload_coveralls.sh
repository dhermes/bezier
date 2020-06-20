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
# Build the bezier docs.

set -e -x

if [[ -z "${CIRCLE_BRANCH}" ]]; then
    echo "CIRCLE_BRANCH environment variable should be set by the caller."
    exit 1
fi

if [[ "${CIRCLE_BRANCH}" != "main" ]]; then
    echo "Coverage upload only happens on main"
    echo "Currently on ${CIRCLE_BRANCH}, doing nothing"
    exit
fi

.nox/cover/bin/python -m pip install coveralls
.nox/cover/bin/coveralls
