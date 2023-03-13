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

SCRIPT_FI=$(readlink -f ${0})
MANYLINUX_DIR=$(dirname ${SCRIPT_FI})
SCRIPTS_DIR=$(dirname ${MANYLINUX_DIR})
REPO_ROOT=$(dirname ${SCRIPTS_DIR})
DOCKER_IMAGE=quay.io/pypa/manylinux2014_i686
PRE_CMD=linux32

docker run \
    --rm \
    --volume ${REPO_ROOT}:/io \
    ${DOCKER_IMAGE} \
    ${PRE_CMD} \
    /io/scripts/manylinux/build-wheels.sh
