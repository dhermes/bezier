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

# NOTE: This is assumed to be running within an activated virtual environment.
if [[ "$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')" != "3.8" ]]; then
    echo "Python 3.8 required."
    exit 1
fi

SCRIPT_FI=$(readlink -f ${0})
SCRIPTS_DIR=$(dirname ${SCRIPT_FI})
REPO_ROOT=$(dirname ${SCRIPTS_DIR})
LOCAL_WHEELHOUSE="${REPO_ROOT}/scripts/manylinux/fixed_wheels"
DOCKER_IMAGE=quay.io/pypa/manylinux2010_x86_64
DUMMY_IMAGE_NAME=bezier-manylinux
# Variables within the container.
PY_ROOT="/opt/python/cp38-cp38"
# NOTE: This path determines the hash of ``libbezier.so`` (i.e. the builds
#       are "deterministic" but not under relocation).
BEZIER_ROOT="/io"
WHEELHOUSE="/wheelhouse"

# 0. Build the `manylinux` wheel (repaired with `auditwheel`).
CURRENT_CONTAINER_ID=$(docker ps --no-trunc --filter "id=$(hostname)" --format '{{ .ID }}')
if [[ "${CI}" == "true" || "$(echo "${CURRENT_CONTAINER_ID}" | wc -w)" == "1" ]]; then
    # See: https://circleci.com/docs/2.0/building-docker-images/#mounting-folders
    # Create a dummy container which will hold a volume. This **assumes** ``/``
    # exists (i.e. the directory) containing ``BEZIER_ROOT``. If not, then
    # ``--volume "$(dirname "${BEZIER_ROOT}")"`` is required as well.
    docker create --name "${DUMMY_IMAGE_NAME}" "${DOCKER_IMAGE}" /bin/true
    # Copy source tree into this volume.
    docker cp "${REPO_ROOT}" "${DUMMY_IMAGE_NAME}":"${BEZIER_ROOT}"
    # Rely on this dummy container.
    SRC_VOLUME_ARG="--volumes-from=${DUMMY_IMAGE_NAME}"
    WHEELHOUSE_VOLUME_ARG=""
else
    # Running on host.
    SRC_VOLUME_ARG="--volume=${REPO_ROOT}:${BEZIER_ROOT}"
    WHEELHOUSE_VOLUME_ARG="--volume=${LOCAL_WHEELHOUSE}:${WHEELHOUSE}"
fi

docker run \
    --rm \
    --env PY_ROOT="${PY_ROOT}" \
    --env WHEELHOUSE="${WHEELHOUSE}" \
    --env BEZIER_ROOT="${BEZIER_ROOT}" \
    "${SRC_VOLUME_ARG}" \
    "${WHEELHOUSE_VOLUME_ARG}" \
    "${DOCKER_IMAGE}" \
    "${BEZIER_ROOT}/scripts/manylinux/build-wheel-for-doctest.sh"

if [[ "${WHEELHOUSE_VOLUME_ARG}" == "" ]]; then
    # Copy built wheel(s) back into this container.
    docker cp "${DUMMY_IMAGE_NAME}":"${WHEELHOUSE}" "${LOCAL_WHEELHOUSE}"
fi

# 1. Install the `manylinux` wheel
python -m pip install bezier \
    --no-index \
    --find-links "${LOCAL_WHEELHOUSE}"
