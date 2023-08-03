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
if [[ "$(python -c 'import sys; print(".".join(map(str, sys.version_info[:2])))')" != "3.11" ]]; then
    echo "Python 3.11 required." >&2
    exit 1
fi

SCRIPT_FI=$(readlink -f "${0}")
SCRIPTS_DIR=$(dirname "${SCRIPT_FI}")
REPO_ROOT=$(dirname "${SCRIPTS_DIR}")
LOCAL_WHEELHOUSE=$(mktemp -d)
DOCKER_IMAGE=quay.io/pypa/manylinux2014_x86_64
DUMMY_CONTAINER_NAME="bezier-manylinux-$(openssl rand -hex 6)"
# Variables within the container.
PY_ROOT="/opt/python/cp311-cp311"
# NOTE: This path determines the hash of `libbezier.so` (i.e. the builds
#       are "deterministic" but not under relocation).
BEZIER_ROOT="/var/code/bezier"
WHEELHOUSE="/var/code/wheelhouse"  # Under same path as `${BEZIER_ROOT}`

# 1. Create a dummy container which will hold a volume.
docker create \
    --volume "$(dirname "${BEZIER_ROOT}")" \
    --name "${DUMMY_CONTAINER_NAME}" \
    "${DOCKER_IMAGE}" \
    /bin/true

# 2. Copy source tree into this volume.
docker cp "${REPO_ROOT}" "${DUMMY_CONTAINER_NAME}":"${BEZIER_ROOT}"

# 3. Build the `manylinux` wheel (repaired with `auditwheel`).
docker run \
    --rm \
    --env PY_ROOT="${PY_ROOT}" \
    --env WHEELHOUSE="${WHEELHOUSE}" \
    --env BEZIER_ROOT="${BEZIER_ROOT}" \
    --volumes-from="${DUMMY_CONTAINER_NAME}" \
    "${DOCKER_IMAGE}" \
    "${BEZIER_ROOT}/scripts/manylinux/build-wheel-for-doctest.sh"

# 4. Copy built wheel(s) back onto the host.
docker cp "${DUMMY_CONTAINER_NAME}":"${WHEELHOUSE}" "${LOCAL_WHEELHOUSE}"

# 5. Install the `manylinux` wheel
python -m pip install bezier \
    --no-index \
    --find-links "${LOCAL_WHEELHOUSE}/wheelhouse"

# 6. Clean up the fixed wheels and dummy container
docker container rm --volumes "${DUMMY_CONTAINER_NAME}"
rm -fr "${LOCAL_WHEELHOUSE}"
