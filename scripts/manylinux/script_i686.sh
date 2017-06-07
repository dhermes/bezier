#!/bin/bash

SCRIPT_FI=$(readlink -f ${0});
MANYLINUX_DIR=$(dirname ${SCRIPT_FI});
SCRIPTS_DIR=$(dirname ${MANYLINUX_DIR});
REPO_ROOT=$(dirname ${SCRIPTS_DIR});
DOCKER_IMAGE=quay.io/pypa/manylinux1_i686
PRE_CMD=linux32

docker run --rm -v ${REPO_ROOT}:/io ${DOCKER_IMAGE} ${PRE_CMD} /io/scripts/manylinux/build-wheels.sh
