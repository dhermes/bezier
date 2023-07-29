#!/bin/bash

set -e -x

GIT_ROOT=$(git rev-parse --show-toplevel)
cd "${GIT_ROOT}"

nox --session libbezier-release --reuse-existing-virtualenvs

export CIBW_BEFORE_BUILD='pip install numpy'
export CIBW_BUILD='cp38-*arm64* cp39-*arm64* cp310-*arm64* cp311-*arm64*'
export CIBW_TEST_REQUIRES=pytest
export CIBW_TEST_COMMAND='pytest {project}/tests/unit'
BEZIER_INSTALL_PREFIX="$(pwd)/.nox/.cache/libbezier-release/usr"
export BEZIER_INSTALL_PREFIX

cibuildwheel --platform macos
