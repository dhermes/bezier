#!/bin/bash
set -e -x

# Some helpful links:
# - https://docs.docker.com/engine/installation/linux/ubuntu/
# - https://github.com/pypa/python-manylinux-demo/blob/master/.travis.yml
# - https://github.com/pypa/python-manylinux-demo/blob/master/travis/build-wheels.sh

PKG_NAME="bezier";
# Install a system package required by our library
yum install -y gcc-gfortran

# Compile wheels
for PYBIN in /opt/python/*/bin; do
    # H/T: https://stackoverflow.com/a/229606/1068170
    if [[ "${PYBIN}" == *"26"* ]]; then
        echo "Python 2.6 (via ${PYBIN}) unsupported";
        echo "=====================================";
        continue;
    elif [[ "${PYBIN}" == *"33"* ]]; then
        echo "Python 3.3 (via ${PYBIN}) unsupported";
        echo "=====================================";
        continue;
    elif [[ "${PYBIN}" == *"34"* ]]; then
        echo "Python 3.4 (via ${PYBIN}) unsupported";
        echo "=====================================";
        continue;
    else
        echo "Building with version: ${PYBIN}";
    fi
    "${PYBIN}/pip" install -r /io/scripts/manylinux/dev-requirements.txt
    "${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/${PKG_NAME}*.whl; do
    auditwheel repair "${whl}" -w /io/wheelhouse/
done
