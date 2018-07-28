FROM dhermes/python-multi

# Install the current versions of nox and NumPy.
RUN python3.6 -m pip install --no-cache-dir \
  colorlog==2.10.0 \
  nox-automation==0.19.1 \
  numpy==1.15.0 \
  py==1.5.4 \
  six==1.11.0 \
  virtualenv==16.0.0

# Install `gfortran` (for Fortran extensions), `libatlas-base-dev`,
# `libblas-dev`, `liblapack-dev` (for SciPy) and `lcov` for
# Fortran code coverage.
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    gfortran \
    libatlas-base-dev \
    libblas-dev \
    liblapack-dev \
    lcov \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -f /var/cache/apt/archives/*.deb

# Build NumPy and SciPy wheels for PyPy since it takes a bit of time.
RUN for PYPY in pypy pypy3; do \
  set -ex \
    && virtualenv --python=${PYPY} pypy-env \
    && pypy-env/bin/python -m pip install --upgrade pip wheel \
    && mkdir /wheelhouse-${PYPY} \
    && pypy-env/bin/python -m pip wheel --wheel-dir=/wheelhouse-${PYPY} numpy==1.15.0 \
    && pypy-env/bin/python -m pip install /wheelhouse-${PYPY}/numpy*.whl \
    && pypy-env/bin/python -m pip wheel --wheel-dir=/wheelhouse-${PYPY} scipy==1.1.0 \
    && rm -fr pypy-env \
  ; done

# Combine the version specific wheelhouses into one directory.
RUN mkdir /wheelhouse \
  && mv /wheelhouse-pypy/* /wheelhouse \
  && mv /wheelhouse-pypy3/* /wheelhouse \
  && rm -fr /wheelhouse-pypy /wheelhouse-pypy3
