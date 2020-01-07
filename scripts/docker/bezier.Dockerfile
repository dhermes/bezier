FROM dhermes/python-multi

# Install the current versions of nox and NumPy.
RUN python3.8 -m pip install --no-cache-dir \
  "colorlog == 4.1.0" \
  "nox == 2019.11.9" \
  "numpy == 1.18.1" \
  "py == 1.8.1" \
  "virtualenv == 16.7.9"

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

# Build NumPy and SciPy wheels for PyPy 3 since it takes a bit of time.
ENV WHEELHOUSE=/wheelhouse
RUN mkdir ${WHEELHOUSE}
RUN set -ex \
  && virtualenv --python=pypy3 pypy3-env \
  && pypy3-env/bin/python -m pip install --upgrade pip wheel \
  && pypy3-env/bin/python -m pip wheel --wheel-dir=${WHEELHOUSE} "numpy == 1.18.1" \
  && pypy3-env/bin/python -m pip install ${WHEELHOUSE}/numpy*.whl \
  && pypy3-env/bin/python -m pip wheel --wheel-dir=${WHEELHOUSE} "scipy == 1.4.1" \
  && rm -fr pypy3-env
