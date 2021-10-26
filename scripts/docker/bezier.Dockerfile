FROM dhermes/python-multi:latest

# Install the current versions of CMake, nox and NumPy.
RUN python3.10 -m pip install --no-cache-dir \
  "cmake == 3.21.3" \
  "nox == 2021.10.1" \
  "numpy == 1.21.3" \
  "argcomplete==1.12.3" \
  "backports.entry-points-selectable==1.1.0" \
  "colorlog==6.5.0" \
  "distlib==0.3.3" \
  "filelock==3.3.1" \
  "packaging==21.0" \
  "platformdirs==2.4.0" \
  "py==1.10.0" \
  "pyparsing==3.0.1" \
  "six==1.16.0" \
  "virtualenv==20.9.0"

# Install `gfortran` (for building the Fortran code used by the binary
# extension), `libatlas-base-dev`, `libblas-dev`, `liblapack-dev` (for SciPy)
# and `lcov` for Fortran code coverage.
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
  && pypy3-env/bin/python -m pip wheel --wheel-dir=${WHEELHOUSE} "numpy == 1.21.3" \
  && pypy3-env/bin/python -m pip install ${WHEELHOUSE}/numpy*.whl \
  && pypy3-env/bin/python -m pip wheel --wheel-dir=${WHEELHOUSE} "scipy == 1.7.1" \
  && rm -fr pypy3-env

# Install Docker CLI (used to build `manylinux` wheel for `nox --session doctest`).
# Via: https://docs.docker.com/install/linux/docker-ce/ubuntu/
# Includes a diagnostic-only use of `apt-key fingerprint`.
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    apt-transport-https \
    ca-certificates \
    curl \
    gnupg-agent \
    software-properties-common \
  && curl -fsSL https://download.docker.com/linux/ubuntu/gpg | apt-key add - \
  && apt-key fingerprint 0EBFCD88 \
  && add-apt-repository \
    "deb [arch=amd64] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) \
    stable" \
  && apt-get update \
  && apt-get install -y docker-ce-cli \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -f /var/cache/apt/archives/*.deb
