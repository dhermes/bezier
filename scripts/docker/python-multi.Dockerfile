FROM ubuntu:16.04

ENV DEBIAN_FRONTEND noninteractive

# Ensure local Python is preferred over distribution Python.
ENV PATH /usr/local/bin:$PATH

# http://bugs.python.org/issue19846
# At the moment, setting "LANG=C" on a Linux system fundamentally breaks
# Python 3.
ENV LANG C.UTF-8

# Install dependencies.
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    git \
    libreadline-dev \
    libssl-dev \
    ssh \
    tcl \
    tcl-dev \
    tk \
    tk-dev \
    wget \
    lcov \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && rm -rf /var/lib/apt/lists/* \
  && rm -f /var/cache/apt/archives/*.deb

# Install the desired versions of Python.
RUN for PYTHON_VERSION in 2.7.13 3.5.4 3.6.2; do \
  set -ex \
    && wget --no-check-certificate -O python-${PYTHON_VERSION}.tar.xz "https://www.python.org/ftp/python/${PYTHON_VERSION%%[a-z]*}/Python-$PYTHON_VERSION.tar.xz" \
    && wget --no-check-certificate -O python-${PYTHON_VERSION}.tar.xz.asc "https://www.python.org/ftp/python/${PYTHON_VERSION%%[a-z]*}/Python-$PYTHON_VERSION.tar.xz.asc" \
    && export GNUPGHOME="$(mktemp -d)" \
    && gpg --keyserver ha.pool.sks-keyservers.net --recv-keys \
      # 2.7.13 (REF: https://github.com/docker-library/python/blob/c9954b06c8b178d7888bc1626bed5a14e43a9203/2.7/stretch/Dockerfile#L16)
      C01E1CAD5EA2C4F0B8E3571504C367C218ADD4FF \
      # 3.5.4 (REF: https://github.com/docker-library/python/blob/6ebbaa8a56cdf4021c78e87b3872be3861ac072a/3.5/jessie/Dockerfile#L22)
      97FC712E4C024BBEA48A61ED3A5CA953F73C700D \
      # 3.6.2 (REF: https://github.com/docker-library/python/blob/c9954b06c8b178d7888bc1626bed5a14e43a9203/3.6/stretch/Dockerfile#L22)
      0D96DF4D4110E5C43FBFB17F2D347EA6AA65421D \
    && gpg --batch --verify python-${PYTHON_VERSION}.tar.xz.asc python-${PYTHON_VERSION}.tar.xz \
    && rm -r "$GNUPGHOME" python-${PYTHON_VERSION}.tar.xz.asc \
    && mkdir -p /usr/src/python-${PYTHON_VERSION} \
    && tar -xJC /usr/src/python-${PYTHON_VERSION} --strip-components=1 -f python-${PYTHON_VERSION}.tar.xz \
    && rm python-${PYTHON_VERSION}.tar.xz \
    && cd /usr/src/python-${PYTHON_VERSION} \
    && ./configure \
      --enable-shared \
      # This works only on Python 2.7 and throws a warning on every other
      # version, but seems otherwise harmless.
      --enable-unicode=ucs4 \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
  ; done \
  && rm -rf /usr/src/python* \
  && rm -rf ~/.cache/

# REF: https://github.com/docker-library/pypy/blob/19f27e6e89012de63a5cb39e5e1d011651fd3a61/2/Dockerfile
ENV PYPY_VERSION 5.8.0
RUN set -ex; \
  \
# this "case" statement is generated via "update.sh"
    dpkgArch="$(dpkg --print-architecture)"; \
    case "${dpkgArch##*-}" in \
# amd64
    amd64) pypyArch='linux64'; sha256='6274292d0e954a2609b15978cde6efa30942ba20aa5d2acbbf1c70c0a54e9b1e' ;; \
# arm32v5
    armel) pypyArch='linux-armel'; sha256='28b7fd0cc7418ffc66c71520728e87941be40ebf4b82675c57e25598a2a702b0' ;; \
# arm32v7
    armhf) pypyArch='linux-armhf-raring'; sha256='ddceca9c5c9a456d4bf1beab177660adffbbdf255a922244e1cc05f20318be46' ;; \
# i386
    i386) pypyArch='linux32'; sha256='a0b125a5781f7e5ddfc3baca46503b14f4ee6a0e234e8d72bfcf3afdf4120bef' ;; \
    *) echo >&2 "error: current architecture ($dpkgArch) does not have a corresponding PyPy $PYPY_VERSION binary release"; exit 1 ;; \
  esac; \
  \
  wget -O pypy.tar.bz2 "https://bitbucket.org/pypy/pypy/downloads/pypy2-v${PYPY_VERSION}-${pypyArch}.tar.bz2"; \
  echo "$sha256 *pypy.tar.bz2" | sha256sum -c; \
  tar -xjC /usr/local --strip-components=1 -f pypy.tar.bz2; \
  rm pypy.tar.bz2; \
  \
  pypy --version

# Install pip on Python 3.6 only.
# If the environment variable is called "PIP_VERSION", pip explodes with
# "ValueError: invalid truth value '<VERSION>'"
ENV PYTHON_PIP_VERSION 9.0.1
RUN wget --no-check-certificate -O /tmp/get-pip.py 'https://bootstrap.pypa.io/get-pip.py' \
  && python3.6 /tmp/get-pip.py "pip==$PYTHON_PIP_VERSION" \
  && rm /tmp/get-pip.py \

  # we use "--force-reinstall" for the case where the version of pip we're trying to install is the same as the version bundled with Python
  # ("Requirement already up-to-date: pip==8.1.2 in /usr/local/lib/python3.6/site-packages")
  # https://github.com/docker-library/python/pull/143#issuecomment-241032683
  && pip3 install --no-cache-dir --upgrade --force-reinstall "pip==$PYTHON_PIP_VERSION" \

  # then we use "pip list" to ensure we don't have more than one pip version installed
  # https://github.com/docker-library/python/pull/100
  && [ "$(pip list |tac|tac| awk -F '[ ()]+' '$1 == "pip" { print $2; exit }')" = "$PYTHON_PIP_VERSION" ]

CMD ["python3.6"]
