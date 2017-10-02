FROM dhermes/python-multi

# Install the current versions of nox and NumPy.
RUN pip install --no-cache-dir \
  colorlog==2.10.0 \
  nox-automation==0.18.1 \
  numpy==1.13.3 \
  py==1.4.34 \
  six==1.11.0 \
  virtualenv==15.1.0

# Build a NumPy wheel for PyPy since it takes a bit of time
RUN virtualenv --python=pypy pypy-env \
  && pypy-env/bin/pip install --upgrade pip wheel \
  && mkdir /wheelhouse \
  && pypy-env/bin/pip wheel --wheel-dir=/wheelhouse numpy==1.13.3 \
  && rm -fr pypy-env

# Install `gfortran` and `lcov`.
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    gfortran \
    lcov
