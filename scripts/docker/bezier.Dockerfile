FROM dhermes/python-multi:2017.08

# Install the current versions of nox and NumPy.
RUN pip install --no-cache-dir \
  colorlog==2.10.0 \
  nox-automation==0.17.0 \
  numpy==1.13.1 \
  py==1.4.34 \
  six==1.10.0 \
  virtualenv==15.1.0

# Install `gfortran`.
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    gfortran
