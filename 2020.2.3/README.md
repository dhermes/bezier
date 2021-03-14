# `2020.2.3`

## CircleCI

Build `1722`

### Workflow

| Filename      | Step                                                    |
|---------------|---------------------------------------------------------|
| `1722-00.txt` | Spin up Environment                                     |
| `1722-01.txt` | Checkout code                                           |
| `1722-02.txt` | Setup a remote Docker engine                            |
| `1722-03.txt` | Check that all Cython generated files have been updated |
| `1722-04.txt` | Unit tests in Python 3.6                                |
| `1722-05.txt` | Unit tests in Python 3.7                                |
| `1722-06.txt` | Unit tests in pypy3                                     |
| `1722-07.txt` | Unit tests AND line coverage in Python 3.8              |
| `1722-08.txt` | Fortran unit tests                                      |
| `1722-09.txt` | Functional tests in Python 3.8                          |
| `1722-10.txt` | Run all doctests                                        |
| `1722-11.txt` | Build docs                                              |
| `1722-12.txt` | Lint code for style issues                              |
| `1722-13.txt` | Check that test case examples are valid for JSON schema |
| `1722-14.txt` | Upload coverage to coveralls                            |

## Travis CI

Build `913` / `645805500`

### Build Matrix

| Job ID  | `PY_VERSION` |
|---------|--------------|
| `913.1` | `3.6`        |
| `913.2` | `3.7`        |
| `913.3` | `3.8`        |

## AppVeyor

Build `1.0.1404.master`

### Build Matrix

| Job ID             | `NOX_SESSION`                       |
|--------------------|-------------------------------------|
| `xxyq1d3jsgv5ws8i` | `unit-3.6-32`                       |
| `yamdne12hctkeoa2` | `unit-3.6`                          |
| `63dqysa8v1008mjr` | `unit-3.7-32`                       |
| `dtfn0fmvhaho7ggw` | `unit-3.7`                          |
| `8kd1jm5ys7wnt67i` | `unit-3.8-32`                       |
| `gxvg4jr034lu95vp` | `cover`                             |
| `pogm7wkcaxdafobw` | `functional-3.8`                    |
| `uft7eb4840u5s7ni` | `doctest`                           |
