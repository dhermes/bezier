# `2020.5.19`

## CircleCI

Build `1800`

### Workflow

| Filename      | Step                                                    |
|---------------|---------------------------------------------------------|
| `1800-00.txt` | Spin up Environment                                     |
| `1800-01.txt` | Preparing Environment Variables                         |
| `1800-02.txt` | Checkout code                                           |
| `1800-03.txt` | Setup a remote Docker engine                            |
| `1800-04.txt` | Check that all Cython generated files have been updated |
| `1800-05.txt` | Unit tests in Python 3.6                                |
| `1800-06.txt` | Unit tests in Python 3.7                                |
| `1800-07.txt` | Unit tests in pypy3                                     |
| `1800-08.txt` | Unit tests AND line coverage in Python 3.8              |
| `1800-09.txt` | Fortran unit tests                                      |
| `1800-10.txt` | Functional tests in Python 3.8                          |
| `1800-11.txt` | Run all doctests                                        |
| `1800-12.txt` | Build docs                                              |
| `1800-13.txt` | Lint code for style issues                              |
| `1800-14.txt` | Check that test case examples are valid for JSON schema |
| `1800-15.txt` | Upload coverage to coveralls                            |

## Travis CI

Build `976` / `689100504`

### Build Matrix

| Job ID  | `PY_VERSION` |
|---------|--------------|
| `976.1` | `3.6`        |
| `976.2` | `3.7`        |
| `976.3` | `3.8`        |

## AppVeyor

Build `1.0.1467.master`

### Build Matrix

| Job ID             | `NOX_SESSION`    |
|--------------------|------------------|
| `e1u0f19vbs8k58t1` | `unit-3.6-32`    |
| `px3kqs8ce025qldw` | `unit-3.6`       |
| `m411s8ollyb3rxl4` | `unit-3.7-32`    |
| `5q8wpshe4h98t8e1` | `unit-3.7`       |
| `9n4kw45a4qra6ft9` | `unit-3.8-32`    |
| `7r14756oxc8jihfd` | `cover`          |
| `ltv6jr7hc6r88xr5` | `functional-3.8` |
| `npye5afbbaw5smjb` | `doctest`        |
