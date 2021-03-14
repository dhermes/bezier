# `0.10.0`

## CircleCI

Build `1479`

### Workflow

| Filename      | Step                                                    |
|---------------|---------------------------------------------------------|
| `1479-00.txt` | Spin up Environment                                     |
| `1479-01.txt` | Checkout code                                           |
| `1479-02.txt` | Check that all Cython generated files have been updated |
| `1479-03.txt` | Unit tests in Python 3.6                                |
| `1479-04.txt` | Unit tests in pypy3                                     |
| `1479-05.txt` | Unit tests AND line coverage in Python 3.7              |
| `1479-06.txt` | Fortran unit tests                                      |
| `1479-07.txt` | Functional tests in Python 3.7                          |
| `1479-08.txt` | Run all doctests                                        |
| `1479-09.txt` | Build docs                                              |
| `1479-10.txt` | Lint code for style issues                              |
| `1479-11.txt` | Check that test case examples are valid for JSON schema |
| `1479-12.txt` | Check the journaled `build_ext` commands                |
| `1479-13.txt` | Upload coverage to coveralls                            |

## Travis CI

Build `671` / `571656947`

### Build Matrix

| Job ID  | `PY_VERSION` |
|---------|--------------|
| `671.1` | `3.6`        |
| `671.2` | `3.7`        |

## AppVeyor

Build `1.0.1147.master`

### Build Matrix

| Job ID             | `NOX_SESSION`                       |
|--------------------|-------------------------------------|
| `k7mhfrbxcg5xs6pl` | `unit-3.6-32`                       |
| `3x9qxt36p9px5xl5` | `unit-3.6`                          |
| `59y12m1x0orjqo89` | `unit-3.7-32`                       |
| `uo969cn2gn6fkxux` | `cover`                             |
| `gdhubo9vtvjuide6` | `functional-3.7`                    |
| `buvf8xsaiuy380ju` | `doctest`                           |
| `7ljts8c85rvjxx5i` | `check_journal(machine='appveyor')` |
