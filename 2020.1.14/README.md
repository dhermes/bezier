# `2020.1.14`

## CircleCI

Build `1609`

### Workflow

| Filename      | Step                                                    |
|---------------|---------------------------------------------------------|
| `1609-00.txt` | Spin up Environment                                     |
| `1609-01.txt` | Checkout code                                           |
| `1609-02.txt` | Check that all Cython generated files have been updated |
| `1609-03.txt` | Unit tests in Python 3.6                                |
| `1609-04.txt` | Unit tests in Python 3.7                                |
| `1609-05.txt` | Unit tests in pypy3                                     |
| `1609-06.txt` | Unit tests AND line coverage in Python 3.8              |
| `1609-07.txt` | Fortran unit tests                                      |
| `1609-08.txt` | Functional tests in Python 3.8                          |
| `1609-09.txt` | Run all doctests                                        |
| `1609-10.txt` | Build docs                                              |
| `1609-11.txt` | Lint code for style issues                              |
| `1609-12.txt` | Check that test case examples are valid for JSON schema |
| `1609-13.txt` | Check the journaled `build_ext` commands                |
| `1609-14.txt` | Upload coverage to coveralls                            |

## Travis CI

Build `810` / `636950692`

### Build Matrix

| Job ID  | `PY_VERSION` |
|---------|--------------|
| `810.1` | `3.6`        |
| `810.2` | `3.7`        |
| `810.3` | `3.8`        |

## AppVeyor

Build `1.0.1295.master`

### Build Matrix

| Job ID             | `NOX_SESSION`                       |
|--------------------|-------------------------------------|
| `g5vl98ufojcoeda1` | `unit-3.6-32`                       |
| `p857stfp5m8dj5jt` | `unit-3.6`                          |
| `vh70352kac8s2m25` | `unit-3.7-32`                       |
| `6fbuh9hi8cwv0e7u` | `unit-3.7`                          |
| `m71hah8dsefluw67` | `unit-3.8-32`                       |
| `c7kqldvjoybf6coo` | `cover`                             |
| `p8m73530jys0ogj6` | `functional-3.8`                    |
| `u2ria5cnuketmyri` | `doctest`                           |
| `nj1ml8aq8jas20v0` | `check_journal(machine='appveyor')` |
