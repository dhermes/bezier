# `0.11.0`

## CircleCI

Build `1588`

### Workflow

| Filename      | Step                                                    |
|---------------|---------------------------------------------------------|
| `1588-00.txt` | Spin up Environment                                     |
| `1588-01.txt` | Checkout code                                           |
| `1588-02.txt` | Check that all Cython generated files have been updated |
| `1588-03.txt` | Unit tests in Python 3.6                                |
| `1588-04.txt` | Unit tests in Python 3.7                                |
| `1588-05.txt` | Unit tests in pypy3                                     |
| `1588-06.txt` | Unit tests AND line coverage in Python 3.8              |
| `1588-07.txt` | Fortran unit tests                                      |
| `1588-08.txt` | Functional tests in Python 3.8                          |
| `1588-09.txt` | Run all doctests                                        |
| `1588-10.txt` | Build docs                                              |
| `1588-11.txt` | Lint code for style issues                              |
| `1588-12.txt` | Check that test case examples are valid for JSON schema |
| `1588-13.txt` | Check the journaled `build_ext` commands                |
| `1588-14.txt` | Upload coverage to coveralls                            |

## Travis CI

Build `784` / `635799581`

### Build Matrix

| Job ID  | `PY_VERSION` |
|---------|--------------|
| `784.1` | `3.6`        |
| `784.2` | `3.7`        |
| `784.3` | `3.8`        |

## AppVeyor

Build `1.0.1269.master`

### Build Matrix

| Job ID             | `NOX_SESSION`                       |
|--------------------|-------------------------------------|
| `3xkhknuswrq6mydb` | `unit-3.6-32`                       |
| `e2y53ec97un99k6i` | `unit-3.6`                          |
| `c6xfnjfqxeiocx9y` | `unit-3.7-32`                       |
| `bwxuf6afw017mw9a` | `unit-3.7`                          |
| `hscrxnqqha8lrekc` | `unit-3.8-32`                       |
| `lg471nidpa8sl220` | `cover`                             |
| `124axf4hyhn9ghd8` | `functional-3.8`                    |
| `xjq4k3vxnr9q8qqb` | `doctest`                           |
| `c8ulq9xf6tw35x23` | `check_journal(machine='appveyor')` |
