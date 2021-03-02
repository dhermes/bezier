# `0.9.0`

## CircleCI

Build `1405`

### Workflow

| Filename      | Step                                                    |
|---------------|---------------------------------------------------------|
| `1405-00.txt` | Spin up Environment                                     |
| `1405-01.txt` | Checkout code                                           |
| `1405-02.txt` | Check that all Cython generated files have been updated |
| `1405-03.txt` | Unit tests in Python 2.7                                |
| `1405-04.txt` | Unit tests in Python 3.6                                |
| `1405-05.txt` | Unit tests in pypy                                      |
| `1405-06.txt` | Unit tests AND line coverage in Python 3.7              |
| `1405-07.txt` | Fortran unit tests                                      |
| `1405-08.txt` | Functional tests in Python 3.7                          |
| `1405-09.txt` | Run all doctests                                        |
| `1405-10.txt` | Build docs                                              |
| `1405-11.txt` | Lint code for style issues                              |
| `1405-12.txt` | Check the journaled `build_ext` commands                |
| `1405-13.txt` | Upload coverage to coveralls                            |

## Travis CI

Build `602` / `448688540`

### Build Matrix

| Job ID  | `PY_VERSION` |
|---------|--------------|
| `602.1` | `2.7`        |
| `602.2` | `3.6`        |
| `602.3` | `3.7`        |

## AppVeyor

Build `1.0.1078.master`

### Build Matrix

| Job ID             | `NOX_SESSION`                       |
|--------------------|-------------------------------------|
| `8c27uuqnyf9nthej` | `unit-2.7-32`                       |
| `09o7ke9251xcg3od` | `unit-2.7`                          |
| `oty7w5o3wm2sepop` | `unit-3.6-32`                       |
| `eo6v2ngjf8y586fl` | `unit-3.6`                          |
| `2092dq7o3syhwwmu` | `unit-3.7-32`                       |
| `4ivdqgiviankvf7t` | `cover`                             |
| `xctpbqrq4bwt7rkl` | `functional-3.7`                    |
| `6jnjfq47kxgvaxbg` | `doctest`                           |
| `4f22pqhpewuy9767` | `check_journal(machine='appveyor')` |
