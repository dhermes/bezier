# `2021.2.12`

## CircleCI

Build `1900`

### Workflow

| Filename      | Step                                                    |
|---------------|---------------------------------------------------------|
| `1900-00.txt` | Spin up Environment                                     |
| `1900-01.txt` | Preparing Environment Variables                         |
| `1900-02.txt` | Checkout code                                           |
| `1900-03.txt` | Setup a remote Docker engine                            |
| `1900-04.txt` | Check that all Cython generated files have been updated |
| `1900-05.txt` | Unit tests in Python 3.7                                |
| `1900-06.txt` | Unit tests in Python 3.8                                |
| `1900-07.txt` | Unit tests in pypy3                                     |
| `1900-08.txt` | Unit tests AND line coverage in Python 3.9              |
| `1900-09.txt` | Functional tests in Python 3.9                          |
| `1900-10.txt` | Run all doctests                                        |
| `1900-11.txt` | Build docs                                              |
| `1900-12.txt` | Lint code for style issues                              |
| `1900-13.txt` | Check that test case examples are valid for JSON schema |
| `1900-14.txt` | Upload coverage to coveralls                            |

## GitHub Actions (macOS)

Run `10` / `564975428`

### Build Matrix

| Run ID       | `matrix.python-version` |
|--------------|-------------------------|
| `1895890169` | `3.7`                   |
| `1895890174` | `3.8`                   |
| `1895890182` | `3.9`                   |

## AppVeyor

Build `1.0.1539.main`

### Build Matrix

| Job ID             | `NOX_SESSION`    |
|--------------------|------------------|
| `4614ljqd8kiljd70` | `unit-3.7-32`    |
| `e48dtste9u0u83g7` | `unit-3.7`       |
| `86a57eec8psphsua` | `unit-3.8-32`    |
| `lbtyb4q35b3v28k5` | `unit-3.8`       |
| `5jbjv6k0t6var1a2` | `unit-3.9-32`    |
| `3oi5to71jeamfiep` | `cover`          |
| `kc28cs8xahh0we53` | `functional-3.9` |
| `b8gn9el9l96rvo86` | `doctest`        |
