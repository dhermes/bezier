## [`bezier`][1] Dockerfile

> Provides Python, NumPy, gfortran and nox. For automated testing of bezier Python package.

`bezier` is a Python library for interacting with B&#xe9;zier Curves and Surfaces.

This repository is intended to be used as a fully-functional environment for installing / running / testing `bezier`. Provides:

- CPython 2.7, 3.5, 3.6 and PyPy 2.7
- NumPy pre-installed (needed for installing extensions)
- `/wheelhouse` directory with NumPy and SciPy pre-built for PyPy
- `gfortran` compiler
- `nox` test runner
- `lcov` line coverage tool
- `libatlas-dev`, `libblas-dev` and `liblapack-dev` (needed by SciPy)

## Commands

This is a collection of useful `docker` commands relevant to
this directory. Some of them are just generic `docker` commands
for making sure the local system is in a good state (useful
for those who don't use `docker` on a regular basis).

### Specific

```
$ docker build \
>   --file bezier.Dockerfile \
>   --tag dhermes/bezier:latest \
>   .
$ docker run \
>   --rm \
>   --tty \
>   --interactive \
>   --volume $(git rev-parse --show-toplevel):/var/code/bezier/ \
>   dhermes/bezier:latest \
>   /bin/bash
```

### Generic

```
$ docker system prune  # Clean-up
$ docker image  prune  # Clean-up
$ docker container ls  # Running containers
$ docker images        # List installed/cached images
$ docker image rm ${IMAGE_ID}
$ docker image rm ${REPOSITORY}:${TAG}
```

[1]: https://hub.docker.com/r/dhermes/bezier/
