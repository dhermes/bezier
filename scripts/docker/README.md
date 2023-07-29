## [`bezier`][1] Dockerfile

> Provides Python, NumPy, CMake, gfortran and nox. For automated testing of
> `bezier` Python package.

`bezier` is a Python library for interacting with B&#xe9;zier Curves and
Triangles.

This repository is intended to be used as a fully-functional environment for
installing / running / testing `bezier`. Provides:

-   CPython 3.9, 3.10, 3.11 and PyPy 3.9 (version 7.3.11)
-   NumPy pre-installed (needed for installing the binary extension)
-   `/wheelhouse` directory with NumPy and SciPy pre-built for PyPy 3
-   `gfortran` compiler
-   `nox` test runner
-   `lcov` line coverage tool
-   `libatlas-base-dev`, `libblas-dev` and `liblapack-dev` (needed by SciPy)
-   `cmake` for building `libbezier`
-   `docker` for running the Docker CLI (assuming `/var/run/docker.sock` has
    been shared from the host OS)

## Commands

This is a collection of useful `docker` commands relevant to
this directory. Some of them are just generic `docker` commands
for making sure the local system is in a good state (useful
for those who don't use `docker` on a regular basis).

### Specific

```
$ docker build \
>     --file bezier.Dockerfile \
>     --tag dhermes/bezier:latest \
>     .
$ docker run \
>     --rm \
>     --tty \
>     --interactive \
>     --volume /var/run/docker.sock:/var/run/docker.sock \
>     --volume $(git rev-parse --show-toplevel):/var/code/bezier/ \
>     --env MATPLOTLIBRC=tests \
>     --workdir /var/code/bezier/ \
>     dhermes/bezier:latest \
>     /bin/bash
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
