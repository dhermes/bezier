This is a collection of useful `docker` commands relevant to
this directory. Some of them are just generic `docker` commands
for making sure the local system is in a good state (useful
for those who don't use `docker` on a regular basis).

## Specific

```
$ docker build \
>   --file python-multi.Dockerfile \
>   --tag dhermes/python-multi:latest \
>   .
$ docker build \
>   --file bezier.Dockerfile \
>   --tag dhermes/bezier:latest \
>   .
$ docker run \
>   --rm \
>   --tty \
>   --interactive \
>   --volume $(pwd):/var/code/bezier/
>   dhermes/bezier:latest \
>   /bin/bash
```

## Generic

```
$ docker system prune  # Clean-up
$ docker image  prune  # Clean-up
$ docker container ls  # Running containers
$ docker images        # List installed/cached images
$ docker image rm ${IMAGE_ID}
$ docker image rm ${REPOSITORY}:${TAG}
```
