# Unit Tests for Fortran Interface of `libbezier`

This is so that the Fortran source can evolve without having to have an
equivalent Python function for every piece of code.

Usage:

```
$ make
Makefile for Fortran unit and functional tests

Usage:
   make unit          Build and run unit tests
   make valgrind      Run unit tests with valgrind
   make lcov          Run unit tests with coverage
   make functional    Build and run functional tests
   make clean         Clean generated files

```
