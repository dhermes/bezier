# Copyright 2017 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function

import os
import sys

import nox
import nox.command


NUMPY = 'numpy'
BASE_DEPS = (
    'mock >= 1.3.0',
    NUMPY,
    'pytest',
)
NOX_DIR = os.path.abspath(os.path.dirname(__file__))
DOCS_DEPS = (
    '--requirement',
    os.path.join(NOX_DIR, 'docs', 'requirements.txt'),
)
SINGLE_INTERP = 'python3.6'
PYPY = 'pypy'
PYPY_NUMPY = 'git+https://bitbucket.org/pypy/numpy.git'


def get_path(*names):
    return os.path.join(NOX_DIR, *names)


def pypy_setup(local_deps):
    local_deps = list(local_deps)
    local_deps.remove(NUMPY)
    local_deps.append(PYPY_NUMPY)
    local_deps = tuple(local_deps)

    if 'MATPLOTLIBRC' not in os.environ:
        reason = 'MATPLOTLIBRC=test/ must be set'
        print(reason, file=sys.stderr)
        raise nox.command.CommandFailed(reason=reason)

    return local_deps


@nox.session
@nox.parametrize('python_version', ['2.7', '3.5', '3.6', PYPY])
def unit_tests(session, python_version):
    if python_version == PYPY:
        session.interpreter = PYPY
        local_deps = pypy_setup(BASE_DEPS)
    else:
        session.interpreter = 'python{}'.format(python_version)
        local_deps = BASE_DEPS + ('scipy',)

    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Run py.test against the unit tests.
    run_args = ['py.test'] + session.posargs + [get_path('tests')]
    session.run(*run_args)


@nox.session
def cover(session):
    session.interpreter = 'python2.7'

    # Install all test dependencies.
    local_deps = BASE_DEPS + ('scipy', 'pytest-cov', 'coverage')
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Run py.test with coverage against the unit tests.
    run_args = ['py.test', '--cov=bezier', '--cov=tests']
    run_args += session.posargs
    run_args += [
        get_path('tests'),
        get_path('functional_tests', 'test_segment_box.py'),
    ]
    session.run(*run_args)


@nox.session
@nox.parametrize('python_version', ['2.7', '3.5', '3.6', PYPY])
def functional(session, python_version):
    local_deps = BASE_DEPS
    if python_version == PYPY:
        session.interpreter = PYPY
        local_deps = pypy_setup(local_deps)
    else:
        session.interpreter = 'python{}'.format(python_version)

    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Run py.test against the functional tests.
    run_args = ['py.test'] + session.posargs + [get_path('functional_tests')]
    session.run(*run_args)


@nox.session
def docs(session):
    session.interpreter = SINGLE_INTERP

    # Install all dependencies.
    session.install(*DOCS_DEPS)
    # Install this package.
    session.install('.')

    # Run the script for building docs.
    command = get_path('scripts', 'build_docs.sh')
    session.run(command)


def get_doctest_args(session):
    run_args = [
        'sphinx-build', '-W',
        '-b', 'doctest',
        '-d', get_path('docs', 'build', 'doctrees'),
        get_path('docs'),
        get_path('docs', 'build', 'doctest'),
    ]
    run_args += session.posargs
    return run_args


@nox.session
def doctest(session):
    session.interpreter = SINGLE_INTERP
    if 'NO_IMAGES' not in os.environ:
        reason = 'NO_IMAGES=True must be set'
        print(reason, file=sys.stderr)
        raise nox.command.CommandFailed(reason=reason)

    # Install all dependencies.
    session.install(*DOCS_DEPS)
    # Install this package.
    session.install('.')

    # Run the script for building docs and running doctests.
    run_args = get_doctest_args(session)
    session.run(*run_args)


@nox.session
def docs_images(session):
    session.interpreter = SINGLE_INTERP
    if 'MATPLOTLIBRC' not in os.environ:
        reason = 'MATPLOTLIBRC=docs/ must be set'
        print(reason, file=sys.stderr)
        raise nox.command.CommandFailed(reason=reason)

    # Install all dependencies.
    local_deps = DOCS_DEPS + ('matplotlib >= 2.0.0', 'seaborn')
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Run the script for generating images for docs.
    run_args = get_doctest_args(session)
    session.run(*run_args)


@nox.session
def lint(session):
    session.interpreter = SINGLE_INTERP
    if 'PYTHONPATH' not in os.environ:
        reason = 'PYTHONPATH env. var. must point to functional_tests/'
        print(reason, file=sys.stderr)
        raise nox.command.CommandFailed(reason=reason)

    # Install all dependencies.
    local_deps = BASE_DEPS + (
        'docutils',
        'flake8',
        'flake8-import-order',
        'matplotlib',
        'Pygments',
        'pylint',
    )
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Run the script to check that the README is valid.
    check_path = get_path('scripts', 'check_readme.py')
    session.run('python', check_path)
    # Run the script to check that setup.py is valid.
    setup_file = get_path('setup.py')
    session.run(
        'python', setup_file, 'check', '--metadata',
        '--restructuredtext', '--strict')
    # Run flake8 over the code to check import order.
    session.run(
        'flake8',
        '--import-order-style=google',
        '--application-import-names=bezier,tests',
        get_path('src', 'bezier'),
        get_path('tests'),
    )
    # Run Pylint over the library source.
    session.run(
        'pylint', '--rcfile', 'pylintrc',
        '--max-module-lines=2463',
        get_path('src', 'bezier'),
    )
    # Run Pylint over the tests source.
    session.run(
        'pylint', '--rcfile', 'pylintrc',
        '--disable=missing-docstring',
        '--disable=protected-access',
        '--disable=too-many-public-methods',
        '--max-module-lines=2324',
        get_path('functional_tests'),
        get_path('tests'),
    )
