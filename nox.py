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

import os

import nox


BASE_DEPS = (
    'mock >= 1.3.0',
    'numpy',
    'pytest',
)
DOCS_DEPS = (
    os.path.join('docs', 'requirements.txt'),
)
NOX_DIR = os.path.abspath(os.path.dirname(__file__))


def get_path(*names):
    return os.path.join(NOX_DIR, *names)


@nox.session
@nox.parametrize('python_version', ['2.7', '3.5', '3.6'])
def unit_tests(session, python_version):
    session.interpreter = 'python{}'.format(python_version)

    # Install all test dependencies.
    local_deps = BASE_DEPS + ('scipy',)
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

    # Run py.test against the unit tests.
    run_args = ['py.test', '--cov=bezier', '--cov=tests']
    run_args += session.posargs
    run_args += [
        get_path('tests'),
        get_path('functional_tests', 'test_segment_box.py'),
    ]
    session.run(*run_args)


@nox.session
def functional(session):
    session.interpreter = 'python2.7'

    # Install all test dependencies.
    session.install(*BASE_DEPS)
    # Install this package.
    session.install('.')

    # Run py.test against the unit tests.
    run_args = ['py.test'] + session.posargs + [get_path('functional_tests')]
    session.run(*run_args)
