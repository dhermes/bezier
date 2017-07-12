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

import glob
import os
import sys

import nox


NUMPY = 'numpy'
MOCK_DEP = 'mock >= 1.3.0'
BASE_DEPS = (
    MOCK_DEP,
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

    env = {'MATPLOTLIBRC': 'test'}
    return local_deps, env


@nox.session
@nox.parametrize('python_version', ['2.7', '3.5', '3.6', PYPY])
def unit_tests(session, python_version):
    if python_version == PYPY:
        session.interpreter = PYPY
        local_deps, env = pypy_setup(BASE_DEPS)
    else:
        session.interpreter = 'python{}'.format(python_version)
        local_deps = BASE_DEPS + ('scipy',)
        env = None

    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Run py.test against the unit tests.
    run_args = ['py.test'] + session.posargs + [get_path('tests')]
    session.run(*run_args, env=env)


def functional_env():
    return {'PYTHONPATH': 'functional_tests'}


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
    session.run(*run_args, env=functional_env())


@nox.session
@nox.parametrize('python_version', ['2.7', '3.5', '3.6', PYPY])
def functional(session, python_version):
    local_deps = BASE_DEPS
    if python_version == PYPY:
        session.interpreter = PYPY
        local_deps, env = pypy_setup(local_deps)
    else:
        session.interpreter = 'python{}'.format(python_version)
        env = {}

    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Run py.test against the functional tests.
    run_args = ['py.test'] + session.posargs + [get_path('functional_tests')]
    env.update(functional_env())
    session.run(*run_args, env=env)


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
    # Install all dependencies.
    local_deps = DOCS_DEPS + (MOCK_DEP,)
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Run the script for building docs and running doctests.
    run_args = get_doctest_args(session)
    session.run(*run_args, env={'NO_IMAGES': 'True'})


@nox.session
def docs_images(session):
    session.interpreter = SINGLE_INTERP
    # Install all dependencies.
    local_deps = DOCS_DEPS
    local_deps += ('matplotlib >= 2.0.0', MOCK_DEP, 'seaborn', 'pytest')
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    # Use custom RC-file for matplotlib.
    env = {'MATPLOTLIBRC': 'docs'}

    # Run the script for generating images for docs.
    run_args = get_doctest_args(session)
    session.run(*run_args, env=env)

    # Run the functional tests with --save-plot.
    fnl_tests_glob = get_path('functional_tests', 'test_*.py')
    modules_to_run = glob.glob(fnl_tests_glob)
    # Generate images for ``curve_intersections.json``.
    curve_curve_filename = get_path(
        'functional_tests', 'make_curve_intersection_images.py')
    modules_to_run.append(curve_curve_filename)

    for filename in modules_to_run:
        session.run('python', filename, '--save-plot', env=env)


@nox.session
def lint(session):
    session.interpreter = SINGLE_INTERP
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
        env=functional_env(),
    )
