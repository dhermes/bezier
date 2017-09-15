# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     https://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function

import glob
import os
import tempfile

import nox


NUMPY = 'numpy'
MOCK_DEP = 'mock >= 1.3.0'
SEABORN_DEP = 'seaborn >= 0.8'
BASE_DEPS = (
    MOCK_DEP,
    NUMPY,
    'pytest',
)
NOX_DIR = os.path.abspath(os.path.dirname(__file__))
WHEELHOUSE = os.environ.get('WHEELHOUSE')
DOCS_DEPS = (
    '--requirement',
    os.path.join(NOX_DIR, 'docs', 'requirements.txt'),
)
SINGLE_INTERP = 'python3.6'
PYPY = 'pypy'
JOURNAL_PATHS = {
    'circleci': os.path.join('.circleci', 'expected_journal.txt'),
    'travis-osx': os.path.join('scripts', 'osx', 'travis_journal.txt'),
}


def get_path(*names):
    return os.path.join(NOX_DIR, *names)


def pypy_setup(local_deps, session):
    if WHEELHOUSE is not None:
        # Remove NumPy from dependencies.
        local_deps = list(local_deps)
        local_deps.remove(NUMPY)
        local_deps = tuple(local_deps)
        # Install from the pre-built wheel.
        session.install(
            '--use-wheel',
            '--no-index',
            '--find-links',
            WHEELHOUSE,
            NUMPY,
        )

    env = {'MATPLOTLIBRC': 'test'}
    return local_deps, env


@nox.session
@nox.parametrize('check', [True, False])
def update_generated(session, check):
    name = 'update-{}'.format(check)
    # NOTE: ``nox`` requires virtualenv_dirname to be lowercase.
    session.virtualenv_dirname = name.lower()

    # Update Cython generated source code.
    session.interpreter = SINGLE_INTERP

    # Install all dependencies.
    session.install('Cython')

    if check:
        command = get_path('scripts', 'remove_cython_files.py')
        session.run('python', command)

    pyx_glob = get_path('src', 'bezier', '*.pyx')
    for pyx_module in glob.glob(pyx_glob):
        session.run('cython', pyx_module)

    command = get_path('scripts', 'clean_cython.py')
    c_glob = get_path('src', 'bezier', '*.c')
    for c_source in glob.glob(c_glob):
        session.run(
            'python', command,
            '--filename', c_source,
            '--virtualenv-dirname', session.virtualenv_dirname,
        )

    if check:
        command = get_path('scripts', 'cython_update_check.py')
        session.run('python', command)


@nox.session
@nox.parametrize('python_version', ['2.7', '3.5', '3.6', PYPY])
def unit_tests(session, python_version):
    if python_version == PYPY:
        session.interpreter = PYPY
        local_deps, env = pypy_setup(BASE_DEPS, session)
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
    if python_version == PYPY:
        session.interpreter = PYPY
        local_deps, env = pypy_setup(BASE_DEPS, session)
    else:
        session.interpreter = 'python{}'.format(python_version)
        local_deps = BASE_DEPS
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
    local_deps += ('matplotlib >= 2.0.0', MOCK_DEP, SEABORN_DEP, 'pytest')
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
    # Generate images for ``curve_intersections.json`` and
    # ``surface_intersections.json``.
    modules_to_run.extend((
        get_path('functional_tests', 'make_curve_curve_images.py'),
        get_path('functional_tests', 'make_surface_surface_images.py')
    ))

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

    # Run the script to check that the README and other docs are valid.
    check_path = get_path('scripts', 'check_doc_templates.py')
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
        '--max-module-lines=2582',
        get_path('src', 'bezier'),
    )
    # Run Pylint over the tests source.
    session.run(
        'pylint', '--rcfile', 'pylintrc',
        '--disable=missing-docstring',
        '--disable=protected-access',
        '--disable=too-many-public-methods',
        '--max-module-lines=2607',
        get_path('functional_tests'),
        get_path('tests'),
        env=functional_env(),
    )


@nox.session
@nox.parametrize('target', ['memory', 'time'])
def benchmark(session, target):
    session.interpreter = SINGLE_INTERP

    if target == 'memory':
        local_deps = (NUMPY, 'psutil', 'memory_profiler')
        test_fi = get_path('benchmarks', 'memory', 'test_curves.py')
        run_args = ('python', test_fi)
    elif target == 'time':
        local_deps = BASE_DEPS + ('pytest-benchmark',)
        test_fi = get_path('benchmarks', 'time', 'test_curves.py')
        run_args = ['py.test'] + session.posargs + [test_fi]

    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    session.install('.')

    session.run(*run_args, env=functional_env())


@nox.session
@nox.parametrize('machine', ['circleci', 'travis-osx'])
def check_journal(session, machine):
    session.virtualenv_dirname = 'journal-{}'.format(machine)
    session.interpreter = SINGLE_INTERP

    # Get a temporary file where the journal will be written.
    filehandle, journal_filename = tempfile.mkstemp(suffix='-journal.txt')
    os.close(filehandle)

    # Set the journal environment variable and install ``bezier``.
    session.install(NUMPY)  # Install requirement(s).
    env = {'BEZIER_JOURNAL': journal_filename}
    session.run('pip', 'install', '.', env=env)

    # Compare the expected file to the actual results.
    expected_journal = get_path(JOURNAL_PATHS[machine])
    session.run('diff', journal_filename, expected_journal)
