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
import glob
import os
import shutil
import sys
import tempfile

import nox
import py.path

NUMPY_DEP = 'numpy >= 1.14.0'
IS_MAC_OS_X = sys.platform == 'darwin'
SCIPY_DEP = 'scipy >= 1.0.0'
MOCK_DEP = 'mock >= 1.3.0'
SEABORN_DEP = 'seaborn >= 0.8'
MATPLOTLIB_DEP = 'matplotlib >= 2.2.0'
BASE_DEPS = (NUMPY_DEP, 'pytest')
NOX_DIR = os.path.abspath(os.path.dirname(__file__))
WHEELHOUSE = os.environ.get('WHEELHOUSE')
DOCS_DEPS = (
    '--requirement', os.path.join(NOX_DIR, 'docs', 'requirements.txt')
)
SINGLE_INTERP = 'python3.6'
PYPY = 'pypy'
JOURNAL_PATHS = {
    'appveyor': os.path.join('appveyor', 'expected_journal.txt'),
    'circleci': os.path.join('.circleci', 'expected_journal.txt'),
    'travis-osx': os.path.join('scripts', 'osx', 'travis_journal.txt'),
}


def get_path(*names):
    return os.path.join(NOX_DIR, *names)


def pypy_setup(local_deps, session):
    if WHEELHOUSE is not None:
        # Remove NumPy from dependencies.
        local_deps = list(local_deps)
        local_deps.remove(NUMPY_DEP)
        local_deps = tuple(local_deps)
        # Install NumPy and SciPy from pre-built wheels.
        session.install(
            '--use-wheel',
            '--no-index',
            '--find-links',
            WHEELHOUSE,
            NUMPY_DEP,
            SCIPY_DEP,
        )
    return local_deps


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
            'python',
            command,
            '--filename',
            c_source,
            '--virtualenv-dirname',
            session.virtualenv_dirname,
        )
    if check:
        command = get_path('scripts', 'cython_update_check.py')
        session.run('python', command)


@nox.session
@nox.parametrize('py', ['2.7', '3.5', '3.6', PYPY])
def unit(session, py):
    if py == PYPY:
        session.interpreter = PYPY
        local_deps = pypy_setup(BASE_DEPS, session)
    else:
        session.interpreter = 'python{}'.format(py)
        local_deps = BASE_DEPS + (SCIPY_DEP,)
    if py in ('2.7', PYPY):
        local_deps += (MOCK_DEP,)
    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    session.install('.')
    # Run py.test against the unit tests.
    run_args = ['py.test'] + session.posargs + [get_path('tests', 'unit')]
    session.run(*run_args)


@nox.session
def cover(session):
    session.interpreter = SINGLE_INTERP
    # Install all test dependencies.
    local_deps = BASE_DEPS + (SCIPY_DEP, 'pytest-cov', 'coverage')
    session.install(*local_deps)
    # Install this package.
    session.install('.')
    # Run py.test with coverage against the unit tests.
    run_args = ['py.test', '--cov=bezier', '--cov=tests.unit']
    run_args += session.posargs
    run_args += [get_path('tests', 'unit')]
    session.run(*run_args)


@nox.session
@nox.parametrize('py', ['2.7', '3.5', '3.6', PYPY])
def functional(session, py):
    if py == PYPY:
        session.interpreter = PYPY
        local_deps = pypy_setup(BASE_DEPS, session)
    else:
        session.interpreter = 'python{}'.format(py)
        local_deps = BASE_DEPS
    if py in ('2.7', PYPY):
        local_deps += (MOCK_DEP,)
    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    session.install('.')
    # Run py.test against the functional tests.
    run_args = (
        ['py.test'] + session.posargs + [get_path('tests', 'functional')]
    )
    session.run(*run_args)


@nox.session
def docs(session):
    session.interpreter = SINGLE_INTERP
    # Install all dependencies.
    session.install(*DOCS_DEPS)
    # Install this package.
    env = {'BEZIER_NO_EXTENSIONS': 'True'}
    session.run('pip', 'install', '.', env=env)
    # Run the script for building docs.
    command = get_path('scripts', 'build_docs.sh')
    session.run(command, env=env)


def get_doctest_args(session):
    run_args = [
        'sphinx-build',
        '-W',
        '-b',
        'doctest',
        '-d',
        get_path('docs', 'build', 'doctrees'),
        get_path('docs'),
        get_path('docs', 'build', 'doctest'),
    ]
    run_args += session.posargs
    return run_args


@nox.session
def doctest(session):
    session.interpreter = SINGLE_INTERP
    # Install all dependencies.
    local_deps = DOCS_DEPS
    session.install(*local_deps)
    # Install this package.
    if IS_MAC_OS_X:
        command = get_path('scripts', 'osx', 'nox-install-for-doctest.sh')
        session.run(command)
    else:
        session.install('.')
    # Run the script for building docs and running doctests.
    run_args = get_doctest_args(session)
    session.run(*run_args)


@nox.session
def docs_images(session):
    session.interpreter = SINGLE_INTERP
    # Install all dependencies.
    local_deps = DOCS_DEPS
    local_deps += (MATPLOTLIB_DEP, SEABORN_DEP, 'pytest')
    session.install(*local_deps)
    # Install this package.
    session.install('.')
    # Use custom RC-file for matplotlib.
    env = {'MATPLOTLIBRC': 'docs', 'GENERATE_IMAGES': 'True'}
    # Run the script for generating images for docs.
    run_args = get_doctest_args(session)
    session.run(*run_args, env=env)
    # Run the functional tests with --save-plot.
    fnl_tests_glob = get_path('tests', 'functional', 'test_*.py')
    modules_to_run = glob.glob(fnl_tests_glob)
    # Generate images for ``curve_intersections.json`` and
    # ``surface_intersections.json``.
    modules_to_run.extend(
        (
            get_path('tests', 'functional', 'make_curve_curve_images.py'),
            get_path('tests', 'functional', 'make_surface_surface_images.py'),
        )
    )
    # Make sure that the root directory is on the Python path so that
    # ``tests`` is import-able.
    env['PYTHONPATH'] = get_path()
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
        MATPLOTLIB_DEP,
        'Pygments',
        'pylint',
        SEABORN_DEP,
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
        'python',
        setup_file,
        'check',
        '--metadata',
        '--restructuredtext',
        '--strict',
    )
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
        'pylint',
        '--rcfile',
        'pylintrc',
        '--max-module-lines=2891',
        get_path('src', 'bezier'),
    )
    # Run Pylint over the tests source.
    session.run(
        'pylint',
        '--rcfile',
        'pylintrc',
        '--disable=missing-docstring',
        '--disable=protected-access',
        '--disable=too-many-public-methods',
        '--max-module-lines=2368',
        get_path('tests'),
    )


@nox.session
@nox.parametrize('target', ['memory', 'time'])
def benchmark(session, target):
    session.interpreter = SINGLE_INTERP
    if target == 'memory':
        local_deps = (NUMPY_DEP, 'psutil', 'memory_profiler')
        test_fi1 = get_path('benchmarks', 'memory', 'test_curves.py')
        test_fi2 = get_path('benchmarks', 'memory', 'test_surfaces.py')
        all_run_args = [['python', test_fi1], ['python', test_fi2]]
    elif target == 'time':
        local_deps = BASE_DEPS + ('pytest-benchmark',)
        test_dir = get_path('benchmarks', 'time')
        all_run_args = [['py.test'] + session.posargs + [test_dir]]
    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    session.install('.')
    # NOTE: We need `tests` to be import-able.
    for run_args in all_run_args:
        session.run(*run_args, env={'PYTHONPATH': '.'})


@nox.session
@nox.parametrize('machine', ['appveyor', 'circleci', 'travis-osx'])
def check_journal(session, machine):
    session.virtualenv_dirname = 'journal-{}'.format(machine)
    session.interpreter = SINGLE_INTERP
    # Get a temporary file where the journal will be written.
    filehandle, journal_filename = tempfile.mkstemp(suffix='-journal.txt')
    os.close(filehandle)
    # Set the journal environment variable and install ``bezier``.
    session.install(NUMPY_DEP)  # Install requirement(s).
    env = {'BEZIER_JOURNAL': journal_filename}
    session.run('pip', 'install', '.', env=env)
    # Compare the expected file to the actual results.
    session.run(
        'python',
        get_path('scripts', 'post_process_journal.py'),
        '--journal-filename',
        journal_filename,
        '--machine',
        machine,
    )
    expected_journal = get_path(JOURNAL_PATHS[machine])
    diff_tool = get_path('scripts', 'diff.py')
    session.run('python', diff_tool, journal_filename, expected_journal)


@nox.session
def clean(session):
    """Clean up build files.

    Cleans up all artifacts that might get created during
    other ``nox`` sessions.

    There is no need for the session to create a ``virtualenv``
    here (we are just pretending to be ``make``).
    """
    # No need to create a virtualenv.
    session.virtualenv = False
    clean_dirs = (
        get_path('.cache'),
        get_path('.coverage'),
        get_path('.pytest_cache'),
        get_path('__pycache__'),
        get_path('benchmarks', 'memory', '__pycache__'),
        get_path('benchmarks', 'time', '__pycache__'),
        get_path('build'),
        get_path('dist'),
        get_path('docs', '__pycache__'),
        get_path('docs', 'build'),
        get_path('scripts', 'osx', 'dist_wheels'),
        get_path('scripts', 'osx', 'fixed_wheels'),
        get_path('src', 'bezier.egg-info'),
        get_path('src', 'bezier', '__pycache__'),
        get_path('src', 'bezier', 'extra-dll'),
        get_path('src', 'bezier', 'lib'),
        get_path('tests', '__pycache__'),
        get_path('tests', 'functional', '__pycache__'),
        get_path('tests', 'unit', '__pycache__'),
        get_path('wheelhouse'),
    )
    clean_globs = (
        get_path('.coverage'),
        get_path('*.mod'),
        get_path('*.pyc'),
        get_path('src', 'bezier', '*.pyc'),
        get_path('src', 'bezier', '*.pyd'),
        get_path('src', 'bezier', '*.so'),
        get_path('src', 'bezier', 'quadpack', '*.o'),
        get_path('src', 'bezier', '*.o'),
        get_path('tests', '*.pyc'),
        get_path('tests', 'functional', '*.pyc'),
        get_path('tests', 'unit', '*.pyc'),
    )
    for dir_path in clean_dirs:
        session.run(shutil.rmtree, dir_path, ignore_errors=True)
    for glob_path in clean_globs:
        for filename in glob.glob(glob_path):
            session.run(os.remove, filename)


@nox.session
def fortran_unit(session):
    session.interpreter = SINGLE_INTERP
    session.install('lcov_cobertura', 'pycobertura')
    if py.path.local.sysfind('make') is None:
        session.skip('`make` must be installed')
    if py.path.local.sysfind('gfortran') is None:
        session.skip('`gfortran` must be installed')
    if py.path.local.sysfind('lcov') is None:
        session.skip('`lcov` must be installed')
    test_dir = get_path('tests', 'fortran')
    lcov_filename = os.path.join(test_dir, 'coverage.info')
    session.chdir(test_dir)
    session.run('make', 'unit')
    session.chdir(NOX_DIR)
    session.run(
        'lcov',
        '--capture',
        '--directory',
        test_dir,
        '--output-file',
        lcov_filename,
    )
    session.run(
        'python',
        get_path('scripts', 'report_lcov.py'),
        '--lcov-filename',
        lcov_filename,
    )
    session.chdir(test_dir)
    session.run('make', 'clean')
