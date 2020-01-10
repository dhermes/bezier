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
import pathlib
import shutil
import sys
import tempfile

import nox
import py.path


nox.options.error_on_external_run = True
nox.options.error_on_missing_interpreters = True

IS_MACOS = sys.platform == "darwin"
ON_APPVEYOR = os.environ.get("APPVEYOR") == "True"
DEPS = {
    "black": "black >= 19.10b0",
    "coverage": "coverage",
    "Cython": "Cython >= 0.29.14",
    "docutils": "docutils",
    "flake8": "flake8",
    "flake8-import-order": "flake8-import-order",
    "jsonschema": "jsonschema >= 3.2.0",
    "lcov_cobertura": "lcov_cobertura",
    "matplotlib": "matplotlib >= 3.1.2",
    "numpy": "numpy >= 1.18.1",
    "pycobertura": "pycobertura",
    "Pygments": "Pygments",
    "pylint": "pylint >= 2.4.4",
    "pytest": "pytest >= 5.3.2",
    "pytest-cov": "pytest-cov",
    "scipy": "scipy >= 1.4.1",
    "sympy": "sympy >= 1.5.1",
    "seaborn": "seaborn >= 0.9.0",
}
BASE_DEPS = (DEPS["numpy"], DEPS["pytest"])
NOX_DIR = os.path.abspath(os.path.dirname(__file__))
DOCS_DEPS = (
    "--requirement",
    os.path.join(NOX_DIR, "docs", "requirements.txt"),
)
DEFAULT_INTERPRETER = "3.8"
PYPY = "pypy3"
ALL_INTERPRETERS = ("3.6", "3.6-32", "3.7", "3.7-32", "3.8", "3.8-32", PYPY)
# Constants used for checking the journal of commands.
APPVEYOR = "appveyor"
CIRCLE_CI = "circleci"
TRAVIS_MACOS = "travis-macos"
JOURNAL_PATHS = {
    APPVEYOR: os.path.join("appveyor", "expected_journal.txt"),
    CIRCLE_CI: os.path.join(".circleci", "expected_journal.txt"),
    TRAVIS_MACOS: os.path.join("scripts", "macos", "travis_journal.txt"),
}


def get_path(*names):
    return os.path.join(NOX_DIR, *names)


def is_wheelhouse(directory):
    if directory is None:
        return False

    as_path = pathlib.Path(directory)
    if not as_path.is_dir():
        return False

    wheels = list(as_path.glob("*.whl"))
    # NOTE: This could also be done by using `next()` instead of `list()` and
    #       catching a `StopIteration` if empty.
    return len(wheels) > 0


def get_wheelhouse():
    candidates = (
        os.environ.get("WHEELHOUSE"),
        pathlib.Path.home() / "wheelhouse",
        pathlib.Path(os.path.abspath(os.sep)) / "wheelhouse",
    )

    for directory in candidates:
        if is_wheelhouse(directory):
            return directory

    return None


def pypy_setup(local_deps, session):
    wheelhouse = get_wheelhouse()
    if wheelhouse is not None:
        # Remove NumPy from dependencies.
        local_deps = list(local_deps)
        local_deps.remove(DEPS["numpy"])
        local_deps = tuple(local_deps)
        # Install NumPy and SciPy from pre-built wheels. Don't use ``DEPS``
        # to specify version range for NumPy and SciPy.
        session.install(
            "--no-index", "--find-links", wheelhouse, "numpy", "scipy"
        )
    return local_deps


def install_bezier(session, debug=False, env=None):
    if env is None:
        env = {}
    if debug:
        env["DEBUG"] = "True"
    session.install(".", env=env)


@nox.session(py=DEFAULT_INTERPRETER)
@nox.parametrize("check", [True, False])
def update_generated(session, check):
    # Install all dependencies.
    session.install(DEPS["Cython"])
    if check:
        command = get_path("scripts", "remove_cython_files.py")
        session.run("python", command)

    pyx_file = get_path("src", "python", "bezier", "_speedup.pyx")
    session.run("cython", pyx_file)

    command = get_path("scripts", "clean_cython.py")
    c_glob = get_path("src", "python", "bezier", "*.c")
    for c_source in glob.glob(c_glob):
        session.run(
            "python",
            command,
            "--filename",
            c_source,
            "--virtualenv-dirname",
            os.path.basename(session.virtualenv.location),
        )
    if check:
        command = get_path("scripts", "cython_update_check.py")
        session.run("python", command)


@nox.session(py=ALL_INTERPRETERS)
def unit(session):
    interpreter = session.virtualenv.interpreter
    unit_deps = BASE_DEPS + (DEPS["sympy"],)
    if interpreter == PYPY:
        local_deps = pypy_setup(unit_deps, session)
    else:
        local_deps = unit_deps + (DEPS["scipy"],)

    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    install_bezier(session, debug=True)
    # Run pytest against the unit tests.
    run_args = ["pytest"] + session.posargs + [get_path("tests", "unit")]
    session.run(*run_args)


@nox.session(py=DEFAULT_INTERPRETER)
def cover(session):
    # Install all test dependencies.
    local_deps = BASE_DEPS + (
        DEPS["scipy"],
        DEPS["sympy"],
        DEPS["pytest-cov"],
        DEPS["coverage"],
    )
    session.install(*local_deps)
    # Install this package.
    install_bezier(session, debug=True)
    # Run pytest with coverage against the unit tests.
    run_args = ["pytest", "--cov=bezier", "--cov=tests.unit"]
    run_args += session.posargs
    run_args += [get_path("tests", "unit")]
    session.run(*run_args)


@nox.session(py=ALL_INTERPRETERS)
def functional(session):
    interpreter = session.virtualenv.interpreter
    if interpreter == PYPY:
        local_deps = pypy_setup(BASE_DEPS, session)
    else:
        local_deps = BASE_DEPS

    # Install all test dependencies.
    session.install(*local_deps)
    # Install this package.
    install_bezier(session, debug=True)
    # Run pytest against the functional tests.
    run_args = ["pytest"] + session.posargs + [get_path("tests", "functional")]
    session.run(*run_args)


@nox.session(py=DEFAULT_INTERPRETER)
def docs(session):
    # Install all dependencies.
    session.install(*DOCS_DEPS)
    # Install this package.
    install_bezier(session, env={"BEZIER_NO_EXTENSION": "True"})
    # Run the script for building docs.
    command = get_path("scripts", "build_docs.sh")
    session.run(command, external=True)


def get_doctest_args(session):
    run_args = [
        "sphinx-build",
        "-W",
        "-b",
        "doctest",
        "-d",
        get_path("docs", "build", "doctrees"),
        get_path("docs"),
        get_path("docs", "build", "doctest"),
    ]
    run_args += session.posargs
    return run_args


@nox.session(py=DEFAULT_INTERPRETER)
def doctest(session):
    # Install all dependencies.
    session.install(DEPS["sympy"], *DOCS_DEPS)
    # Install this package.
    if IS_MACOS:
        command = get_path("scripts", "macos", "nox-install-for-doctest.sh")
        session.run(command, external=True)
    else:
        install_bezier(session)
    # Run the script for building docs and running doctests.
    run_args = get_doctest_args(session)
    session.run(*run_args)


@nox.session(py=DEFAULT_INTERPRETER)
def docs_images(session):
    # Install all dependencies.
    local_deps = DOCS_DEPS
    local_deps += (DEPS["matplotlib"], DEPS["seaborn"], DEPS["pytest"])
    session.install(*local_deps)
    # Install this package.
    install_bezier(session)
    # Use custom RC-file for matplotlib.
    env = {"MATPLOTLIBRC": "docs", "GENERATE_IMAGES": "True"}
    # Run the script for generating images for docs.
    run_args = get_doctest_args(session)
    session.run(*run_args, env=env)
    # Run the functional tests with --save-plot.
    fnl_tests_glob = get_path("tests", "functional", "test_*.py")
    modules_to_run = glob.glob(fnl_tests_glob)
    # Generate images for ``curve_intersections.json`` and
    # ``surface_intersections.json``.
    modules_to_run.extend(
        (
            get_path("tests", "functional", "make_curve_curve_images.py"),
            get_path("tests", "functional", "make_surface_surface_images.py"),
        )
    )
    # Make sure that the root directory is on the Python path so that
    # ``tests`` is import-able.
    env["PYTHONPATH"] = get_path()
    for filename in modules_to_run:
        session.run("python", filename, "--save-plot", env=env)


@nox.session(py=DEFAULT_INTERPRETER)
def lint(session):
    # Install all dependencies.
    local_deps = BASE_DEPS + (
        DEPS["black"],
        DEPS["docutils"],
        DEPS["flake8"],
        DEPS["flake8-import-order"],
        DEPS["matplotlib"],
        DEPS["Pygments"],
        DEPS["pylint"],
        DEPS["seaborn"],
    )
    session.install(*local_deps)
    # Install this package.
    install_bezier(session)
    # Run the script to check that the README and other docs are valid.
    check_path = get_path("scripts", "check_doc_templates.py")
    session.run("python", check_path)
    # Run the script to check that setup.py is valid.
    setup_file = get_path("setup.py")
    session.run(
        "python",
        setup_file,
        "check",
        "--metadata",
        "--restructuredtext",
        "--strict",
    )
    # Run `black --check` over all Python files
    check_black = get_path("scripts", "black_check_all_files.py")
    session.run("python", check_black)
    # Run flake8 over the code to check import order.
    session.run(
        "flake8",
        "--import-order-style=google",
        "--application-import-names=bezier,tests",
        get_path("src", "python", "bezier"),
        get_path("tests"),
    )
    # Run Pylint over the library source.
    session.run(
        "pylint",
        "--rcfile",
        "pylintrc",
        "--max-module-lines=3035",
        get_path("src", "python", "bezier"),
    )
    # Run Pylint over the tests source.
    session.run(
        "pylint",
        "--rcfile",
        "pylintrc",
        "--disable=missing-docstring",
        "--disable=protected-access",
        "--disable=too-many-public-methods",
        "--disable=import-outside-toplevel",
        "--disable=arguments-out-of-order",
        "--max-module-lines=2473",
        get_path("tests"),
    )

    if py.path.local.sysfind("clang-format") is not None:
        filenames = glob.glob(get_path("docs", "abi", "*.c"))
        filenames.append(get_path("src", "fortran", "include", "bezier.h"))
        filenames.extend(
            glob.glob(get_path("src", "fortran", "include", "bezier", "*.h"))
        )
        session.run(
            "clang-format", "-i", "-style=file", *filenames, external=True
        )


@nox.session(py=DEFAULT_INTERPRETER)
def blacken(session):
    """Run black code formatter."""
    session.install(DEPS["black"])
    check_black = get_path("scripts", "blacken_all_files.py")
    session.run("python", check_black)


@nox.session(py=DEFAULT_INTERPRETER)
@nox.parametrize("machine", [APPVEYOR, CIRCLE_CI, TRAVIS_MACOS])
def check_journal(session, machine):
    if machine == APPVEYOR and not ON_APPVEYOR:
        session.skip("Not currently running in AppVeyor.")
    if machine == CIRCLE_CI and os.environ.get("CIRCLECI") != "true":
        session.skip("Not currently running in CircleCI.")
    if machine == TRAVIS_MACOS:
        if os.environ.get("TRAVIS") != "true":
            session.skip("Not currently running in Travis.")
        if os.environ.get("TRAVIS_OS_NAME") != "osx":
            session.skip("Running in Travis, but not in a macOS job.")

    # Get a temporary file where the journal will be written.
    filehandle, journal_filename = tempfile.mkstemp(suffix="-journal.txt")
    os.close(filehandle)
    # Set the journal environment variable and install ``bezier``. Making sure
    # to limit to a single build job so commands are always in serial.
    session.install(DEPS["numpy"])  # Install requirement(s).
    env = {"BEZIER_JOURNAL": journal_filename, "NPY_NUM_BUILD_JOBS": "1"}
    install_bezier(session, env=env)
    # Compare the expected file to the actual results.
    session.run(
        "python",
        get_path("scripts", "post_process_journal.py"),
        "--journal-filename",
        journal_filename,
        "--machine",
        machine,
    )
    expected_journal = get_path(JOURNAL_PATHS[machine])
    diff_tool = get_path("scripts", "diff.py")
    session.run("python", diff_tool, journal_filename, expected_journal)


@nox.session(py=DEFAULT_INTERPRETER)
def fortran_unit(session):
    session.install(DEPS["lcov_cobertura"], DEPS["pycobertura"])
    if py.path.local.sysfind("make") is None:
        session.skip("`make` must be installed")
    if py.path.local.sysfind("gfortran") is None:
        session.skip("`gfortran` must be installed")
    if py.path.local.sysfind("lcov") is None:
        session.skip("`lcov` must be installed")
    test_dir = get_path("tests", "fortran")
    lcov_filename = os.path.join(test_dir, "coverage.info")
    session.chdir(test_dir)
    session.run("make", "unit", external=True)
    session.chdir(NOX_DIR)
    session.run(
        "lcov",
        "--capture",
        "--directory",
        test_dir,
        "--output-file",
        lcov_filename,
        external=True,
    )
    session.run(
        "python",
        get_path("scripts", "report_lcov.py"),
        "--lcov-filename",
        lcov_filename,
    )
    session.chdir(test_dir)
    session.run("make", "clean", external=True)


@nox.session(py=DEFAULT_INTERPRETER)
def validate_functional_test_cases(session):
    # Install all dependencies.
    session.install(DEPS["jsonschema"])

    session.run(
        "python", get_path("scripts", "validate_functional_test_cases.py")
    )


@nox.session(py=False)
def clean(session):
    """Clean up build files.

    Cleans up all artifacts that might get created during
    other ``nox`` sessions.

    There is no need for the session to create a ``virtualenv``
    here (we are just pretending to be ``make``).
    """
    clean_dirs = (
        get_path(".cache"),
        get_path(".coverage"),
        get_path(".pytest_cache"),
        get_path("__pycache__"),
        get_path("build"),
        get_path("dist"),
        get_path("docs", "__pycache__"),
        get_path("docs", "build"),
        get_path("scripts", "macos", "__pycache__"),
        get_path("scripts", "macos", "dist_wheels"),
        get_path("scripts", "macos", "fixed_wheels"),
        get_path("src", "python", "bezier.egg-info"),
        get_path("src", "python", "bezier", "__pycache__"),
        get_path("src", "python", "bezier", "extra-dll"),
        get_path("src", "python", "bezier", "lib"),
        get_path("tests", "__pycache__"),
        get_path("tests", "functional", "__pycache__"),
        get_path("tests", "unit", "__pycache__"),
        get_path("wheelhouse"),
    )
    clean_globs = (
        get_path(".coverage"),
        get_path("*.mod"),
        get_path("*.pyc"),
        get_path("src", "python", "bezier", "*.pyc"),
        get_path("src", "python", "bezier", "*.pyd"),
        get_path("src", "python", "bezier", "*.so"),
        get_path("src", "fortran", "quadpack", "*.o"),
        get_path("src", "fortran", "*.o"),
        get_path("tests", "*.pyc"),
        get_path("tests", "functional", "*.pyc"),
        get_path("tests", "unit", "*.pyc"),
    )
    for dir_path in clean_dirs:
        session.run(shutil.rmtree, dir_path, ignore_errors=True)
    for glob_path in clean_globs:
        for filename in glob.glob(glob_path):
            session.run(os.remove, filename)
