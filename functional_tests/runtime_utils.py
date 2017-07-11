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

"""Utilities for running functional tests as scripts."""


import argparse
import contextlib
import inspect
import io
import json
import os
import types

try:
    import matplotlib.pyplot as plt
except ImportError:
    plt = None
import numpy as np
try:
    import seaborn  # pylint: disable=unused-import
except ImportError:
    seaborn = None
import six

from bezier import _helpers


_FNL_TESTS_DIR = os.path.dirname(__file__)
_DOCS_DIR = os.path.abspath(
    os.path.join(_FNL_TESTS_DIR, '..', 'docs'))
IMAGES_DIR = os.path.join(_DOCS_DIR, 'images')


def _start_line(func):
    """Get the start line (in source) of a function.

    Args:
        func (~types.FunctionType): A Python function object.

    Returns:
        int: The start line (in source).
    """
    _, line = inspect.getsourcelines(func)
    return line


def get_parser():
    """Create a command line argument parser.

    Returns:
        argparse.ArgumentParser: An argument parser for functional tests.
    """
    parser = argparse.ArgumentParser(
        description='Run functional tests.')
    parser.add_argument('--save-plot', dest='save_plot',
                        action='store_true')
    return parser


def real_roots(coeffs):
    """Get real roots of a polynomial.

    Args:
        coeffs (List[Float]): List of polynomial coefficients.

    Returns:
        numpy.ndarray: The (sorted) real roots of the polynomial.
    """
    all_roots = np.roots(coeffs)
    filtered = all_roots[all_roots.imag == 0.0].real
    return np.sort(filtered)


def _convert_float(value):
    """Convert an "exact" value to a ``float``.

    Also works recursively if ``value`` is a list.

    Assumes a value is one of the following:

    * an integer
    * a string in C "%a" hex format for an IEEE-754 double precision number
    * a string fraction of the format "N/D"

    Args:
        value (Union[int, str, list]): Values to be converted.

    Returns:
        Union[float, list]: The converted value (or list of values).
    """
    if isinstance(value, list):
        return [_convert_float(element) for element in value]
    elif isinstance(value, six.integer_types):
        return float(value)
    elif value.startswith('0x') or value.startswith('-0x'):
        return float.fromhex(value)
    else:
        numerator, denominator = value.split('/')
        return float(numerator) / float(denominator)


@contextlib.contextmanager
def no_op_manager():
    """No-op context manager."""
    yield


def id_func(intersection_info):
    """Turn info from ``curve_intersections.json`` into a test ID.

    Args:
        intersection_info (dict): An intersection value loaded from
            ``curve_intersections.json``.

    Returns:
        str: An identifier formatted from the info.
    """
    return 'curves {:d} and {:d} (ID: {:d})'.format(
        intersection_info['curve1'], intersection_info['curve2'],
        intersection_info['id'])


def convert_floats(info, keys):
    """Modify ``info`` in-place to convert strings to floating point numbers.

    Args:
        info (List[dict, ...]): A list of dictionaries to be modified.
        keys (List[str, ...]): The keys within each dictionary that contain
            floating point values to be converted from a "custom" form
            to native Python ``float`` values.
    """
    for element in info:
        for key in keys:
            converted = _convert_float(element[key])
            if isinstance(converted, list):
                converted = np.asfortranarray(converted)

            element[key] = converted


def get_intersections_info():
    """Load curve and intersections info from JSON file.

    Returns:
        Tuple[List[dict, ...], List[dict, ...]]: The lists of

        * curve info dictionaries.
        * intersection info dictionaries.
    """
    filename = os.path.join(_FNL_TESTS_DIR, 'curves.json')
    with io.open(filename, 'r', encoding='utf-8') as file_obj:
        curves = json.load(file_obj)
    convert_floats(curves, keys=['control_points'])

    filename = os.path.join(_FNL_TESTS_DIR, 'curve_intersections.json')
    with io.open(filename, 'r', encoding='utf-8') as file_obj:
        intersections = json.load(file_obj)
    keys = ['intersections', 'curve1_params', 'curve2_params']
    convert_floats(intersections, keys=keys)

    return curves, intersections


class Config(object):
    """Run-time configuration.

    This is a mutable stand-in to allow test set-up to modify
    global state.
    """

    def __init__(self):
        self.running = False
        self.save_plot = False
        self.current_test = None
        self.parser = get_parser()
        self._wiggle = 8

    def _run(self, mod_globals):
        """Run all tests, in source order.

        Args:
            mod_globals (dict): The globals from a module.
        """
        found = []
        for key, value in six.iteritems(mod_globals):
            if not key.lower().startswith('test'):
                continue
            if not isinstance(value, types.FunctionType):
                continue
            found.append(value)

        found.sort(key=_start_line)

        for func in found:
            self.current_test = func.__name__
            func()

    @contextlib.contextmanager
    def wiggle(self, wiggle):
        """Make the context use a temporary wiggle room.

        Args:
            wiggle (int): The temporary amount of wiggle room.

        Yields:
            Config: the current configuration.
        """
        old_wiggle = self._wiggle
        try:
            self._wiggle = wiggle
            yield self
        finally:
            self._wiggle = old_wiggle

    def assert_close(self, approximated, exact):
        """Assert two values are close, with the local configuration in place.

        Args:
            approximated (float): The value that was computed.
            exact (float): The expected value.
        """
        msg = '{} ~= {} to {:d} bits'.format(
            approximated.hex(), exact.hex(), self._wiggle)
        assert _helpers.ulps_away(
            exact, approximated, num_bits=self._wiggle), msg

    def run(self, mod_globals):
        """Run all tests, in source order.

        Args:
            mod_globals (dict): The globals from a module.
        """
        running = self.running
        save_plot = self.save_plot
        try:
            self.running = True
            args = self.parser.parse_args()
            self.save_plot = args.save_plot
            self._run(mod_globals)
        finally:
            self.running = running
            self.save_plot = save_plot
            self.current_test = None

    def save_fig(self, extra=''):
        """Save the current figure.

        Uses the ``current_test`` for the filename and puts it
        in the ``${GIT_ROOT}/docs/images`` directory.

        Args:
            extra (Optional[str]): Extra information to put in the filename.
                Filename defaults to ``{current_test}.png`` but if ``extra``
                is passed, it will be ``{current_test}{extra}.png``.
        """
        filename = '{}{}.png'.format(self.current_test, extra)
        path = os.path.join(IMAGES_DIR, filename)
        plt.savefig(path, bbox_inches='tight')
        print('Saved {}'.format(filename))
