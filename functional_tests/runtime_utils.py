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

import bezier
from bezier import _helpers


if seaborn is not None:
    seaborn.set()  # Required in `seaborn >= 0.8`
FNL_TESTS_DIR = os.path.abspath(os.path.dirname(__file__))
_DOCS_DIR = os.path.abspath(
    os.path.join(FNL_TESTS_DIR, '..', 'docs'))
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

    * :data:`None`
    * an integer
    * a string in C "%a" hex format for an IEEE-754 double precision number
    * a string fraction of the format "N/D"
    * a list of one of the accepted types (incl. a list)

    Args:
        value (Union[int, str, list]): Values to be converted.

    Returns:
        Union[float, list]: The converted value (or list of values).
    """
    if value is None:
        return None
    elif isinstance(value, list):
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


def curve_id_func(intersection_info):
    """Turn info from ``curve_intersections.json`` into a test ID.

    Args:
        intersection_info (dict): An intersection value loaded from
            ``curve_intersections.json``.

    Returns:
        str: An identifier formatted from the info.
    """
    return 'curves {!r} and {!r} (ID: {:d})'.format(
        intersection_info['curve1'], intersection_info['curve2'],
        intersection_info['id'])


def surface_id_func(intersection_info):
    """Turn info from ``surface_intersections.json`` into a test ID.

    Args:
        intersection_info (dict): An intersection value loaded from
            ``surface_intersections.json``.

    Returns:
        str: An identifier formatted from the info.
    """
    return 'surfaces {!r} and {!r} (ID: {:d})'.format(
        intersection_info['surface1'], intersection_info['surface2'],
        intersection_info['id'])


def convert_floats(info, keys):
    """Modify ``info`` in-place to convert strings to floating point numbers.

    Args:
        info (List[dict]): A list of dictionaries to be modified.
        keys (List[str]): The keys within each dictionary that contain
            floating point values to be converted from a "custom" form
            to native Python ``float`` values.
    """
    for element in info:
        for key in keys:
            if key not in element:
                continue

            converted = _convert_float(element[key])
            if isinstance(converted, list):
                converted = np.asfortranarray(converted)

            element[key] = converted


def curve_intersections_info():
    """Load curve and intersections info from JSON file.

    Returns:
        Tuple[Dict[str, dict], List[dict]]: The

        * mapping of curve info dictionaries.
        * list of intersection info dictionaries.
    """
    filename = os.path.join(FNL_TESTS_DIR, 'curves.json')
    with io.open(filename, 'r', encoding='utf-8') as file_obj:
        curve_json = json.load(file_obj)
    curves = {id_: CurveInfo.from_json(id_, info)
              for id_, info in six.iteritems(curve_json)}

    filename = os.path.join(FNL_TESTS_DIR, 'curve_intersections.json')
    with io.open(filename, 'r', encoding='utf-8') as file_obj:
        intersections = json.load(file_obj)
    keys = ['intersections', 'curve1_params', 'curve2_params']
    convert_floats(intersections, keys=keys)

    return curves, intersections


def surface_intersections_info():
    """Load surface and intersections info from JSON file.

    Returns:
        Tuple[Dict[str, dict], List[dict]]: The

        * mapping of surface info dictionaries.
        * list of intersection info dictionaries.
    """
    filename = os.path.join(FNL_TESTS_DIR, 'surfaces.json')
    with io.open(filename, 'r', encoding='utf-8') as file_obj:
        surfaces = json.load(file_obj)
    convert_floats(six.itervalues(surfaces), keys=['control_points'])

    filename = os.path.join(FNL_TESTS_DIR, 'surface_intersections.json')
    with io.open(filename, 'r', encoding='utf-8') as file_obj:
        intersections = json.load(file_obj)
    keys = ['intersections', 'start_params', 'end_params']
    convert_floats(intersections, keys=keys)

    return surfaces, intersections


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


class CurveInfo(object):
    r"""Information about a curve from ``curves.json``.

    These are expected to have three keys:

    * ``control_points``: A list of ``x-y`` coordinates of the control points
      in the curve. The coordinates themselves can be integers, stringified
      fractions or stringified IEEE-754 values (``%a`` format).
    * ``note`` (optional): Description of the curve / curve segment.
    * ``implicitized`` (optional): The algebraic curve that contains
      this curve as a segment. (Only provided if the curve comes from
      rational control points.) For example, for "curve 2"

      .. math::

         x(s) = \frac{9 - 8 s}{8}
         y(s) = \frac{(2 s - 1)^2}{2}

      so we have

      .. math::

         8 x = 9 - 8 s \Longrightarrow 8 s - 4 = 5 - 8 x
         2 y = (2 s - 1)^2 \Longrightarrow 32 y = (8 s - 4)^2
         \Longrightarrow 32 y = (5 - 8 x)^2
         \Longrightarrow 25 - 32 y - 80 x + 64 x^2 = 0

      and this implicitized algebraic curve corresponds to

      .. code-block:: python

         [
             [ 25, 0, 0],
             [-32, 0, 1],
             [-80, 1, 0],
             [ 64, 2, 0],
         ]

    This representation is a list of triples of integers
    ``coefficient, degree_x, degree_y`` where ``degree_x`` and ``degree_y``
    are non-negative. Though it's not required, we also have
    ``coefficient != 0`` and the triples are sorted by total degree
    (``degree_x + degree_y``) and then by ``degree_x`` to break ties.

    In addition, each curve comes with an ID from a dictionary, i.e.
    ``curves.json`` uses ID keys to identify the curves, rather than just
    having a list of curve info.

    Args:
        id_ (str): The ID of the curve.
        control_points (numpy.ndarray): The control points.
        implicitized (Optional[List[List[int]]]): The coefficient triples
            defining the algebraic curve that contains the curve.
        note (Optional[str]): A note about the curve (e.g. what is it
            related to).
    """

    def __init__(self, id_, control_points, implicitized=None, note=None):
        self.id_ = id_
        self.control_points = control_points
        self.curve = bezier.Curve.from_nodes(control_points, _copy=False)
        self.implicitized = implicitized
        self.note = note

    @classmethod
    def from_json(cls, id_, info):
        """Convert JSON curve info into ``CurveInfo``.

        This involves parsing the dictionary and converting some stringified
        values (rationals and IEEE-754) to Python ``float``-s.

        Args:
            id_ (str): The ID of the curve.
            info (dict): The JSON data of the curve.
        """
        control_points = info.pop('control_points')
        control_points = np.asfortranarray(_convert_float(control_points))
        implicitized = info.pop('implicitized', None)
        note = info.pop('note', None)
        return cls(id_, control_points, implicitized=implicitized, note=note)
