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

"""Utilities for running functional tests.

.. |eacute| unicode:: U+000E9 .. LATIN SMALL LETTER E WITH ACUTE
   :trim:
"""

import argparse
import contextlib
import enum
import inspect
import io
import json
import os
import types

import numpy as np
import six

import bezier
import bezier.curve
from tests import utils as base_utils

SPACING = np.spacing  # pylint: disable=no-member
FNL_TESTS_DIR = os.path.abspath(os.path.dirname(__file__))
_DOCS_DIR = os.path.abspath(os.path.join(FNL_TESTS_DIR, "..", "..", "docs"))
IMAGES_DIR = os.path.join(_DOCS_DIR, "images")


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
    parser = argparse.ArgumentParser(description="Run functional tests.")
    parser.add_argument("--save-plot", dest="save_plot", action="store_true")
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

    elif value.startswith("0x") or value.startswith("-0x"):
        return float.fromhex(value)

    else:
        numerator, denominator = value.split("/")
        return float(numerator) / float(denominator)


@contextlib.contextmanager
def no_op_manager():
    """No-op context manager."""
    yield


def curve_intersections_info():
    """Load curve and intersections info from JSON file.

    Returns:
        Tuple[Dict[str, dict], List[dict]]: The

        * mapping of curve info dictionaries.
        * list of intersection info dictionaries.
    """
    filename = os.path.join(FNL_TESTS_DIR, "curves.json")
    with io.open(filename, "r", encoding="utf-8") as file_obj:
        curve_json = json.load(file_obj)
    curves = {
        id_: CurveInfo.from_json(id_, info)
        for id_, info in six.iteritems(curve_json)
    }
    filename = os.path.join(FNL_TESTS_DIR, "curve_intersections.json")
    with io.open(filename, "r", encoding="utf-8") as file_obj:
        intersections_json = json.load(file_obj)
    intersections = [
        CurveIntersectionInfo.from_json(info, curves)
        for info in intersections_json
    ]
    return curves, intersections


def surface_intersections_info():
    """Load surface and intersections info from JSON file.

    Returns:
        Tuple[Dict[str, dict], List[dict]]: The

        * mapping of surface info dictionaries.
        * list of intersection info dictionaries.
    """
    filename = os.path.join(FNL_TESTS_DIR, "surfaces.json")
    with io.open(filename, "r", encoding="utf-8") as file_obj:
        surface_json = json.load(file_obj)
    surfaces = {
        id_: SurfaceInfo.from_json(id_, info)
        for id_, info in six.iteritems(surface_json)
    }
    filename = os.path.join(FNL_TESTS_DIR, "surface_intersections.json")
    with io.open(filename, "r", encoding="utf-8") as file_obj:
        intersections_json = json.load(file_obj)
    intersections = [
        SurfaceIntersectionsInfo.from_json(info, surfaces)
        for info in intersections_json
    ]
    return surfaces, intersections


def _ensure_empty(info):
    """Make sure a JSON info dictionary if empty.

    Args:
        info (dict): Expected to be exhausted.

    Raises:
        ValueError: If there are any keys remaining in ``info``.
    """
    # Make sure we've exhausted the data.
    if info:
        raise ValueError("Unexpected keys remaining in JSON info", info)


def id_func(value):
    """ID function for pytest parametrized tests.

    Args:
        value (Union[.IntersectionStrategy, CurveIntersectionInfo, \
            SurfaceIntersectionsInfo]: Either intersection info or an
            intersection strategy.

    Returns:
        str: The ID for a parameter in a parametrized test.
    """
    if isinstance(value, bezier.curve.IntersectionStrategy):
        return "strategy: {}".format(value.name)

    else:
        return value.test_id


def ulps_away(value1, value2, num_bits=1, eps=0.5 ** 40):
    r"""Determines if ``value1`` is within ``n`` ULPs of ``value2``.

    Uses ``np.spacing`` to determine the unit of least precision (ULP)
    for ``value1`` and then checks that the different between the values
    does not exceed ``n`` ULPs.

    When ``value1 == 0`` or ``value2 == 0``, we instead check that the other
    is exactly equal to ``0.0``.

    Args:
        value1 (float): The first value that being compared.
        value2 (float): The second value that being compared.
        num_bits (Optional[int]): The number of bits allowed to differ.
            Defaults to ``1``.
        eps (Optional[float]): The "zero threshold" to use when one of
            ``value1`` / ``value2`` is zero, but the other is not. Will
            only be used on 32-bit Linux.

    Returns:
        bool: Predicate indicating if the values agree to ``n`` bits.
    """
    if value1 == 0.0:
        if base_utils.IS_LINUX and not base_utils.IS_64_BIT:
            return abs(value2) < eps

        return value2 == 0.0

    if value2 == 0.0:
        if base_utils.IS_LINUX and not base_utils.IS_64_BIT:
            return abs(value1) < eps

        return value1 == 0.0

    local_epsilon = SPACING(value1)
    return abs(value1 - value2) <= num_bits * abs(local_epsilon)


class IncorrectCount(ValueError):
    """Custom exception for a "very bad" answer.

    This should be raised when the **computed** number of intersections
    disagrees with the actual number of intersections.
    """


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
            if not key.lower().startswith("test"):
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
        msg = "{} ~= {} to {:d} bits".format(
            approximated.hex(), exact.hex(), self._wiggle
        )
        assert ulps_away(exact, approximated, num_bits=self._wiggle), msg

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

    def save_fig(self, extra=""):
        """Save the current figure.

        Uses the ``current_test`` for the filename and puts it
        in the ``${GIT_ROOT}/docs/images`` directory.

        Args:
            extra (Optional[str]): Extra information to put in the filename.
                Filename defaults to ``{current_test}.png`` but if ``extra``
                is passed, it will be ``{current_test}{extra}.png``.
        """
        # NOTE: We import the plotting library at runtime to
        #       avoid the cost for users that only want to compute.
        #       The ``matplotlib`` import is a tad expensive.
        import matplotlib.pyplot as plt

        filename = "{}{}.png".format(self.current_test, extra)
        path = os.path.join(IMAGES_DIR, filename)
        plt.savefig(path, bbox_inches="tight")
        print("Saved {}".format(filename))


class CurveInfo(object):  # pylint: disable=too-few-public-methods
    """Information about a curve from ``curves.json``.

    The ``curves.json`` file contains a dictionary where each key is the ID
    of the given curve and each value is a curve. The curves are described
    by the JSON-schema in ``tests/functional/schema/curve.json``.

    Args:
        id_ (str): The ID of the curve.
        control_points (numpy.ndarray): The control points.
        implicitized (Optional[List[List[int]]]): The coefficient triples
            defining the algebraic curve that contains the
            B |eacute| zier curve.
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

        Returns:
            CurveInfo: The curve info parsed from the JSON.
        """
        control_points = info.pop("control_points")
        control_points = np.asfortranarray(_convert_float(control_points))
        implicitized = info.pop("implicitized", None)
        # Optional fields.
        note = info.pop("note", None)
        _ensure_empty(info)
        return cls(id_, control_points, implicitized=implicitized, note=note)


class CurveIntersectionType(enum.Enum):
    """Enum describing curve intersection.

    These values correspond to the ``type`` enum property in the
    ``curve_intersection.json`` JSON-schema.
    """

    coincident = "coincident"
    """Curves lie on the same underlying algebraic curve."""
    no_intersection = "no-intersection"
    """Curves do not intersect."""
    tangent = "tangent"
    """Curves are tangent at the point of intersection."""
    standard = "standard"
    """Intersection is not **any** of the other types."""


# pylint: disable=too-many-instance-attributes


class CurveIntersectionInfo(object):
    """Information about an intersection from ``curve_intersections.json``.

    The ``curve_intersections.json`` file contains a list of intersection
    cases. The intersection cases are described by the JSON-schema in
    ``tests/functional/schema/curve_intersection.json``.

    Args:
        id_ (int): The intersection ID.
        curve1_info (CurveInfo): The curve information for the first curve in
            the intersection.
        curve2_info (CurveInfo): The curve information for the second curve in
            the intersection.
        type_ (CurveIntersectionType): Describes how the curves intersect.
        intersections (numpy.ndarray): ``Nx2`` array of ``x-y`` coordinate
            pairs of intersection points.
        curve1_params (numpy.ndarray): 1D array, the parameters along
            ``curve1`` where the intersections occur (in the same order as
            the rows of ``intersections``). These are typically called
            ``s``-parameters.
        curve2_params (numpy.ndarray): 1D array, the parameters along
            ``curve2`` where the intersections occur (in the same order as
            the rows of ``intersections``). These are typically called
            ``t``-parameters.
        curve1_polys (Optional[List[List[int]]]): The coefficients of the
            polynomials that determine the values in ``curve1_params``.
        curve2_polys (Optional[List[List[int]]]): The coefficients of the
            polynomials that determine the values in ``curve2_params``.
        note (Optional[str]): A note about the intersection (e.g. why it is
            unique / problematic).
    """

    # pylint: disable=too-many-arguments

    def __init__(
        self,
        id_,
        curve1_info,
        curve2_info,
        type_,
        intersections,
        curve1_params,
        curve2_params,
        curve1_polys=None,
        curve2_polys=None,
        note=None,
    ):
        self.id_ = id_
        self.intersections = intersections
        self.type_ = type_
        self.curve1_info = curve1_info
        self.curve1_params = curve1_params
        self.curve1_polys = curve1_polys
        self.curve2_info = curve2_info
        self.curve2_params = curve2_params
        self.curve2_polys = curve2_polys
        self.note = note
        self.num_params = self._verify_dimensions()
        self._verify_data()

    # pylint: enable=too-many-arguments

    def _verify_dimensions(self):
        """Verify that all the dimensions are the same.

        Returns:
            int: The number of parameters / intersections expected.

        Raises:
            ValueError: If one of the values is not the "expected" shape.
        """
        if self.curve1_params.ndim != 1:
            raise ValueError(
                "Expected 1-dimensional data for ``curve1_params``."
            )

        # Unpack into one value now that we know 1D.
        num_params, = self.curve1_params.shape
        shape = (num_params,)
        if self.curve2_params.shape != shape:
            msg = "Expected shape {} for ``curve2_params``.".format(shape)
            raise ValueError(msg)

        shape = (2, num_params)
        if self.intersections.shape != shape:
            msg = "Expected shape {} for ``intersections``.".format(shape)
            raise ValueError(msg)

        if self.curve1_polys is not None:
            if len(self.curve1_polys) != num_params:
                raise ValueError(
                    "Unexpected number of ``curve1_polys``",
                    len(self.curve1_polys),
                    "Expected",
                    num_params,
                )

        if self.curve2_polys is not None:
            if len(self.curve2_polys) != num_params:
                raise ValueError(
                    "Unexpected number of ``curve2_polys``",
                    len(self.curve2_polys),
                    "Expected",
                    num_params,
                )

        return num_params

    def _verify_data(self):
        """Verify assumptions about the data.

        * The intersections are sorted by s-value.

        Raises:
            ValueError: If the assumptions are not met.
        """
        sorted_s = np.sort(self.curve1_params)
        if not np.all(sorted_s == self.curve1_params):
            raise ValueError(
                "Expected s-parameters (``curve1_params``) to be "
                "in ascending order."
            )

    @property
    def test_id(self):
        """str: The ID for this intersection in functional tests."""
        return "curves {!r} and {!r} (ID: {:d})".format(
            self.curve1_info.id_, self.curve2_info.id_, self.id_
        )

    @property
    def img_filename(self):
        """str: Filename to use when saving images for this intersection."""
        return "curves{}_and_{}".format(
            self.curve1_info.id_, self.curve2_info.id_
        )

    @property
    def params(self):
        """Get the parameters and intersection points.

        No verification is done here, rather, it's done in the
        constructor / in :meth:`from_json`.

        Returns:
            Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]: The

            * ``s``-parameters (1D ``N`` array)
            * ``t``-parameters (1D ``N`` array)
            * intersection points (2D ``N x 2`` array)
        """
        return self.curve1_params, self.curve2_params, self.intersections

    # pylint: disable=missing-return-type-doc

    @property
    def curve1(self):
        """The first B |eacute| zier curve in the intersection.

        Returns:
            ~bezier.curve.Curve: The first B |eacute| zier curve.
        """
        return self.curve1_info.curve

    @property
    def curve2(self):
        """The second B |eacute| zier curve in the intersection.

        Returns:
            ~bezier.curve.Curve: The second B |eacute| zier curve.
        """
        return self.curve2_info.curve

    # pylint: enable=missing-return-type-doc

    @property
    def nodes1(self):
        """The first B |eacute| zier nodes in the intersection.

        Returns:
            numpy.ndarray: The first B |eacute| zier nodes.
        """
        return self.curve1_info.control_points

    @property
    def nodes2(self):
        """The second B |eacute| zier nodes in the intersection.

        Returns:
            numpy.ndarray: The second B |eacute| zier nodes.
        """
        return self.curve2_info.control_points

    @classmethod
    def from_json(cls, info, curves):
        """Convert JSON curve intersection info into ``CurveIntersectionInfo``.

        This involves parsing the dictionary and converting some stringified
        values (rationals and IEEE-754) to Python ``float``-s.

        Args:
            info (dict): The JSON data of the curve intersection.
            curves (Dict[str, CurveInfo]): An already parsed dictionary of
                curve information.

        Returns:
            CurveIntersectionInfo: The intersection info parsed from the JSON.
        """
        id_ = info.pop("id")
        curve1_info = curves[info.pop("curve1")]
        curve2_info = curves[info.pop("curve2")]
        type_ = CurveIntersectionType(info.pop("type"))
        intersections = np.asfortranarray(
            _convert_float(info.pop("intersections"))
        )
        curve1_params = np.asfortranarray(
            _convert_float(info.pop("curve1_params"))
        )
        curve2_params = np.asfortranarray(
            _convert_float(info.pop("curve2_params"))
        )
        if intersections.size == 0:
            intersections = intersections.reshape((2, 0), order="F")
        # Optional fields.
        curve1_polys = info.pop("curve1_polys", None)
        curve2_polys = info.pop("curve2_polys", None)
        note = info.pop("note", None)
        _ensure_empty(info)
        return cls(
            id_,
            curve1_info,
            curve2_info,
            type_,
            intersections,
            curve1_params,
            curve2_params,
            curve1_polys=curve1_polys,
            curve2_polys=curve2_polys,
            note=note,
        )


# pylint: enable=too-many-instance-attributes


class SurfaceInfo(object):  # pylint: disable=too-few-public-methods
    """Information about a surface from ``surfaces.json``.

    The ``surfaces.json`` file contains a dictionary where each key is the ID
    of the given surface and each value is a surface. The surfaces are
    described by the JSON-schema in ``tests/functional/schema/surface.json``.

    Args:
        id_ (str): The ID of the surface.
        control_points (numpy.ndarray): The control points.
        note (Optional[str]): A note about the surface (e.g. what is it
            related to).
    """

    def __init__(self, id_, control_points, note=None):
        self.id_ = id_
        self.control_points = control_points
        self.surface = bezier.Surface.from_nodes(control_points, _copy=False)
        self.note = note

    @classmethod
    def from_json(cls, id_, info):
        """Convert JSON surface info into ``SurfaceInfo``.

        This involves parsing the dictionary and converting some stringified
        values (rationals and IEEE-754) to Python ``float``-s.

        Args:
            id_ (str): The ID of the surface.
            info (dict): The JSON data of the surface.

        Returns:
            SurfaceInfo: The surface info parsed from the JSON.
        """
        control_points = info.pop("control_points")
        control_points = np.asfortranarray(_convert_float(control_points))
        # Optional fields.
        note = info.pop("note", None)
        _ensure_empty(info)
        return cls(id_, control_points, note=note)


# pylint: disable=too-few-public-methods


class SurfaceIntersectionInfo(object):
    """Basic wrapper indicating an intersection is one of the two surfaces.

    Args:
        first (bool): Flag indicating if the first or second surface in an
            intersection is contained in the other.
    """

    def __init__(self, first):
        self.first = first
        # Will be set later by the `SurfaceIntersectionsInfo` constructor.
        self.parent = None


# pylint: enable=too-few-public-methods


# pylint: disable=too-few-public-methods


class CurvedPolygonInfo(object):
    """Information about a single curved polygon intersection.

    The ``surface_intersections.json`` file contains surface-surface
    intersections, each of which may contain multiple curved polygons
    as the intersected area (e.g. the "6Q"-"7Q" intersection splits into two
    disjoint regions). Such a curved polygon is described by the JSON-schema
    in ``tests/functional/schema/curved_polygon.json``.

    Args:
        nodes (Optional[numpy.ndarray]): ``Nx2`` array of ``x-y``
            coordinate pairs of intersection points.
        edge_list (List[int]): List of edge indices among
            ``{0, 1, 2, 3, 4, 5}``. If the index is in ``{0, 1, 2}``, the
            intersection occurs on the first surface and on the edge
            index (one of 0, 1 or 2) **of that surface**. If the index is in
            ``{3, 4, 5}``, the intersection is on the second surface.
            An intersection requires **two** edges, but this one is the edge
            that the boundary (of the intersection) continues along.
        start_params (numpy.ndarray): 1D array, the parameters along
            ``surface1`` where the intersections occur (in the same order as
            the rows of ``intersections``). These are typically called
            ``s``-parameters.
        end_params (numpy.ndarray): 1D array, the parameters along
            ``surface2`` where the intersections occur (in the same order as
            the rows of ``intersections``). These are typically called
            ``t``-parameters.
        start_param_polys (Optional[List[List[int]]]): The coefficients of the
            polynomials that determine the values in ``start_params``.
        end_param_polys (Optional[List[List[int]]]): The coefficients of the
            polynomials that determine the values in ``end_params``.
    """

    # pylint: disable=too-many-arguments

    def __init__(
        self,
        nodes,
        edge_list,
        start_params,
        end_params,
        start_param_polys=None,
        end_param_polys=None,
    ):
        # Will be set later by the `SurfaceIntersectionsInfo` constructor.
        self.parent = None
        self.nodes = nodes
        self.edge_list = edge_list
        self.start_params = start_params
        self.start_param_polys = start_param_polys
        self.end_params = end_params
        self.end_param_polys = end_param_polys
        self.num_nodes = self._verify()

    # pylint: enable=too-many-arguments

    def _verify_polynomials(self, num_nodes, start=True):
        """Verify a list of polynomial coefficients.

        Args:
            num_nodes (int): The (expected) number of nodes.
            start (Optional[bool]): Indicates if start or end polynomials
                should be verified.

        Raises:
            ValueError: If the sizes are not as expected.
            ValueError: If the types are not as expected.
        """
        if start:
            polynomials = self.start_param_polys
            name = "start"
        else:
            polynomials = self.end_param_polys
            name = "end"
        if polynomials is None:
            return

        err_msg1 = "Unexpected number of {} parameter polynomials.".format(
            name
        )
        template2 = "{} parameter polynomial should have integer coefficients."
        err_msg2 = template2.format(name)
        if len(polynomials) != num_nodes:
            raise ValueError(err_msg1)

        for polynomial in polynomials:
            if not all(isinstance(coeff, int) for coeff in polynomial):
                raise ValueError(err_msg2, polynomial)

    def _verify(self):
        """Verify the state.

        Returns:
            int: The number of nodes.

        Raises:
            ValueError: If the sizes are not as expected.
            ValueError: If the types are not as expected.
        """
        num_nodes, cols = self.nodes.shape
        if cols != 2:
            raise ValueError("Unexpected shape of nodes")

        np_edges = np.asfortranarray(self.edge_list, dtype=int)
        if np_edges.shape != (num_nodes,):
            raise ValueError("Unexpected shape of edge list")

        if not np.all(np_edges == self.edge_list):
            raise ValueError("Edge pairs were expected to be integers.")

        if self.start_params.shape != (num_nodes,):
            raise ValueError("Unexpected shape of start parameters")

        if self.end_params.shape != (num_nodes,):
            raise ValueError("Unexpected shape of end parameters")

        self._verify_polynomials(num_nodes, start=True)
        self._verify_polynomials(num_nodes, start=False)
        return num_nodes

    @classmethod
    def from_json(cls, info):
        """Parse and convert JSON curved polygon info.

        This is a curved polygon arising from surface-surface intersection,
        so the data is tailored to how each edge of the curved polygon relates
        to the original two surfaces.

        This involves parsing the dictionary and converting some stringified
        values (rationals and IEEE-754) to Python ``float``-s.

        Args:
            info (Union[dict, bool]): The JSON data of the curved polygon or
                a boolean indicating if the first or second surface is the
                intersection (i.e. one is fully contained in the other).

        Returns:
            Union[SurfaceIntersectionInfo, CurvedPolygonInfo]: A basic
            object containing surface information (if one of the surfaces is
            contained in the other) or the curved polygon info parsed from
            the JSON.
        """
        if isinstance(info, bool):
            return SurfaceIntersectionInfo(info)

        else:
            nodes = np.asfortranarray(_convert_float(info.pop("nodes")))
            if nodes.size == 0:
                nodes = nodes.reshape((2, 0), order="F")
            start_params = np.asfortranarray(
                _convert_float(info.pop("start_params"))
            )
            end_params = np.asfortranarray(
                _convert_float(info.pop("end_params"))
            )
            edge_list = info.pop("edge_list")
            # Optional fields.
            start_param_polys = info.pop("start_param_polys", None)
            end_param_polys = info.pop("end_param_polys", None)
            _ensure_empty(info)
            return cls(
                nodes,
                edge_list,
                start_params,
                end_params,
                start_param_polys=start_param_polys,
                end_param_polys=end_param_polys,
            )


# pylint: enable=too-few-public-methods


class SurfaceIntersectionsInfo(object):
    """Information about an intersection from ``surface_intersections.json``.

    The ``surface_intersections.json`` file contains a list of intersection
    cases. The intersection cases are described by the JSON-schema in
    ``tests/functional/schema/surface_intersection.json``.

    Args:
        id_ (int): The intersection ID.
        surface1_info (SurfaceInfo): The surface information for the first
            surface in the intersection.
        surface2_info (SurfaceInfo): The surface information for the second
            surface in the intersection.
        intersections (List[Union[SurfaceIntersectionInfo, \
            CurvedPolygonInfo]]): The info for each intersection region.
        note (Optional[str]): A note about the intersection(s) (e.g. why
            unique / problematic).
    """

    def __init__(
        self, id_, surface1_info, surface2_info, intersections, note=None
    ):
        self.id_ = id_
        self.surface1_info = surface1_info
        self.surface2_info = surface2_info
        self.intersections = intersections
        self.note = note
        self._set_parent()

    def _set_parent(self):
        """Set the current instance as parent in each child."""
        for info in self.intersections:
            info.parent = self

    @property
    def test_id(self):
        """str: The ID for this intersection in functional tests."""
        return "surfaces {!r} and {!r} (ID: {:d})".format(
            self.surface1_info.id_, self.surface2_info.id_, self.id_
        )

    @property
    def img_filename(self):
        """str: Filename to use when saving images for this intersection."""
        return "surfaces{}_and_{}".format(
            self.surface1_info.id_, self.surface2_info.id_
        )

    # pylint: disable=missing-return-type-doc

    @property
    def surface1(self):
        """The first B |eacute| zier surface in the intersection.

        Returns:
            ~bezier.surface.Surface: The first B |eacute| zier surface.
        """
        return self.surface1_info.surface

    @property
    def surface2(self):
        """The second B |eacute| zier surface in the intersection.

        Returns:
            ~bezier.surface.Surface: The second B |eacute| zier surface.
        """
        return self.surface2_info.surface

    # pylint: enable=missing-return-type-doc

    @classmethod
    def from_json(cls, info, surfaces):
        """Parse and convert JSON surface intersection info.

        This involves parsing the dictionary and converting some stringified
        values (rationals and IEEE-754) to Python ``float``-s.

        Args:
            info (dict): The JSON data of the surface intersection.
            surfaces (Dict[str, SurfaceInfo]): An already parsed dictionary of
                surface information.

        Returns:
            SurfaceIntersectionsInfo: The intersection info parsed from
                the JSON.
        """
        id_ = info.pop("id")
        surface1_info = surfaces[info.pop("surface1")]
        surface2_info = surfaces[info.pop("surface2")]
        intersections = info.pop("intersections")
        # Optional fields.
        note = info.pop("note", None)
        _ensure_empty(info)
        # Convert ``intersections`` JSON to ``CurvedPolygonInfo``.
        intersections = [
            CurvedPolygonInfo.from_json(curved_polygon)
            for curved_polygon in intersections
        ]
        return cls(id_, surface1_info, surface2_info, intersections, note=note)
