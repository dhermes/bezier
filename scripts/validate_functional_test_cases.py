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

import json
import os
import pathlib
import sys

import jsonschema


HERE = pathlib.Path(__file__).resolve().parent
REPO_ROOT = HERE.parent
DATA_DIR = REPO_ROOT / "tests" / "functional"
SCHEMA_DIR = DATA_DIR / "schema"


def _verify_map(map_filename, schema_filename, name):
    """Verifies a map with string keys and values of a given schema.

    Args:
        map_filename (pathlib.Path): The filename containing the map to verify.
        schema_filename (pathlib.Path): The filename containing the schema for
            the values in the map
        name (str): The name of the schema.

    Returns:
        bool: Indicated if there were any failures.
    """
    with open(map_filename, "r") as file_obj:
        object_map = json.load(file_obj)

    with open(schema_filename, "r") as file_obj:
        schema = json.load(file_obj)

    # NOTE: We need to set a custom resolver for ``$ref`` to the local
    #       filesystem.
    #       See: https://github.com/Julian/jsonschema/issues/313
    resolver = jsonschema.RefResolver(
        base_uri=f"file://{SCHEMA_DIR}{os.path.sep}", referrer=schema
    )
    failed = False
    for object_id, info in object_map.items():
        try:
            jsonschema.validate(info, schema, resolver=resolver)
        except jsonschema.ValidationError:
            print(f"{name} {object_id} does not adhere to the schema.")
            failed = True

    return failed


def verify_curves(exit_status):
    """Verify that all curves adhere to the curve schema.

    Args:
        exit_status (int): The already set exit status.

    Returns:
        int: The newly updated exit status (will add a flag in the bit
            corresponding to 2**0 == 1).
    """
    curves_file = DATA_DIR / "curves.json"
    curve_schema_file = SCHEMA_DIR / "curve.json"
    failed = _verify_map(curves_file, curve_schema_file, "Curve")
    if failed:
        return exit_status | 1

    return exit_status


def verify_triangles(exit_status):
    """Verify that all triangles adhere to the triangle schema.

    Args:
        exit_status (int): The already set exit status.

    Returns:
        int: The newly updated exit status (will add a flag in the bit
            corresponding to 2**1 == 2).
    """
    triangles_file = DATA_DIR / "triangles.json"
    triangle_schema_file = SCHEMA_DIR / "triangle.json"
    failed = _verify_map(triangles_file, triangle_schema_file, "Triangle")
    if failed:
        return exit_status | 2

    return exit_status


def _verify_list(list_filename, schema_filename, name):
    """Verifies a list with elements of a given schema.

    Args:
        list_filename (pathlib.Path): The filename containing the list to
            verify.
        schema_filename (pathlib.Path): The filename containing the schema for
            the elements in the list.
        name (str): The name of the schema.

    Returns:
        bool: Indicated if there were any failures.
    """
    with open(list_filename, "r") as file_obj:
        elements = json.load(file_obj)

    with open(schema_filename, "r") as file_obj:
        schema = json.load(file_obj)

    # NOTE: We need to set a custom resolver for ``$ref`` to the local
    #       filesystem.
    #       See: https://github.com/Julian/jsonschema/issues/313
    resolver = jsonschema.RefResolver(
        base_uri=f"file://{SCHEMA_DIR}{os.path.sep}", referrer=schema
    )
    failed = False
    for element in elements:
        id_ = element.get("id", "<unknown>")
        try:
            jsonschema.validate(element, schema, resolver=resolver)
        except jsonschema.ValidationError:
            print(f"{name} {id_} does not adhere to the schema.")
            failed = True

    return failed


def verify_curve_intersections(exit_status):
    """Verify that all curve intersections adhere to the schema.

    Args:
        exit_status (int): The already set exit status.

    Returns:
        int: The newly updated exit status (will add a flag in the bit
            corresponding to 2**2 == 4).
    """
    intersections_file = DATA_DIR / "curve_intersections.json"
    intersections_schema_file = SCHEMA_DIR / "curve_intersection.json"
    failed = _verify_list(
        intersections_file, intersections_schema_file, "Curve Intersection"
    )
    if failed:
        return exit_status | 4

    return exit_status


def verify_triangle_intersections(exit_status):
    """Verify that all triangles intersections adhere to the schema.

    Args:
        exit_status (int): The already set exit status.

    Returns:
        int: The newly updated exit status (will add a flag in the bit
            corresponding to 2**3 == 8).
    """
    intersections_file = DATA_DIR / "triangle_intersections.json"
    intersections_schema_file = SCHEMA_DIR / "triangle_intersection.json"
    failed = _verify_list(
        intersections_file, intersections_schema_file, "Triangle Intersection"
    )
    if failed:
        return exit_status | 8

    return exit_status


def main():
    """Main entrypoint for this script."""
    exit_status = verify_curves(0)
    exit_status = verify_curve_intersections(exit_status)
    exit_status = verify_triangles(exit_status)
    exit_status = verify_triangle_intersections(exit_status)
    sys.exit(exit_status)


if __name__ == "__main__":
    main()
