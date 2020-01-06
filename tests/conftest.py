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

"""pytest shared testing configuration.

This

* Gets slow imports out of the way before running tests, so
  as not to have the cost of imports reflected in
  ``pytest --durations=N``.
"""

import pytest

try:
    from bezier import _HAS_SPEEDUP as HAS_SPEEDUP
except ImportError:  # pragma: NO COVER
    HAS_SPEEDUP = False


def pytest_addoption(parser):
    parser.addoption(
        "--ignore-slow",
        dest="ignore_slow",
        action="store_true",
        help="Ignore slow tests.",
    )


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow to run")


def pytest_collection_modifyitems(config, items):
    if not config.getoption("--ignore-slow"):
        return

    if HAS_SPEEDUP:
        return

    skip_slow = pytest.mark.skip(reason="--ignore-slow ignores the slow tests")
    for item in items:
        if "slow" in item.keywords:
            item.add_marker(skip_slow)
