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

# NOTE: This module is written in the ``pytest`` style, i.e. it does not use
#       ``unittest``. The other test modules in this project use ``unittest``,
#       they were written before ``pytest`` was used in this project (the
#       original test runner was ``nose``).

import unittest.mock
import warnings

import numpy as np


class TestSurface:
    @staticmethod
    def _get_target_class():
        from bezier import _legacy

        return _legacy.Surface

    def _make_one(self, *args, **kwargs):
        klass = self._get_target_class()
        return klass(*args, **kwargs)

    def test_warn(self, monkeypatch):
        from bezier import _legacy

        warn_mock = unittest.mock.Mock(spec=())
        monkeypatch.setattr(warnings, "warn", warn_mock)

        nodes = np.empty((1, 3), order="F")
        self._make_one(nodes, 1, copy=False, verify=False)

        warn_mock.assert_called_once_with(
            _legacy.DEPRECATION_MSG, DeprecationWarning
        )
