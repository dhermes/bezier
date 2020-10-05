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

"""Legacy members that have been replaced or renamed."""

import warnings

from bezier import triangle


DEPRECATION_MSG = (
    "`bezier.Surface` is deprecated and has been replaced by "
    "`bezier.Triangle`. It will be removed in a future release."
)


class Surface(triangle.Triangle):
    """Legacy alias for B |eacute| zier Surface class.

    This will be deprecated in a future release.
    """

    def __init__(self, *args, **kwargs):
        warnings.warn(DEPRECATION_MSG, DeprecationWarning)
        super().__init__(*args, **kwargs)
