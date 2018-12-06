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

import os
import platform
import sys

# See: https://docs.python.org/3/library/platform.html#cross-platform
if sys.maxsize == 2 ** 63 - 1:
    IS_64_BIT = True
elif sys.maxsize == 2 ** 31 - 1:  # pragma: NO COVER
    IS_64_BIT = False
else:  # pragma: NO COVER
    raise ImportError("Unexpected maxsize", sys.maxsize)

IS_MACOS = sys.platform == "darwin"
IS_WINDOWS = os.name == "nt"
IS_LINUX = sys.platform in ("linux", "linux2")
IS_PYPY = platform.python_implementation() == "PyPy"
