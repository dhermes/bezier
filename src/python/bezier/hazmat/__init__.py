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

"""The ``bezier.hazmat`` subpackage is a "hazardous materials" layer.

The modules contained here provide low-level functionality for the public
API of the ``bezier`` package. They are not subject to any API stability
guarantees. Understanding these low-level details is not required to be an
everyday user of ``bezier``, but documenting them will allow examples and
narrative for those interested on learning more.

In many cases, the code is a pure-Python implementation of functions provided
in the :doc:`binary extension <../binary-extension>` speedup.
"""
