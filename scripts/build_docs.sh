#!/bin/bash
#
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
#
# Build the bezier docs.

set -e

rm -r docs/python/reference
OPTIONS="members,inherited-members,undoc-members,show-inheritance"
SPHINX_APIDOC_OPTIONS="${OPTIONS}" sphinx-apidoc \
    --separate \
    --force \
    --module-first \
    --output-dir docs/python/reference \
    src/python/bezier
# Remove unused modules.rst
rm docs/python/reference/modules.rst
# Rewrite main package RST
python scripts/rewrite_package_rst.py

# If anything has changed
if [[ -n "$(git diff -- docs/)" ]]; then
    echo "sphinx-apidoc generated changes that are not checked in to version control."
    exit 1
fi

sphinx-build -W \
    -b html \
    -d docs/build/doctrees \
    docs \
    docs/build/html
echo "Build finished. The HTML pages are in docs/build/html."
