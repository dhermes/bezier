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
"""Add some features to ``sphinx.ext.doctest``.

Does so by monkey-patching the ``option_spec`` in

* ``DoctestDirective``
* ``TestcodeDirective``
* ``TestoutputDirective``

with some extra options:

* ``:linux-only:``
* ``:mac-os-x-only:``
* ``:windows-skip:``
* ``:windows-only:``

Also monkey-patches ``TestDirective.run`` to honor these directives
and skip a test based on the options.

.. note::

   This works with version(s) 1.6.4 of Sphinx, but may not work with
   other versions (i.e. the monkey-patch relies on some knowledge of
   the implementation).
"""

import doctest
import os
import sys

import docutils.parsers.rst
import sphinx.ext.doctest

IS_WINDOWS = os.name == 'nt'
IS_LINUX = sys.platform in ('linux', 'linux2')
IS_MAC_OS_X = sys.platform == 'darwin'
LINUX_ONLY = 'linux-only'
MAC_OS_X_ONLY = 'mac-os-x-only'
WINDOWS_ONLY = 'windows-only'
WINDOWS_SKIP = 'windows-skip'
OPTION_NAMES = (LINUX_ONLY, MAC_OS_X_ONLY, WINDOWS_ONLY, WINDOWS_SKIP)
OLD_RUN = sphinx.ext.doctest.TestDirective.run
SKIP_FLAG = doctest.OPTIONFLAGS_BY_NAME['SKIP']


def custom_run(directive):
    """Custom over-ride for :meth:`.TestDirective.run`.

    ``directive`` acts like ``self`` when this function is bound to a class.

    Helps to skip tests based on the directive options:

    * ``windows-only``
    * ``windows-skip``

    Args:
        directive (sphinx.ext.doctest.TestDirective): The currently active
            Sphinx directive.

    Returns:
        docutils.nodes.Element: The element to be added.
    """
    node, = OLD_RUN(directive)
    num_options = sum(1 for name in OPTION_NAMES if name in directive.options)
    if num_options > 1:
        raise RuntimeError(
            'At most one option can be used among', *OPTION_NAMES
        )

    if LINUX_ONLY in directive.options:
        if not IS_LINUX:
            node['options'][SKIP_FLAG] = True
    if MAC_OS_X_ONLY in directive.options:
        if not IS_MAC_OS_X:
            node['options'][SKIP_FLAG] = True
    if WINDOWS_ONLY in directive.options:
        if not IS_WINDOWS:
            node['options'][SKIP_FLAG] = True
    if WINDOWS_SKIP in directive.options:
        if IS_WINDOWS:
            node['options'][SKIP_FLAG] = True
    return [node]


def setup(app):
    """Set-up this extension.

    Args:
        app (sphinx.application.Sphinx): A running Sphinx app.
    """
    sphinx.ext.doctest.TestDirective.run = custom_run
    directive_types = (
        sphinx.ext.doctest.DoctestDirective,
        sphinx.ext.doctest.TestcodeDirective,
        sphinx.ext.doctest.TestoutputDirective,
    )
    for directive in directive_types:
        option_spec = directive.option_spec
        for option in OPTION_NAMES:
            if option in option_spec:
                raise RuntimeError('Unexpected option in option spec', option)

            option_spec[option] = docutils.parsers.rst.directives.flag
