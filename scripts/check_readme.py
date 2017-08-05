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

"""Check that the current README.rst is built from the template."""


from __future__ import print_function

import difflib
import functools
import json
import os
import re


_SCRIPTS_DIR = os.path.dirname(__file__)
_ROOT_DIR = os.path.abspath(os.path.join(_SCRIPTS_DIR, '..'))
TEMPLATE_FILE = os.path.join(_ROOT_DIR, 'README.rst.template')
TEMPLATES_FILE = os.path.join(_ROOT_DIR, 'README.templates.json')
INDEX_FILE = os.path.join(_ROOT_DIR, 'docs', 'index.rst')
README_FILE = os.path.join(_ROOT_DIR, 'README.rst')
RTD_VERSION = 'latest'
REVISION = 'master'

PLAIN_CODE_BLOCK = '.. code-block:: python'
SPHINX_CODE_BLOCK1 = """\
.. testsetup:: getting-started

   import sys

   import mock
   import numpy as np

   import bezier

   # Fake the matplotlib/seaborn imports.
   plt_mod = mock.Mock(spec=['figure', 'show'])
   plt_mod.show.return_value = None
   sys.modules['matplotlib.pyplot'] = plt_mod
   mpl_mod = mock.Mock(pyplot=plt_mod, spec=[])
   sys.modules['matplotlib'] = mpl_mod
   seaborn_mod = mock.Mock(spec=['set'])
   seaborn_mod.set.return_value = None
   sys.modules['seaborn'] = seaborn_mod

.. doctest:: getting-started"""
SPHINX_CODE_BLOCK2 = """\
.. doctest:: getting-started
   :options: +NORMALIZE_WHITESPACE"""
SPHINX_CODE_BLOCK3 = '.. doctest:: getting-started'
TEST_CLEANUP = """\
.. testcleanup:: getting-started

   sys.modules.pop('matplotlib')
   sys.modules.pop('matplotlib.pyplot')
   sys.modules.pop('seaborn')

"""
INLINE_MATH_EXPR = re.compile(r':math:`(?P<math>.*?)`')
MOD_EXPR = re.compile(r':mod:`(?P<value>.*) <(?P<module>.*)>`')
DOC_EXPR = re.compile(r':doc:`(?P<value>.*) <(?P<path>.*)>`')
TOCTREE = """\
.. toctree::
   :hidden:
   :maxdepth: 4

   Bezier Package <reference/bezier>
   curve-curve-intersection
   algorithm-helpers
   development

"""
BERNSTEIN_BASIS_SPHINX = """\
.. math::

   b_{j, n} = \\binom{n}{j} s^j (1 - s)^{n - j}"""
BEZIER_DEFN_SPHINX = """\
.. math::

   B(s) = \\sum_{j = 0}^n b_{j, n} \\cdot v_j."""
SUM_TO_UNITY_SPHINX = """\
.. math::

   b_{0, n} + b_{1, n} + \\cdots + b_{n, n} =
       \\left(s + (1 - s)\\right)^n = 1."""
PYPI_IMG = """
.. |pypi| image:: https://img.shields.io/pypi/v/bezier.svg
   :target: https://pypi.python.org/pypi/bezier
   :alt: PyPI Latest"""
VERSIONS_IMG = """
.. |versions| image:: https://img.shields.io/pypi/pyversions/bezier.svg
   :target: https://pypi.python.org/pypi/bezier
   :alt: Package Versions"""
ZENODO_IMG = """
.. |zenodo| image:: https://zenodo.org/badge/73047402.svg
   :target: https://zenodo.org/badge/latestdoi/73047402
   :alt: Zenodo DOI for ``bezier``"""
JOSS_IMG = """
.. |JOSS| image:: http://joss.theoj.org/papers/10.21105/joss.00267/status.svg
   :target: https://dx.doi.org/10.21105/joss.00267
   :alt: "Journal of Open Source Science" DOI for ``bezier``"""


def inline_math(match):
    """Convert Sphinx inline math to plain reST literal.

    Args:
        match (_sre.SRE_Match): A match (from ``re``) to be used
            in substitution.

    Returns:
        str: The ``match`` converted to literal (i.e. code font).
    """
    return '``{}``'.format(match.group('math'))


def mod_replace(match, sphinx_modules):
    """Convert Sphinx ``:mod:`` to plain reST link.

    Args:
        match (_sre.SRE_Match): A match (from ``re``) to be used
            in substitution.
        sphinx_modules (list): List to be track the modules that have been
            encountered.

    Returns:
        str: The ``match`` converted to a link.
    """
    sphinx_modules.append(match.group('module'))
    return '`{}`_'.format(match.group('value'))


def doc_replace(match, sphinx_docs):
    """Convert Sphinx ``:doc:`` to plain reST link.

    Args:
        match (_sre.SRE_Match): A match (from ``re``) to be used
            in substitution.
        sphinx_docs (list): List to be track the documents that have been
            encountered.

    Returns:
        str: The ``match`` converted to a link.
    """
    sphinx_docs.append(match.group('path'))
    return '`{}`_'.format(match.group('value'))


def get_diff(value1, value2, name1, name2):
    """Get a diff between two strings.

    Args:
        value1 (str): First string to be compared.
        value2 (str): Second string to be compared.
        name1 (str): Name of the first string.
        name2 (str): Name of the second string.

    Returns:
        str: The full diff.
    """
    lines1 = [line + '\n' for line in value1.splitlines()]
    lines2 = [line + '\n' for line in value2.splitlines()]
    diff_lines = difflib.context_diff(
        lines1, lines2, fromfile=name1, tofile=name2)
    return ''.join(diff_lines)


def readme_verify():
    """Populate the template and compare to ``README``.

    Raises:
        ValueError: If the current README doesn't agree with the expected
            value computed from the template.
        ValueError: If the ``sphinx_modules`` encountered are not as expected.
        ValueError: If the ``sphinx_docs`` encountered are not as expected.
    """
    with open(TEMPLATE_FILE, 'r') as file_obj:
        template = file_obj.read()

    with open(TEMPLATES_FILE, 'r') as file_obj:
        templates_info = json.load(file_obj)

    img_prefix = templates_info['img_prefix'].format(revision=REVISION)
    extra_links = templates_info['extra_links'].format(
        rtd_version=RTD_VERSION, revision=REVISION)
    docs_img = templates_info['docs_img'].format(rtd_version=RTD_VERSION)
    bernstein_basis = templates_info['bernstein_basis'].format(
        img_prefix=img_prefix)
    bezier_defn = templates_info['bezier_defn'].format(img_prefix=img_prefix)
    sum_to_unity = templates_info['sum_to_unity'].format(img_prefix=img_prefix)
    expected = template.format(
        code_block1=PLAIN_CODE_BLOCK,
        code_block2=PLAIN_CODE_BLOCK,
        code_block3=PLAIN_CODE_BLOCK,
        testcleanup='',
        toctree='',
        bernstein_basis=bernstein_basis,
        bezier_defn=bezier_defn,
        sum_to_unity=sum_to_unity,
        img_prefix=img_prefix,
        extra_links=extra_links,
        docs='|docs| ',
        docs_img=docs_img,
        pypi='\n\n|pypi| ',
        pypi_img=PYPI_IMG,
        versions='|versions|\n\n',
        versions_img=VERSIONS_IMG,
        rtd_version=RTD_VERSION,
        coveralls_branch='master',
        revision=REVISION,
        zenodo='|zenodo|',
        zenodo_img=ZENODO_IMG,
        joss=' |JOSS|',
        joss_img=JOSS_IMG,
    )

    # Apply regular expressions to convert Sphinx "roles" to plain reST.
    expected = INLINE_MATH_EXPR.sub(inline_math, expected)

    sphinx_modules = []
    to_replace = functools.partial(
        mod_replace, sphinx_modules=sphinx_modules)
    expected = MOD_EXPR.sub(to_replace, expected)
    if sphinx_modules != ['bezier.curve', 'bezier.surface']:
        raise ValueError('Unexpected sphinx_modules', sphinx_modules)

    sphinx_docs = []
    to_replace = functools.partial(
        doc_replace, sphinx_docs=sphinx_docs)
    expected = DOC_EXPR.sub(to_replace, expected)
    if sphinx_docs != ['reference/bezier', 'development']:
        raise ValueError('Unexpected sphinx_docs', sphinx_docs)

    # Actually get the stored contents.
    with open(README_FILE, 'r') as file_obj:
        contents = file_obj.read()

    if contents != expected:
        err_msg = '\n' + get_diff(
            contents, expected, 'README.rst.actual', 'README.rst.expected')
        raise ValueError(err_msg)
    else:
        print('README contents are as expected.')


def docs_index_verify():
    """Populate the template and compare to ``docs/index.rst``.

    Raises:
        ValueError: If the current ``index.rst`` doesn't agree with the
            expected value computed from the template.
    """
    with open(TEMPLATE_FILE, 'r') as file_obj:
        template = file_obj.read()

    img_prefix = ''
    extra_links = ''
    docs_img = ''
    expected = template.format(
        code_block1=SPHINX_CODE_BLOCK1,
        code_block2=SPHINX_CODE_BLOCK2,
        code_block3=SPHINX_CODE_BLOCK3,
        testcleanup=TEST_CLEANUP,
        toctree=TOCTREE,
        bernstein_basis=BERNSTEIN_BASIS_SPHINX,
        bezier_defn=BEZIER_DEFN_SPHINX,
        sum_to_unity=SUM_TO_UNITY_SPHINX,
        img_prefix=img_prefix,
        extra_links=extra_links,
        docs='',
        docs_img=docs_img,
        pypi='\n\n|pypi| ',
        pypi_img=PYPI_IMG,
        versions='|versions|\n\n',
        versions_img=VERSIONS_IMG,
        rtd_version=RTD_VERSION,
        coveralls_branch='master',
        revision=REVISION,
        zenodo='|zenodo|',
        zenodo_img=ZENODO_IMG,
        joss=' |JOSS|',
        joss_img=JOSS_IMG,
    )

    with open(INDEX_FILE, 'r') as file_obj:
        contents = file_obj.read()

    if contents != expected:
        err_msg = '\n' + get_diff(
            contents, expected, 'index.rst.actual', 'index.rst.expected')
        raise ValueError(err_msg)
    else:
        print('docs/index.rst contents are as expected.')


def main():
    """Verify both ``README.rst`` and ``docs/index.rst``."""
    readme_verify()
    docs_index_verify()


if __name__ == '__main__':
    main()
