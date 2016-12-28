# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Setup file for bezier."""

import os

from setuptools import find_packages
from setuptools import setup


VERSION = '0.2.1'
PACKAGE_ROOT = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(PACKAGE_ROOT, 'README.rst.template')) as file_obj:
    TEMPLATE = file_obj.read()

README = TEMPLATE.format(
    pypi='',
    pypi_img='',
    versions='',
    versions_img='',
    rtd_version=VERSION,
    coveralls_branch=VERSION,
)

REQUIREMENTS = (
    'enum34',
    'matplotlib >= 1.5.3',
    'numpy >= 1.11.2',
    'six >= 1.9.0',
)
DESCRIPTION = (
    u'Helper for B\u00e9zier Curves, Triangles, and Higher Order Objects')


setup(
    name='bezier',
    version=VERSION,
    description=DESCRIPTION,
    author='Danny Hermes',
    author_email='daniel.j.hermes@gmail.com',
    long_description=README,
    scripts=(),
    url='https://github.com/dhermes/bezier',
    packages=find_packages(),
    license='Apache 2.0',
    platforms='Posix; MacOS X; Windows',
    include_package_data=True,
    zip_safe=True,
    install_requires=REQUIREMENTS,
    classifiers=(
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Topic :: Internet',
    ),
)
