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

"""Convert ``lcov`` report into "cobertura" format.

Using the converted report, print a table of the report and
exit with status code 0 if all files have 100% coverage.
"""

import argparse
import sys
import tempfile

import pycobertura
import lcov_cobertura


THRESHOLD = 1.0  # 100%


def report_coverage(lcov_filename):
    with open(lcov_filename, "r") as file_obj:
        contents = file_obj.read()
    converter = lcov_cobertura.LcovCobertura(contents)
    cobertura_xml = converter.convert()

    with tempfile.NamedTemporaryFile(mode="w+") as file_obj:
        file_obj.write(cobertura_xml)
        file_obj.seek(0)
        report = file_obj.name
        cobertura = pycobertura.Cobertura(report)

    reporter = pycobertura.TextReporter(cobertura)
    print(reporter.generate())
    # The status code will be the number of files under the
    # threshold.
    return sum(
        cobertura.line_rate(source_file) < THRESHOLD
        for source_file in cobertura.files()
    )


def main():
    parser = argparse.ArgumentParser(
        description="Convert lcov output to cobertura XML and report."
    )
    parser.add_argument(
        "--lcov-filename",
        dest="lcov_filename",
        required=True,
        help="Filename of `lcov` report to be converted.",
    )
    args = parser.parse_args()
    status_code = report_coverage(args.lcov_filename)
    sys.exit(status_code)


if __name__ == "__main__":
    main()
