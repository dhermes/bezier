# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.

"""Utilities for running functional tests as scripts."""


import inspect
import types

import numpy as np
import six


def assert_close(approximated, exact):
    """Assert that two floating point values are close.

    Makes sure the error is isolated to the last 3 bits.

    Args:
        approximated (float): The value that was computed.
        exact (float): The expected value.
    """
    local_epsilon = np.spacing(exact)  # pylint: disable=no-member
    assert abs(approximated - exact) < 8.0 * abs(local_epsilon)


def _start_line(func):
    """Get the start line (in source) of a function.

    Args:
        func (~types.FunctionType): A Python function object.

    Returns:
        int: The start line (in source).
    """
    _, line = inspect.getsourcelines(func)
    return line


def _run(mod_globals):
    """Run all tests, in source order.

    Args:
        mod_globals (dict): The globals from a module.
    """
    found = []
    for key, value in six.iteritems(mod_globals):
        if not key.lower().startswith('test'):
            continue
        if not isinstance(value, types.FunctionType):
            continue
        found.append(value)

    found.sort(key=_start_line)

    for func in found:
        func()


class Config(object):  # pylint: disable=too-few-public-methods
    """Run-time configuration.

    This is a mutable stand-in to allow test set-up to modify
    global state.
    """

    def __init__(self):
        self.running = False

    def run(self, mod_globals):
        """Run all tests, in source order.

        Args:
            mod_globals (dict): The globals from a module.
        """
        running = self.running
        try:
            self.running = True
            _run(mod_globals)
        finally:
            self.running = running
