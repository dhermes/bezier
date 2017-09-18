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

import sys
import unittest

try:
    import bezier
except ImportError:  # pragma: NO COVER
    bezier = None


WRONG_FLAGS_TEMPLATE = """\
Arrays are not Fortran contiguous
array1 flags =
{}
array2 flags =
{}
"""
WRONG_TYPE_TEMPLATE = """\
Arrays have different types
array1({}) =
{!r}
array2({}) =
{!r}
"""
WRONG_SHAPE_TEMPLATE = """\
Arrays have different shapes
array1{} =
{!r}
array2{} =
{!r}
"""
NOT_EQ_TEMPLATE = """\
Arrays not equal
array1 =
{!r}
array2 =
{!r}
"""
if sys.maxsize == 2**63 - 1:
    IS_64_BIT = True
elif sys.maxsize == 2**31 - 1:  # pragma: NO COVER
    IS_64_BIT = False
else:  # pragma: NO COVER
    raise ImportError('Unexpected maxsize', sys.maxsize)


def get_random(seed):
    import numpy as np

    return np.random.RandomState(seed=seed)  # pylint: disable=no-member


def binary_round(value, num_bits):
    # NOTE: This assumes ``value`` is not Inf/-Inf/NaN or
    #       a subnormal number.
    hex_val = value.hex()
    # NOTE: `pre` is either '' or '-'.
    pre, hex_digits = hex_val.split('0x1.')
    hex_digits, post = hex_digits.split('p')
    assert len(hex_digits) == 13

    all_bits = '{:052b}'.format(int(hex_digits, 16))
    assert len(all_bits) == 52

    truncated_bits = all_bits[:num_bits] + '0' * (52 - num_bits)
    truncated_hex = '{:013x}'.format(int(truncated_bits, 2))

    python_hex = pre + '0x1.' + truncated_hex + 'p' + post
    return float.fromhex(python_hex)


def get_random_nodes(shape, seed, num_bits):
    import functools
    import numpy as np

    random_state = get_random(seed)
    nodes = np.asfortranarray(random_state.random_sample(shape))
    # Round the nodes to ``num_bits`` bits to avoid round-off.
    to_vectorize = functools.partial(binary_round, num_bits=num_bits)
    return np.vectorize(to_vectorize)(nodes)


def ref_triangle_uniform_nodes(pts_exponent):
    import numpy as np
    import six

    # Using the exponent means that we will divide by
    # 2**exp, which can be done without roundoff (for small
    # enough exponents).
    pts_per_side = 2**pts_exponent + 1
    total = ((pts_per_side + 1) * pts_per_side) // 2
    result = np.zeros((total, 2), order='F')

    index = 0
    for y_val in six.moves.xrange(pts_per_side):
        remaining = pts_per_side - y_val
        for x_val in six.moves.xrange(remaining):
            result[index, :] = x_val, y_val
            index += 1

    result /= (pts_per_side - 1.0)
    return result


def check_plot_call(test_case, call, expected, **kwargs):
    # Unpack the call as name, positional args, keyword args
    _, positional, keyword = call
    test_case.assertEqual(keyword, kwargs)
    test_case.assertEqual(len(positional), 2)
    test_case.assertEqual(positional[0], expected[:, 0])
    test_case.assertEqual(positional[1], expected[:, 1])


def needs_helpers_speedup(test_class):
    if bezier is None:
        has_speedup = False  # pragma: NO COVER
    else:
        has_speedup = bezier._HAS_HELPERS_SPEEDUP

    decorator = unittest.skipUnless(
        has_speedup, 'No helpers speedup available')
    return decorator(test_class)


def needs_curve_speedup(test_class):
    if bezier is None:
        has_speedup = False  # pragma: NO COVER
    else:
        has_speedup = bezier._HAS_CURVE_SPEEDUP

    decorator = unittest.skipUnless(
        has_speedup, 'No curve speedup available')
    return decorator(test_class)


def needs_surface_speedup(test_class):
    if bezier is None:
        has_speedup = False  # pragma: NO COVER
    else:
        has_speedup = bezier._HAS_SURFACE_SPEEDUP

    decorator = unittest.skipUnless(
        has_speedup, 'No surface speedup available')
    return decorator(test_class)


# pylint: disable=invalid-name
def needs_curve_intersection_speedup(test_class):
    if bezier is None:
        has_speedup = False  # pragma: NO COVER
    else:
        has_speedup = bezier._HAS_CURVE_INTERSECTION_SPEEDUP

    decorator = unittest.skipUnless(
        has_speedup, 'No curve intersection speedup available')
    return decorator(test_class)
# pylint: enable=invalid-name


class NumPyTestCase(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        import numpy as np

        super(NumPyTestCase, self).__init__(*args, **kwargs)
        self.addTypeEqualityFunc(np.ndarray, self.assertArrayEqual)

    def assertArrayEqual(self, arr1, arr2, msg=None):
        import numpy as np

        if (not arr1.flags.f_contiguous or
                not arr2.flags.f_contiguous):  # pragma: NO COVER
            standard_msg = WRONG_FLAGS_TEMPLATE.format(
                arr1.flags, arr2.flags)
            self.fail(self._formatMessage(msg, standard_msg))

        if arr1.dtype is not arr2.dtype:  # pragma: NO COVER
            standard_msg = WRONG_TYPE_TEMPLATE.format(
                arr1.dtype, arr1, arr2.dtype, arr2)
            self.fail(self._formatMessage(msg, standard_msg))

        if arr1.shape != arr2.shape:  # pragma: NO COVER
            standard_msg = WRONG_SHAPE_TEMPLATE.format(
                arr1.shape, arr1, arr2.shape, arr2)
            self.fail(self._formatMessage(msg, standard_msg))

        if not np.all(arr1 == arr2):  # pragma: NO COVER
            standard_msg = NOT_EQ_TEMPLATE.format(arr1, arr2)
            self.fail(self._formatMessage(msg, standard_msg))
