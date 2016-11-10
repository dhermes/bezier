# See the License for the specific language governing permissions and
# limitations under the License.
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

import functools

import numpy as np


def get_random(seed):
    # pylint: disable=no-member
    return np.random.RandomState(seed=seed)
    # pylint: enable=no-member


def binary_round(value, num_bits):
    # NOTE: This assumes ``value`` is not Inf/-Inf/NaN or
    #       a subnormal number.
    hex_val = value.hex()
    pre, hex_digits = hex_val.split('0x1.')
    hex_digits, post = hex_digits.split('p')
    assert len(hex_digits) == 13

    # NOTE: Python produces 0b100101... instead of 100101...
    all_bits = bin(int(hex_digits, 16))[2:]
    all_bits = all_bits.zfill(52)
    assert len(all_bits) == 52

    truncated_bits = all_bits[:num_bits] + '0' * (52 - num_bits)
    truncated_hex = '{:013x}'.format(int(truncated_bits, 2))

    python_hex = pre + '0x1.' + truncated_hex + 'p' + post
    return float.fromhex(python_hex)


def get_random_nodes(shape, seed, num_bits):
    random_state = get_random(seed)
    nodes = random_state.random_sample(shape)
    # Round the nodes to ``num_bits`` bits to avoid round-off.
    to_vectorize = functools.partial(binary_round, num_bits=num_bits)
    return np.vectorize(to_vectorize)(nodes)
