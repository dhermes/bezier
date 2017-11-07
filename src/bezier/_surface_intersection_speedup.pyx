#!python
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

"""Cython "wrapped" interface for `_surface_intersection`."""


import numpy as np

cimport bezier._surface_intersection


def newton_refine(
        double[::1, :] nodes, int degree,
        double x_val, double y_val, double s, double t):
    cdef int num_nodes
    cdef double updated_s, updated_t

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)

    bezier._surface_intersection.newton_refine_surface(
        &num_nodes,
        &nodes[0, 0],
        &degree,
        &x_val,
        &y_val,
        &s,
        &t,
        &updated_s,
        &updated_t,
    )

    return updated_s, updated_t


def locate_point(double[::1, :] nodes, int degree, double x_val, double y_val):
    cdef int num_nodes
    cdef double s_val, t_val

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)

    bezier._surface_intersection.locate_point_surface(
        &num_nodes,
        &nodes[0, 0],
        &degree,
        &x_val,
        &y_val,
        &s_val,
        &t_val,
    )

    if s_val == -1.0:  # LOCATE_MISS
        return None
    else:
        return s_val, t_val
