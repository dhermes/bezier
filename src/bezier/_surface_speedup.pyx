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

"""Cython "wrapped" interface for `_surface_helpers`."""


import numpy as np
from numpy cimport ndarray as ndarray_t

cimport bezier._surface


cpdef enum IntersectionClassification:
    FIRST
    SECOND
    OPPOSED
    TANGENT_FIRST
    TANGENT_SECOND
    IGNORED_CORNER


def de_casteljau_one_round(
        double[::1, :] nodes, int degree,
        double lambda1, double lambda2, double lambda3):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode='fortran'] new_nodes

    num_nodes, dimension = np.shape(nodes)
    new_nodes = np.empty((num_nodes - degree - 1, dimension), order='F')

    bezier._surface.de_casteljau_one_round(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &lambda1,
        &lambda2,
        &lambda3,
        &new_nodes[0, 0],
    )

    return new_nodes


def evaluate_barycentric(
        double[::1, :] nodes, int degree,
        double lambda1, double lambda2, double lambda3):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode='fortran'] point

    num_nodes, dimension = np.shape(nodes)
    point = np.empty((1, dimension), order='F')

    bezier._surface.evaluate_barycentric(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &lambda1,
        &lambda2,
        &lambda3,
        &point[0, 0],
    )

    return point


def evaluate_barycentric_multi(
        double[::1, :] nodes, int degree,
        double[::1, :] param_vals, int dimension):
    cdef int num_nodes, num_vals
    cdef ndarray_t[double, ndim=2, mode='fortran'] evaluated

    # NOTE: We don't check that there are ``dimension`` columns.
    num_nodes, _ = np.shape(nodes)
    # NOTE: We don't check that there are 3 columns.
    num_vals, _ = np.shape(param_vals)
    evaluated = np.empty((num_vals, dimension), order='F')

    bezier._surface.evaluate_barycentric_multi(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &num_vals,
        &param_vals[0, 0],
        &evaluated[0, 0],
    )

    return evaluated


def evaluate_cartesian_multi(
        double[::1, :] nodes, int degree,
        double[::1, :] param_vals, int dimension):
    cdef int num_nodes, num_vals
    cdef ndarray_t[double, ndim=2, mode='fortran'] evaluated

    # NOTE: We don't check that there are ``dimension`` columns.
    num_nodes, _ = np.shape(nodes)
    # NOTE: We don't check that there are 2 columns.
    num_vals, _ = np.shape(param_vals)
    evaluated = np.empty((num_vals, dimension), order='F')

    bezier._surface.evaluate_cartesian_multi(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &num_vals,
        &param_vals[0, 0],
        &evaluated[0, 0],
    )

    return evaluated


def jacobian_both(double[::1, :] nodes, int degree, int dimension):
    cdef int num_nodes
    cdef ndarray_t[double, ndim=2, mode='fortran'] new_nodes

    # NOTE: We don't check that there are ``dimension`` columns.
    num_nodes, _ = np.shape(nodes)
    new_nodes = np.empty((num_nodes - degree - 1, 2 * dimension), order='F')

    bezier._surface.jacobian_both(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &new_nodes[0, 0],
    )

    return new_nodes


def jacobian_det(double[::1, :] nodes, int degree, double[::1, :] st_vals):
    cdef int num_nodes, num_vals
    cdef ndarray_t[double, ndim=1, mode='fortran'] evaluated

    # NOTE: We don't check that there are 2 columns.
    num_nodes, _ = np.shape(nodes)
    # NOTE: We don't check that there are 2 columns.
    num_vals, _ = np.shape(st_vals)

    evaluated = np.empty((num_vals,), order='F')

    bezier._surface.jacobian_det(
        &num_nodes,
        &nodes[0, 0],
        &degree,
        &num_vals,
        &st_vals[0, 0],
        &evaluated[0],
    )

    return evaluated


def specialize_surface(
        double[::1, :] nodes, int degree,
        double[:] weights_a, double[:] weights_b, double[:] weights_c):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode='fortran'] specialized

    num_nodes, dimension = np.shape(nodes)
    specialized = np.empty((num_nodes, dimension), order='F')

    bezier._surface.specialize_surface(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &weights_a[0],
        &weights_b[0],
        &weights_c[0],
        &specialized[0, 0],
    )

    return specialized


def subdivide_nodes(double[::1, :] nodes, int degree):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode='fortran'] nodes_a
    cdef ndarray_t[double, ndim=2, mode='fortran'] nodes_b
    cdef ndarray_t[double, ndim=2, mode='fortran'] nodes_c
    cdef ndarray_t[double, ndim=2, mode='fortran'] nodes_d

    num_nodes, dimension = np.shape(nodes)
    nodes_a = np.empty((num_nodes, dimension), order='F')
    nodes_b = np.empty((num_nodes, dimension), order='F')
    nodes_c = np.empty((num_nodes, dimension), order='F')
    nodes_d = np.empty((num_nodes, dimension), order='F')

    bezier._surface.subdivide_nodes_surface(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &nodes_a[0, 0],
        &nodes_b[0, 0],
        &nodes_c[0, 0],
        &nodes_d[0, 0],
    )

    return nodes_a, nodes_b, nodes_c, nodes_d


def compute_edge_nodes(double[::1, :] nodes, int degree):
    cdef int num_nodes, dimension
    cdef ndarray_t[double, ndim=2, mode='fortran'] nodes1
    cdef ndarray_t[double, ndim=2, mode='fortran'] nodes2
    cdef ndarray_t[double, ndim=2, mode='fortran'] nodes3

    num_nodes, dimension = np.shape(nodes)
    nodes1 = np.empty((degree + 1, dimension), order='F')
    nodes2 = np.empty((degree + 1, dimension), order='F')
    nodes3 = np.empty((degree + 1, dimension), order='F')

    bezier._surface.compute_edge_nodes(
        &num_nodes,
        &dimension,
        &nodes[0, 0],
        &degree,
        &nodes1[0, 0],
        &nodes2[0, 0],
        &nodes3[0, 0],
    )

    return nodes1, nodes2, nodes3
