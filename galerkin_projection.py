import matplotlib.pyplot as plt
import numpy as np
import seaborn

import bezier


# Perturb the origin:
NODES = np.array([
    [-1.0   , -1.0],
    [ 0.0   , -1.0],
    [ 1.0   , -1.0],
    [-1.0   ,  0.0],
    [-0.125 ,  0.0],  # WAS: [0.0, 0.0] in straight-sided case.
    [ 1.0   ,  0.0],
    [-1.0   ,  1.0],
    [ 0.0   ,  1.0],
    [ 1.0   ,  1.0],
    [-0.5   , -1.0],
    [-0.5   , -0.5],
    [ 0.0   , -0.5],
    [-1.0   , -0.5],
    [-0.5625,  0.0],  # WAS: [-0.5, 0.0] in straight-sided case.
    [ 0.5   , -1.0],
    [ 0.5   , -0.5],
    [ 1.0   , -0.5],
    [ 0.4375,  0.0],  # WAS: [0.5, 0.0] in straight-sided case.
    [-0.5   ,  0.5],
    [ 0.0   ,  0.5],
    [-1.0   ,  0.5],
    [-0.5   ,  1.0],
    [ 0.5   ,  0.5],
    [ 1.0   ,  0.5],
    [ 0.5   ,  1.0],
])
TRIANGLES = np.array([
    [0,  9, 1, 10, 11, 4],
    [0, 10, 4, 12, 13, 3],
    [1, 14, 2, 15, 16, 5],
    [1, 15, 5, 11, 17, 4],
    [3, 13, 4, 18, 19, 7],
    [3, 18, 7, 20, 21, 6],
    [4, 17, 5, 22, 23, 8],
    [4, 22, 8, 19, 24, 7]
], dtype=np.int32)
# (ref. tri.) * 0.875 + [-0.25, -0.25]
TARGET_NODES = np.array([
    [-0.25  , -0.25  ],
    [ 0.1875, -0.25  ],
    [ 0.625 , -0.25  ],
    [-0.25  ,  0.1875],
    [ 0.1875,  0.1875],
    [-0.25  ,  0.625 ],
])
NODES_TO_CONTROL_PTS = np.array([
    [ 2, 0,  0, 0, 0,  0],
    [-1, 4, -1, 0, 0,  0],
    [ 0, 0,  2, 0, 0,  0],
    [-1, 0,  0, 4, 0, -1],
    [ 0, 0, -1, 0, 4, -1],
    [ 0, 0,  0, 0, 0,  2],
], dtype=np.float64) / 2.0


def get_all_surfaces():
    for triangle in TRIANGLES:
        local_nodes = NODES[triangle, :]
        control_pts = NODES_TO_CONTROL_PTS.dot(local_nodes)
        surface = bezier.Surface(control_pts)
        assert surface.is_valid
        yield surface


def main():
    # Target is straight-sided so nodes == control points.
    target = bezier.Surface(TARGET_NODES)
    all_ints = []
    for surface in get_all_surfaces():
        local_ints = target.intersect(surface)
        if local_ints:
            intersection, = local_ints
            all_ints.append((surface, intersection))

        ax = target.plot(256)
        surface.plot(256, ax=ax)
        if local_ints:
            intersection.plot(256, ax=ax)
        ax.axis('scaled')
        ax.set_xlim(-1.125, 1.125)
        ax.set_ylim(-1.125, 1.125)
        plt.show()

    print(all_ints)


if __name__ == '__main__':
    main()
