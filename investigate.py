import numpy as np

import bezier
from bezier import _intersection_helpers
import bezier.surface
from bezier import _surface_helpers


GEOMETRIC = _intersection_helpers.IntersectionStrategy.GEOMETRIC
ACCEPTABLE = _surface_helpers._ACCEPTABLE
NODES = {
    '1L': np.asfortranarray([
        [0.0, 0.0],
        [1.0, 0.0],
        [0.0, 1.0],
    ]),
    '1Lrotate1': np.asfortranarray([
        [1.0, 0.0],
        [0.0, 1.0],
        [0.0, 0.0],
    ]),
    '1Lrotate2': np.asfortranarray([
        [0.0, 1.0],
        [0.0, 0.0],
        [1.0, 0.0],
    ]),
    '9L': np.asfortranarray([
        [0.0, 1.0 / 8],
        [7.0 / 8, 0.0],
        [1.0 / 4, 3.0 / 4],
    ]),
    '10Q': np.asfortranarray([
        [1.0, 0.0],
        [1.0 / 2, 1.0 / 4],
        [0.0, 0.0],
        [3.0 / 4, -3.0 / 8],
        [1.0 / 4, -3.0 / 8],
        [1.0 / 2, -3.0 / 4],
    ]),
    '17Q': np.asfortranarray([
        [1.0 / 2, -3.0 / 4],
        [83.0 / 128, -7.0 / 16],
        [51.0 / 64, -1.0 / 8],
        [45.0 / 128, -7.0 / 16],
        [1.0 / 2, 0.0],
        [13.0 / 64, -1.0 / 8],
    ]),
}


def get_intersections(surface1, surface2):
    edges1 = surface1._get_edges()
    edges2 = surface2._get_edges()
    intersections, duplicates = bezier.surface._surface_intersections(
        edges1, edges2, GEOMETRIC)
    _surface_helpers.verify_duplicates(duplicates, intersections)
    return edges1, edges2, intersections


def all_of_it(key1, key2):
    name = '{}-{}'.format(key1, key2)
    print(name)
    nodes1 = NODES[key1]
    nodes2 = NODES[key2]
    surface1 = bezier.Surface.from_nodes(nodes1, _copy=False)
    surface2 = bezier.Surface.from_nodes(nodes2, _copy=False)
    edges1, edges2, intersections = get_intersections(surface1, surface2)

    for intersection in intersections:
        print(intersection.interior_curve)
        assert intersection.interior_curve in ACCEPTABLE

    matches = []
    unused = intersections[::]
    while unused:
        start = unused.pop()
        curr_node = start
        next_node = _surface_helpers.get_next(
            start, intersections, unused)
        edge_ends = [(curr_node, next_node)]
        while next_node is not start:
            curr_node = _surface_helpers.to_front(
                next_node, intersections, unused)
            # NOTE: We also check to break when moving a corner node
            #       to the front. This is because ``intersections``
            #       de-duplicates corners by selecting the one
            #       (of 2 or 4 choices) at the front of segment(s).
            if curr_node is start:
                break
            next_node = _surface_helpers.get_next(
                curr_node, intersections, unused)
            edge_ends.append((curr_node, next_node))
            if len(edge_ends) > 10:
                raise RuntimeError(
                    'Unexpected number of edges', len(edge_ends))

        edge_info = tuple(
            _surface_helpers.ends_to_curve(start_node, end_node)
            for start_node, end_node in edge_ends
        )
        matches.append(edge_info)

    assert len(matches) == 1
    match = matches[0]
    for row in match:
        print(row)


def main():
    separator = '-' * 60
    all_of_it('10Q', '17Q')
    print(separator)
    all_of_it('17Q', '10Q')
    print(separator)
    all_of_it('1L', '9L')
    print(separator)
    all_of_it('1Lrotate1', '9L')
    print(separator)
    all_of_it('1Lrotate2', '9L')
    print(separator)
    all_of_it('9L', '1L')
    print(separator)
    all_of_it('9L', '1Lrotate1')
    print(separator)
    all_of_it('9L', '1Lrotate2')


if __name__ == '__main__':
    main()
