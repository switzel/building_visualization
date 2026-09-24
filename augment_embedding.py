import argparse
import sys
import numpy as np
from itertools import combinations

def parse_simplices_file(path, check_multiple_coordinates = True):
    simplices = []
    fibers = {}
    coordinate = {}
    index = {}
    section = None
    d = 0

    with open(path, 'r') as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('SIMPLICES'):
                section = 'S'
                continue
            if line.startswith('COORDINATES'):
                _, d = line.split()
                d = int(d)
                section = 'C'
                continue
            if line.startswith('#'):
                continue

            if section == 'S':
                simplices.append([int(v) for v in line.split()])
            elif section == 'C':
                parts = line.split()
                coords = tuple(float(v) for v in parts[:d])
                ids = [int(v) for v in parts[d:]]
                fibers[coords] = ids
                for idx, a in enumerate(ids):
                    coordinate[a] = coords
                    if check_multiple_coordinates and a in index:
                        raise ValueError(f'Vertex should not have more than one coordinate {a}')
                    else:
                        index[a] = (idx, len(ids))
            else:
                raise ValueError('Unexpected content before section header')

    return simplices, fibers

def apply_permutations(fibers, permutations):
    result = {}
    for coords, fiber in fibers.items():
        try:
            permutation = permutations[coords]
            permuted_fiber = [fiber[permutation[i]] for i in range(len(permutation))]
            assert(set(fiber) == set(permuted_fiber))
            result[coords] = permuted_fiber
        except KeyError:
            result[coords] = fiber
    return result

def compute_coordinate_index(fibers):
    coordinate = {}
    index = {}
    for coords, ids in fibers.items():
        for idx, a in enumerate(ids):
            coordinate[a] = coords
            if a in index:
                raise ValueError(f'Vertex should not have more than one coordinate {a}')
            else:
                index[a] = (idx, len(ids))
    return coordinate, index

def write_simplices_file(path, simplices, fibers):
    with open(path, 'w') as f:
        if simplices:
            f.write('SIMPLICES\n')
            for simplex in simplices:
                f.write(' '.join(str(int(v)) for v in simplex) + '\n')
        if fibers:
            d = len(next(iter(fibers.keys())))
            f.write(f'COORDINATES {d}\n')
            for coords, vertices in fibers.items():
                assert(len(coords) == d)
                f.write(' '.join(f'{c:.17g}' for c in coords) + '  ')
                f.write(' '.join(f'{d}' for d in vertices) + '\n')

def write_wavefront_file(path, simplices, coordinates):
    n = max(sum(simplices,[])) + 1
    vertex_coordinates = [coordinates[i] for i in range(n)]
    with open(path, 'w') as f:
        for coords in vertex_coordinates:
            f.write('v ' + ' '.join(f'{c:.17g}' for c in coords) + '\n')
        for simplex in simplices:
            f.write('f ' + ' '.join(str(int(v) + 1) for v in simplex) + '\n')

openscad_header='''
scale = 40;
vertex_radius = 1.25/scale;
edge_radius = 1.25/scale;
triangle_thickness = 1/scale;
fn = 20;

function transpose(A) = [for (j = [0:len(A[0])-1]) [for(i = [0:len(A)-1]) A[i][j]]];
module vertex(p)
{
   translate (v=p) sphere(vertex_radius, $fn = fn);
}
module triangle(p,q,r) {
    list = [p,q,r];
listu = [for(i=[0:len(list)-1]) list[i]+[0,0,triangle_thickness/2]];
listd = [for(i=[0:len(list)-1]) list[i]+[0,0,-triangle_thickness/2]];
polyhedron(concat(listu, listd), [
        [0,2,1],[3,4,5],
        [0,1,4,3],[1,2,5,4],[2,0,3,5]
        ],convexity = 1);
}
module edge(p, q){
    v = q - p;    // vector from P to Q
    L = norm(v);  // height of the cylnder = dist(P, Q) 
    c = v / L;    // unit vector: direction from P to Q    
    is_c_vertical = ( 1 - abs(c * [0, 0, 1]) < 1e-6); //is c parallel to z axis?
    u = is_c_vertical ? [1, 0, 0] : cross([0, 0, 1], c); // normal to c and Z axis
    a = u / norm(u); // unit vector normal to c and the Z axis
    b = cross(c, a); // unit vector normal to a and b
    // [a, b, c] is an orthonormal basis, i.e. the rotation matrix; P is the translation
    MT = [a, b, c, p]; // the transformation matrix
    M = transpose(MT); // OpenSCAD wants vectors in columns, so we need to transpose
    multmatrix(M)
    cylinder(h=L, r=edge_radius, $fn = fn);
}
'''

def edges_of_simplices(simplices):
    return list(set().union(*(set(combinations(simplex,2)) for simplex in simplices)))

def write_openscad_file(path, simplices, coordinates):
    n = max(sum(simplices,[])) + 1
    with open(path,'w') as f:
        f.write(openscad_header)
        f.write('vertices = [\n' + ',\n'.join(f'{list(coordinates[i])}' for i in range(n)) + '];\n')
        f.write('module skeleton() {\nunion() {\n')
        f.write('for (i=[0:len(vertices)-1]) { vertex(vertices[i]); }\n')
        edges = edges_of_simplices(simplices)
        for edge in edges:
            f.write(f'edge(vertices[{edge[0]}], vertices[{edge[1]}]);\n')
        f.write('}}\nmodule chambers() {\nunion() {\n')
        for simplex in simplices:
            f.write(f'triangle(vertices[{simplex[0]}], vertices[{simplex[1]}], vertices[{simplex[2]}]);\n')
        f.write('}}\nscale(scale) rotate([-90,0,0]) union(){\nskeleton();\nchambers();\n}\n')


def write_permutations_file(path, fibers):
    with open(path, 'w') as f:
        d = len(next(iter(fibers.keys()), ()))
        f.write(f'COORDINATES {d}\n')
        for coords, fiber in fibers.items():
            assert(len(coords) == d)
            f.write(' '.join(f'{c:.17g}' for c in coords) + ' ' + ' '.join(f' {a}' for a in fiber) + '\n')

def permutation_diff(fibers1, fibers2):
    assert(fibers1.keys() == fibers2.keys())
    result = {}
    for vertex in fibers1.keys():
        fiber1 = fibers1[vertex]
        fiber2 = fibers2[vertex]
        perm = [fiber1.index(a) for a in fiber2]
        result[vertex] = perm
    return result

def z_coordinate_fixed_width_inside_out(index, number):
    if number <= 1:
        return 0.0
    if number % 2 == 0:
        return (-1)**(index + 1) * (((index+2)//2) * (2.0/(number-1)) - (1.0/(number-1)))
    else:
        return (-1)**index * ((index+1)//2) * (2.0/(number-1))

def z_coordinate_log_width_inside_out(index, number):
    scale = len(bin(number)[2:])/4
    if number <= 1:
        return 0.0
    if number % 2 == 0:
        return (-1)**(index + 1) * (((index+2)//2) * (2.0/(number-1)) - (1.0/(number-1))) * scale
    else:
        return (-1)**index * ((index+1)//2) * (2.0/(number-1)) * scale

def z_coordinate_fixed_width_bottom_up(index, number):
    if number <= 1:
        return 0.0
    return (-1 + 2.0*index/(number - 1))

def z_coordinate_log_width_bottom_up(index, number):
    scale = len(bin(number)[2:])/4
    if number <= 1:
        return 0.0
    return (-1 + 2.0*index/(number - 1)) * scale

z_coordinates = {
        'fixed_inout' : z_coordinate_fixed_width_inside_out,
        'log_inout' : z_coordinate_log_width_inside_out,
        'fixed_up' : z_coordinate_fixed_width_bottom_up,
        'log_up' : z_coordinate_log_width_bottom_up,
        }

def transform(vector, trafo = np.array([[1,-.5,0],[0,.5*np.sqrt(3),0],[0,0,1]])):
    result = trafo @ np.array(vector + (1,))
    return tuple(result.tolist()[:-1])

def parse_args():
    parser = argparse.ArgumentParser(description="Process an input file and write to an output file.")
    parser.add_argument("input_file", help="Path to the input file")
    parser.add_argument("output_file", help="Path to the output file")
    parser.add_argument("--transform", action="store_true", help="Transform standard basis to 120 degree angle.")
    parser.add_argument("--wavefront", action="store_true", help="Output as Wavefront obj.")
    parser.add_argument("--openscad", action="store_true", help="Output as openscad.")
    parser.add_argument("--permutations_diff", type=str, help="Create a permutations file for the difference between this file and input_file.")
    parser.add_argument("--permutations_patch", type=str, help="Apply the permutations file to input_file.")
    parser.add_argument("--zstrategy", type=str, default='log_up', help="Which strategy to use for computing z-coordinates.")
    return parser.parse_args()

def main():
    args = parse_args()
    if sum([1 for exclusive in [args.wavefront, args.openscad, args.permutations_diff] if exclusive != None]) > 1:
        ValueError('Only one type of output can be specified')
    simplices, fibers = parse_simplices_file(args.input_file)
    if args.permutations_patch:
        _, permutations = parse_simplices_file(args.permutations_patch, check_multiple_coordinates = False)
        fibers = apply_permutations(fibers, permutations)
    coordinate, index = compute_coordinate_index(fibers)
    if args.transform:
        apply_transform = transform
    else:
        apply_transform = lambda x : x
    z_coordinate = z_coordinates[args.zstrategy]
    new_coordinate = { vx : apply_transform(coord) + (z_coordinate(*index[vx]),) for vx, coord in coordinate.items() }
    if args.wavefront:
        write_wavefront_file(args.output_file, simplices, new_coordinate)
    elif args.permutations_diff:
        other_simplices, other_fibers = parse_simplices_file(args.permutations_diff)
        diff = permutation_diff(other_fibers, fibers)
        write_permutations_file(args.output_file, diff)
    elif args.openscad:
        write_openscad_file(args.output_file, simplices, new_coordinate)
    else:
        write_simplices_file(args.output_file, simplices, fibers)

if __name__ == "__main__":
    sys.exit(main())
