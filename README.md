# Building visualization

Inspired by Raphael Appenzeller's [visualization](https://www.thingiverse.com/thing:6301564) of a chamber neighborhood inside a ~A_2 building, this project collects code to compute larger neighborhoods. Mathematically this is interesting because larger balls differ in different buildings (see for instance Tits, "Spheres of radius 2 in triangle buildings). In practice the usefulness is limited by the complexity of the visualization (I am interested to hear about understandable visualization). The code computes balls in the buildings associated to PGL_3(Q_2) and PGL_3(F_2((t))). It makes use of the fact that they admit chamber-regular lattices (see Ronan, "Triangle geometries").

It consists of two parts:

1. The GAP-code `building_embedding.g` computes balls in a building together with a retraction to $R^2$.
2. The python code `augment_embedding.py` turns this into a map to R^3 by adding a third coordinate to each vertex position, depending on its position in the fiber (and the size of the fiber). It has a few extra functions:

    * An optional linear transformation can be applied with `--transform`. This transformation turns the integer coordinates of the GAP-output into orthogonal coordinates.
    * A permutation as above can be computed with `--permutations_diff` and applied with `--permutations_patch`.
    * The output can be in the above file format, in wavefront obj with `--wavefront`, or as OpenSCAD code with `--openscad` (based on Appenzeller's code).
    * One of various z-height functions can be applied with `--zstrategy`.

Typical usage would be as follows:

1. Compute the balls
  
        gap building_embedding.g
  
    This produces four files `char0_chamber.simp`, `char0_vertex.simp`, `char2_chamber.simp`, `char2_vertex.simp`.

    [The ones with `char0` correspond to the characteristic 0 building, the ones with `char2` correspond to the characteristic 2 building. The `chamber` files contains a piece of the building consisting of all chambers (triangles) that meet a fixed chamber. They are combinatorially equivalent to Appenzeller's original model and to each other. The `vertex` files contain the complex supported on all vertices of distance at most 2 from a fixed vertex. They are distinct complexes in characteristic 0 and characteristic 2. The easiest distinguishing feature is that projection to radius 1 splits in charateristic 2 but not in characteristic 0.]
  
2. Generate a trivial permutation file
  
        python3 augment_embedding.py char0_vertex.simp --permutations_diff char0_vertex.simp char0_vertex.perm
  
3. (Repeatedly) edit the file `char0_vertex.perm` to make the outcome look as instructive as possible. Then compute the visualization
  
        python3 augment_embedding.py --transform --wavefront --permutations_patch char0_vertex.perm char0_vertex.simp char0_vertex.obj
        python3 augment_embedding.py --transform --openscad --permutations_patch char0_vertex.perm char0_vertex.simp char0_vertex.scad
  
    (Alternatively one could also directly edit the `COORDINATES` part of the `.simp` file and then optionally compute the corresponding permutation using

        python3 augment_embedding.py char0_vertex_modified.simp --permutations_diff char0_vertex_original.simp char0_vertex.perm

    .)
4. Use [OpenSCAD](https://openscad.org/) to render the `.scad`-file to `.stl` or similar.

Files produced by Raphael Appenzeller via step 3 are provided as `charX_vertex_appenzeller.perm`.

## File format

The interface between the two scripts is a primitive file format specifically suited for the task. It encodes simplicial complexes with (vertex) maps to R^k, where vertices mapped to the same point are provided as ordered lists: The format is of the form

    SIMPLICES
    0 1 2 3
    0 3 4
    0 2 4 5
    COORDINATES 3
    0.0 0.0 0.0 0
    1.0 0.0 0.0 1 4
    0.0 1.0 0.0 2
    0.0 1.0 1.0 3 5

indicating that the vertices are `0` to `5`, for instance `0 2 4 5` is a (maximal) simplex and that for instance the vertices `1` and `4` are mapped to `1.0 0.0 0.0`. That way it is possible to specify a single coordinate for each vertex but also several (ordered!) vertices that live in the same fiber. The same file format without the `SIMPLICES` section is also used to specify permutations of these fibers. For instance the file

    COORDINATES 3
    0.0 0.0 0.0 0
    1.0 0.0 0.0 1 0
    0.0 1.0 0.0 0
    0.0 1.0 1.0 0 1

would specify permutations that swap the two vertices in the fiber of `1.0 0.0 0.0` but leave the fiber of `0.0 1.0 1.0` unchanged.


