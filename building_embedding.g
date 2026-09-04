LoadPackage("kbmag");

F := FreeGroup("a","b","c");
AssignGeneratorVariables(F);
# Q_2
Gchar0 := F / [a^3, b^3, c^3, (a*b)^2 * (b*a)^-1, (b*c)^2 * (c*b)^-1, (c*a)^2 * (a*c)^-1];
# F_2((t))
Gchar2 := F / [a^3, b^3, c^3, (a*b)^2 * (b*a)^-1, (c^2*b)^2 * (b*c^2)^-1, (a*c)^2 * (c*a)^-1];
Gchar2alt := F / [ a^3, b^3, c^3, (a*b^2)^2 * (b^2*a)^-1, (c*b^2)^2 * (b^2*c)^-1, (a*c^2)^2 * (c^2*a)^-1];

ComputeBuildingBall := function(group, r, arg...)
  local rws, vertex_stabilizers, stab, Vertex, chambers, vertices, chamber_indices, type, chamber, vertex, Filter;
  rws := KBMAGRewritingSystem(group);
  KnuthBendix(rws);

  if Length(arg) = 0 then
    Filter := x -> true;
  else
    Filter := arg[1];
  fi;
  
  vertex_stabilizers := List([[a,b],[b,c],[c,a]], gens -> SubgroupOfKBMAGRewritingSystem(rws, gens));
  
  for stab in vertex_stabilizers do
    KnuthBendixOnCosetsWithSubgroupRewritingSystem(rws,stab); # WithSubgroupRewritingSystem not actually needed.
  od;
  
  Vertex := function(chamber, type)
    return [ReducedCosetRepresentative(rws,vertex_stabilizers[type],chamber),type];
  end;
  
  chambers := Filtered(Concatenation(List([0..r], i -> EnumerateReducedWords(rws,i,i))),Filter);
  
  vertices := [];
  
  for chamber in chambers do
    for type in [1..3] do
      vertex := Vertex(chamber, type);
      if not vertex in vertices then
        Add(vertices,vertex);
      fi;
    od;
  od;

  chamber_indices := List(chambers, chamber -> List([1..3], type -> Position(vertices,Vertex(chamber,type))-1)); # Shift to base 0

  return [chambers, vertices, chamber_indices];
end;

A := TransposedMat(Matrix(Integers,[[0,1,0],[1,0,0],[0,0,1]]));
B := TransposedMat(Matrix(Integers,[[-1,0,0],[-1,1,0],[0,0,1]]));
C := TransposedMat(Matrix(Integers,[[1,-1,1],[0,-1,2],[0,0,1]]));
I := Matrix(Integers,[[1,0,0],[0,1,0],[0,0,1]]); # IdentityMat(3, Integers) behaves strangely
# ABA := Matrix(Integers,[[1,-1],[0,-1]]);

coxeter_base_vertices := List([[0,0,1],[0,1,1],[1,1,1]], m -> Matrix(Integers, [m]));

CoxeterGenerator := function(letter)
  if letter = a or letter = a^-1 then
    return A;
  fi;
  if letter = b or letter = b^-1 then
    return B;
  fi;
  if letter = c or letter = c^-1 then
    return C;
  fi;
end;

CoxeterElement := function(reduced_word)
  return Product([1..Length(reduced_word)],i -> CoxeterGenerator(Subword(reduced_word,i,i)), I ); 
end;

ComputeCoordinates := function(chambers, vertices)
  local coxeter_vertices, fibers, i, j, vertex, word, type, coxeter_vertex;
  coxeter_vertices := [];
  fibers := [];
  
  for i in [1..Length(vertices)] do
    vertex := vertices[i];
    word := vertex[1];
    type := vertex[2];
    coxeter_vertex := coxeter_base_vertices[type] * CoxeterElement(word);
    j := Position(coxeter_vertices,coxeter_vertex);
    if j = fail then
      Add(coxeter_vertices,coxeter_vertex);
      Add(fibers,[i-1]); # Shift to base 0
    else
      Add(fibers[j],i-1); # Shift to base 0
    fi;
  od;
  return [coxeter_vertices, fibers];
end;

WriteSimplicialComplex := function(filename, chamber_indices, coxeter_vertices, fibers)
  local output, chamber_index, vertex_index, i, coxeter_vertex;
  output := OutputTextFile(filename,false);
  AppendTo(output, "SIMPLICES\n");
  for chamber_index in chamber_indices do
    for vertex_index in chamber_index do
      AppendTo(output, vertex_index, " ");
    od;
    AppendTo(output,"\n");
  od;
  AppendTo(output, "COORDINATES 2\n");
  for i in [1..Length(coxeter_vertices)] do
    coxeter_vertex := coxeter_vertices[i];
    AppendTo(output, coxeter_vertex[1][1], " ", coxeter_vertex[1][2], " ");
    for vertex_index in fibers[i] do
      AppendTo(output,vertex_index, " ");
    od;
    AppendTo(output,"\n");
  od;
  CloseStream(output);
end;

ComputeAndWriteBuildingBall := function(filename, group, r, args...)
  local chambersverticeschamber_indices, chambers, vertices, chamber_indices, coxeter_verticesfibers, coxeter_vertices, fibers;
  if Length(args) = 0 then
    chambersverticeschamber_indices := ComputeBuildingBall(group, r);
  else
    chambersverticeschamber_indices := ComputeBuildingBall(group, r, args[1]);
  fi;
  chambers := chambersverticeschamber_indices[1];
  vertices := chambersverticeschamber_indices[2];
  chamber_indices := chambersverticeschamber_indices[3];
  coxeter_verticesfibers := ComputeCoordinates(chambers, vertices);
  coxeter_vertices := coxeter_verticesfibers[1];
  fibers := coxeter_verticesfibers[2];
  WriteSimplicialComplex(filename, chamber_indices, coxeter_vertices, fibers);
end;

# ComputeAndWriteBuildingBall("char0_rad1.simp", Gchar0, 1);
# ComputeAndWriteBuildingBall("char2_rad1.simp", Gchar2, 1);
# ComputeAndWriteBuildingBall("char0_rad2.simp", Gchar0, 2);
# ComputeAndWriteBuildingBall("char2_rad2.simp", Gchar2, 2);
# ComputeAndWriteBuildingBall("char0_rad3.simp", Gchar0, 3);
# ComputeAndWriteBuildingBall("char2_rad4.simp", Gchar2, 4);
# ComputeAndWriteBuildingBall("char0_rad4.simp", Gchar0, 4);
# ComputeAndWriteBuildingBall("char2_rad5.simp", Gchar2, 5);
# ComputeAndWriteBuildingBall("char0_rad5.simp", Gchar0, 5);
# ComputeAndWriteBuildingBall("char2_rad6.simp", Gchar2, 6);
# ComputeAndWriteBuildingBall("char0_rad6.simp", Gchar0, 6);
# ComputeAndWriteBuildingBall("char2_rad3.simp", Gchar2, 3);

chamber_ball := Union(Set(Group(A,B)),Set(Group(B,C)),Set(Group(A,C)));
vertex_ball := Union(chamber_ball, Set(Group(A,C), x -> x * B), Set(Group(B,C), x -> x * A), Set(Group(B,C), x -> x*A*B), Set(Group(A,C), x -> x*B*A));
ComputeAndWriteBuildingBall("char0_chamber.simp", Gchar0, 3, x -> CoxeterElement(x) in chamber_ball);
ComputeAndWriteBuildingBall("char0_vertex.simp", Gchar0, 5, x -> CoxeterElement(x) in vertex_ball);
ComputeAndWriteBuildingBall("char2_chamber.simp", Gchar2, 3, x -> CoxeterElement(x) in chamber_ball);
ComputeAndWriteBuildingBall("char2_vertex.simp", Gchar2, 5, x -> CoxeterElement(x) in vertex_ball);
# ComputeAndWriteBuildingBall("char2alt_chamber.simp", Gchar2, 3, x -> CoxeterElement(x) in chamber_ball);
# ComputeAndWriteBuildingBall("char2alt_vertex.simp", Gchar2, 5, x -> CoxeterElement(x) in vertex_ball);
