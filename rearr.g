rearr67 := rec(domain := ["000", "001", "010", "011", "100", "101", "11"],
               range := ["00", "010", "011", "100", "101", "110", "111"]);

boundary_of_cells := function(prefix1, prefix2)
  local i;
  # We assume prefix1 and prefix2 represent adjacent cells
  i := 0;
  while prefix1[i + 1] = prefix2[i + 1] do
    i := i + 1;
  od;
  return prefix1{[1 .. i]};
end;

sorted_boundary_points_of_list := function(list)
  local m, points, i, max;
  # list should be the domain or range of a rearrangement
  m := Length(list) - 1;
  points := EmptyPlist(m);
  for i in [1 .. m] do
    points[i] := boundary_of_cells(list[i], list[i + 1]);
  od;
  
  max := Maximum(List(points, Length));
  points := List([0 .. max], len -> Filtered(points, w -> Length(w) = len));
  # note that within one word-length we expect the points to be ordered already
  points := Concatenation(points);
  return points;
end;

# Output a list of words w1, w2, w3 ... s.t. f_w1 f_w2 f_w3... = f
factorise_rearrangement := function(f)
  local m, dom_bp, ran_bp;
  # f is a rearrangement, which is represented as two lists where
  #   f.domain contains prefixes for the domain cells and
  #   f.range contains prefixes for the range cells.
  # f.domain[i] maps to f.range[i] for each i
  
  # Let m+1 be the length of f.domain (= length of f.range),
  # so m is the number of boundary points.
  m := Length(f.domain) - 1;
  dom_bp := sorted_boundary_points_of_list(f.domain);
  ran_bp := sorted_boundary_points_of_list(f.range);

  return [dom_bp, ran_bp];
end;

random_cellular_partition := function(nr_cells)
  local nr_splits, list, split, i;
  nr_splits := nr_cells - 1;
  list := [""];
  for split in [1 .. nr_splits] do
    i := Random([1 .. Length(list)]);
    Add(list, Concatenation(list[i], "0"));
    list[i] := Concatenation(list[i], "1");
  od;
  return SortedList(list);
end;

random_rearrangement := function(nr_cells)
  return rec(domain := random_cellular_partition(nr_cells),
             range := random_cellular_partition(nr_cells));
end;
