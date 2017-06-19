rearr67 := rec(domain := ["000", "001", "010", "011", "100", "101", "11"],
               range := ["00", "010", "011", "100", "101", "110", "111"]);

prod := function(R1, R2)
  local dom1, range1, dom2, range2, i, n, j, suffixes, dom_pref, ran_pref, suff;
  # Return R1 * R2
  dom1 := R1.domain;
  range1 := R1.range;
  dom2 := R2.domain;
  range2 := R2.range;
  # Assume domain and range are sorted
  i := 1;
  repeat
    if Length(range1[i]) < Length(dom2[i]) then
      # the R1-range cell contains several R2-domain cells
      n := Length(range1[i]);
      j := i;
      suffixes := [];
      while j <= Length(dom2) and StartsWith(dom2[j], range1[i]) do
        Add(suffixes, dom2[j]{[n + 1 .. Length(dom2[j])]});
        j := j + 1;
      od;
      dom_pref := Remove(dom1, i);
      ran_pref := Remove(range1, i);
      j := 0;
      for suff in suffixes do
        Add(dom1, Concatenation(dom_pref, suff), i + j);
        Add(range1, Concatenation(ran_pref, suff), i + j);
        j := j + 1;
      od;
    elif Length(range1[i]) > Length(dom2[i]) then
      # the R2-domain cell contains several R1-range cells
      n := Length(dom2[i]);
      j := i;
      suffixes := [];
      while j <= Length(range1) and StartsWith(range1[j], dom2[i]) do
        Add(suffixes, range1[j]{[n + 1 .. Length(range1[j])]});
        j := j + 1;
      od;
      dom_pref := Remove(dom2, i);
      ran_pref := Remove(range2, i);
      j := 0;
      for suff in suffixes do
        Add(dom2, Concatenation(dom_pref, suff), i + j);
        Add(range2, Concatenation(ran_pref, suff), i + j);
        j := j + 1;
      od;
    fi;
    i := i + 1;
  until i > Length(range1);

  return rec(domain := dom1, range := range2);
end;

boundary_of_cells := function(prefix1, prefix2)
  local i;
  # We assume prefix1 and prefix2 represent adjacent cells
  i := 0;
  while prefix1[i + 1] = prefix2[i + 1] do
    i := i + 1;
  od;
  return prefix1{[1 .. i]};
end;

boundary_points_of_list := function(list)
  local m, points, i, max;
  # list should be the domain or range of a rearrangement
  m := Length(list) - 1;
  points := EmptyPlist(m);
  for i in [1 .. m] do
    points[i] := boundary_of_cells(list[i], list[i + 1]);
  od;

  # max := Maximum(List(points, Length));
  # points := List([0 .. max], len -> Filtered(points, w -> Length(w) = len));
  # # note that within one word-length we expect the points to be ordered already
  # points := Concatenation(points);
  return points;
end;

# Output a list of words w1, w2, w3 ... s.t. f_w1 f_w2 f_w3... = f
factorise_rearrangement := function(f)
  local m, dom_bp, ran_bp, max, dom_bp_by_depth, ran_bp_sorted, i;
  # f is a rearrangement, which is represented as two lists where
  #   f.domain contains prefixes for the domain cells and
  #   f.range contains prefixes for the range cells.
  # f.domain[i] maps to f.range[i] for each i

  # Let m+1 be the length of f.domain (= length of f.range),
  # so m is the number of boundary points.
  m := Length(f.domain) - 1;
  dom_bp := boundary_points_of_list(f.domain);
  ran_bp := boundary_points_of_list(f.range);

  # Sort dom_bp by depth order
  max := Maximum(List(dom_bp, Length));
  dom_bp_by_depth := List([0 .. max],
                          len -> Filtered(dom_bp, w -> Length(w) = len));
  # note that within one word-length we expect the points to be ordered already
  dom_bp_by_depth := Concatenation(dom_bp_by_depth);

  # Sort ran_bp by the same order
  ran_bp_sorted := EmptyPlist(m);
  for i in [1 .. m] do
    ran_bp_sorted[Position(dom_bp_by_depth, dom_bp[i])] := ran_bp[i];
  od;

  # Go through all the boundary points
  for i in [1 .. m] do
#    ran_bp_sorted[i]
  od;

  return [dom_bp_by_depth, ran_bp_sorted];
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

f_alpha := function(alpha)
  local domain, range, opposite, i, cell, R;
  # alpha should be a string of 0s and 1s
  domain := ["00", "01", "1"];
  range := ["0", "10", "11"];
  domain := List(domain, str -> Concatenation(alpha, str));
  range := List(range, str -> Concatenation(alpha, str));

  opposite := function(char)
    if char = '0' then
      return '1';
    fi;
    return '0';
  end;

  for i in [1 .. Length(alpha)] do
    cell := Concatenation(alpha{[1 .. i-1]}, [opposite(alpha[i])]);
    Add(domain, cell);
    Add(range, cell);
  od;

  Sort(domain);
  Sort(range);
  R := rec(domain := domain, range := range);
  return R;
end;
