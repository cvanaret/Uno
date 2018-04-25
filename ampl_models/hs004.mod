var x {1..2};

minimize obj: 
  (x[1]+1)^3/3 + x[2]
  ;

subject to constr1: 1 <= x[1];
subject to constr2: 0 <= x[2];

let x[1] := 1.125;
let x[2] := 0.125;

solve;

display x;

display obj;

display obj - 8/3;
