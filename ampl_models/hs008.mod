var x {1..2};

minimize obj: 
  -1
  ;

subject to constr1: x[1]^2 + x[2]^2 = 25;
subject to constr2: x[1]*x[2] = 9;

let x[1] := 2;
let x[2] := 1;

#printf "optimal solution as starting point \n";
#let x[1] := 4.60159;
#let x[2] := 1.95584;

solve;

display x;

display obj;

display obj + 1;
