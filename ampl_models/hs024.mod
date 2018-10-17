var x {1..2} >= 0;

minimize obj: 
  ((x[1] - 3)^2 - 9) * x[2]^3 / (27*sqrt(3))
  ;

subject to constr1: x[1]/sqrt(3) - x[2] >= 0;
subject to constr2: x[1] + sqrt(3)*x[2] >= 0;
subject to constr3: -x[1] - sqrt(3)*x[2] >= -6;

let x[1] := 1;
let x[2] := 1/2;

#printf "optimal solution as starting point \n";
#let x[1] := 3;
#let x[2] := 1.73205;

let constr1 := 0.9;
let constr3 := 0.5;

solve;

display x;
display _conname, _con;

display obj;

display obj + 1;
