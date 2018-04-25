var x {1..3};

minimize obj: 
  (x[1] - x[2])^2 + (x[1] + x[2] - 10)^2/9 + (x[3] - 5)^2
  ;

subject to constr1: x[1]^2 + x[2]^2 + x[3]^2 <= 48;
subject to constr2: -4.5 <= x[1] <= 4.5;
subject to constr3: -4.5 <= x[2] <= 4.5;
subject to constr4:   -5 <= x[3] <=   5;

let x[1] := -5;
let x[2] :=  5;
let x[3] :=  0;

#printf "optimal solution as starting point \n";
#let x[1] := 3.650461821;
#let x[2] := 3.65046168;
#let x[3] := 4.6204170507;

display obj;
display constr1.body, constr2.body, constr3.body;

solve;

display x;

display obj;

display obj - 0.9535288567;
