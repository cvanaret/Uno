var x {1..3};

minimize obj: 
  0.2*x[3] - 0.8*x[1]
  ;

subject to constr1: x[2] - exp(x[1]) >= 0;
subject to constr2: x[3] - exp(x[2]) >= 0;
subject to constr3:    0 <= x[1] <= 100;
subject to constr4:    0 <= x[2] <= 100;
subject to constr5:    0 <= x[3] <=  10;

let x[1] := 0;
let x[2] := 1.05;
let x[3] := 2.9;

#printf "optimal solution as starting point \n";
#let x[1] := 0.1841264879;
#let x[2] := 1.202167873;
#let x[3] := 3.327322322;

display obj;

solve;

display x;

display obj;

display obj - 0.5181632741;
