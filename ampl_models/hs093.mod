var x {1..6} >= 0;

minimize obj: 
  0.0204*x[1]*x[4]*(x[1] + x[2] + x[3]) + 
  0.0187*x[2]*x[3]*(x[1] + 1.57*x[2] + x[4]) +
  0.0607*x[1]*x[4]*x[5]^2*(x[1] + x[2] + x[3]) +
  0.0437*x[2]*x[3]*x[6]^2*(x[1] + 1.57*x[2] + x[4])
  ;

subject to constr1: 
    0.001* prod {j in 1..6} x[j] >= 2.07;

subject to constr2: 
    0.00062*x[1]*x[4]*x[5]^2*(x[1] + x[2] + x[3]) +
    0.00058*x[2]*x[3]*x[6]^2*(x[1] + 1.57*x[2] + x[4])
    <= 1;

data;

let x[1] := 5.54;
let x[2] := 4.4;
let x[3] := 12.02;
let x[4] := 11.82;
let x[5] := 0.702;
let x[6] := 0.852;

#printf "optimal solution as starting point \n";
#let x[1] := 5.332666;
#let x[2] := 4.656744;
#let x[3] := 10.43299;
#let x[4] := 12.08230;
#let x[5] :=  0.7526074;
#let x[6] :=  0.87865084;

solve;

display x;

display obj;

display obj - 135.075961;
