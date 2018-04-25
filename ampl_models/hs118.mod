param l{1..15};
param u{1..15};

var x {j in 1..15} >= l[j], <= u[j];

minimize obj: 
  sum {k in 0..4}
  (2.3*x[3*k+1] + 0.0001*x[3*k+1]^2 + 1.7*x[3*k+2] + 0.0001*x[3*k+2]^2 +
   2.2*x[3*k+3] + 0.00015*x[3*k+3]^2)
  ;

subject to constr1 {j in 1..4}: 
   0 <= x[3*j+1] - x[3*j-2] + 7 <= 13;

subject to constr2 {j in 1..4}: 
   0 <= x[3*j+2] - x[3*j-1] + 7 <= 14;

subject to constr3 {j in 1..4}: 
   0 <= x[3*j+3] - x[3*j] + 7 <= 13;

subject to constr4: x[1] + x[2] + x[3] >= 60;
subject to constr5: x[4] + x[5] + x[6] >= 50;
subject to constr6: x[7] + x[8] + x[9] >= 70;
subject to constr7: x[10] + x[11] + x[12] >= 85;
subject to constr8: x[13] + x[14] + x[15] >= 100;

data;

let l[1] :=  8;
let l[2] := 43;
let l[3] :=  3;
let {j in 4..15} l[j] := 0;

let u[1] := 21;
let u[2] := 57;
let u[3] := 16;
let {k in 1..4} u[3*k+1] :=  90;
let {k in 1..4} u[3*k+2] := 120;
let {k in 1..4} u[3*k+3] :=  60;

let x[1] := 20;
let x[2] := 55;
let x[3] := 15;
let x[4] := 20;
let x[5] := 60;
let x[6] := 20;
let x[7] := 20;
let x[8] := 60;
let x[9] := 20;
let x[10] := 20;
let x[11] := 60;
let x[12] := 20;
let x[13] := 20;
let x[14] := 60;
let x[15] := 20;

#printf "optimal solution as starting point \n";
#let x[1] := 8;
#let x[2] := 49;
#let x[3] := 3;
#let x[4] := 1;
#let x[5] := 56;
#let x[6] := 0;
#let x[7] := 1;
#let x[8] := 63;
#let x[9] := 6;
#let x[10] := 3;
#let x[11] := 70;
#let x[12] := 12;
#let x[13] := 5;
#let x[14] := 77;
#let x[15] := 18;

display obj; 

option loqo_options $loqo_options" convex";

solve;

display x;

display obj;

display obj - 664.8204500;
