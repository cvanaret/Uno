param a {1..8};
param t {1..8};
param b;

var x {1..7} >= 0, <= 1.58;

var r {i in 1..8};
var q {i in 1..8};
var s {i in 1..8};

minimize obj: 
  -( sum {j in 1..7} a[j+1]*(t[j+1]-t[j])*cos(x[j]) )^2
#  -r[8]^2
  ;

subject to constr1: q[8] = 1.0e+5;
subject to constr2: s[8] = 1.0e+3;
subject to constr3: q[1] = 0;
subject to constr4: s[1] = 0;
subject to constr5 {i in 2..8}: 
    q[i] = 0.5*(t[i]-t[i-1])^2*(a[i]*sin(x[i-1]) - b) + (t[i]-t[i-1])*s[i-1]
	 + q[i-1];
subject to constr6 {i in 2..8}: 
    s[i] = (t[i]-t[i-1])*(a[i]*sin(x[i-1]) - b) + s[i-1];

#subject to constr7: r[1] = 0;
#subject to constr8 {i in 2..8}: 
#    r[i] = a[i]*(t[i]-t[i-1])*cos(x[i-1]) + r[i-1];

data;

param:a   t :=
  1   0   0
  2  50  25
  3  50  50
  4  75 100
  5  75 150
  6  75 200
  7 100 290
  8 100 380
  ;

param b := 32;

let x[1] := 0.5;
let x[2] := 0.5;
let x[3] := 0.5;
let x[4] := 0.5;
let x[5] := 0.5;
let x[6] := 0.5;
let x[7] := 0.5;

#printf "optimal solution as starting point \n";
#let x[1] := 0.5424603;
#let x[2] := 0.5290159;
#let x[3] := 0.5084506;
#let x[4] := 0.4802693;
#let x[5] := 0.4512352;
#let x[6] := 0.4091878;
#let x[7] := 0.3527847;

solve;

display x;

display obj;

display obj + 0.831079892e+9;
