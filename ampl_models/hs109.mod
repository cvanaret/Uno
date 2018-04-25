# hs109.mod	OOR2-MY-9-26
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 118.

# Number of variables:  9
# Number of constraints:  26
# Objective separable
# Objective nonconvex
# Nonlinear constraints

param a;
param b1;
param b := sin(b1);
param c := cos(b1);

param l {1..9};
param u {1..9};

var x{j in 1..9} >= l[j], <= u[j];

minimize obj:
         3 * x[1] + 1e-6 * x[1]^3 + 2 * x[2] + .522074e-6 * x[2]^3;

s.t. C1: x[4] - x[3] + .55 >= 0;
s.t. C2: x[3] - x[4] + .55 >= 0;
s.t. C3: 2250000 - x[1]^2 - x[8]^2 >= 0;
s.t. C4: 2250000 - x[2]^2 - x[9]^2 >= 0;
s.t. C5: x[5] * x[6] * sin(-x[3] - .25) + x[5] * x[7] * sin(-x[4] - .25) +
	     2 * b * x[5]^2 - a * x[1] + 400 * a = 0;
s.t. C6: x[5] * x[6] * sin(x[3] - .25) + x[6] * x[7] * sin(x[3] - x[4] - .25) +
	     2 * b * x[6]^2 - a * x[2] + 400 * a = 0;
s.t. C7: x[5] * x[7] * sin(x[4] - .25) + x[6] * x[7] * sin(x[4] - x[3] - .25) +
	     2 * b * x[7]^2 + 881.779 * a = 0;
s.t. C8: a * x[8] + x[5] * x[6] * cos(-x[3] - .25) +
	     x[5] * x[7] * cos(-x[4] - .25) - 200 * a - 2 * c * x[5]^2 +
	     .7533e-3 * a * x[5]^2 = 0;
s.t. C9: a * x[9] + x[5] * x[6] * cos(x[3] - .25) +
	     x[6] * x[7] * cos(x[3] - x[4] - .25) - 2 * c * x[6]^2 +
	     .7533e-3 * a * x[6]^2 - 200 * a = 0;
s.t. C10: x[5] * x[7] * cos(x[4] - .25) + x[6] * x[7] * cos(x[4] - x[3] - .25) -
	     2 * c * x[7]^2 + 22.938 * a + .7533e-3 * a * x[7] ^2 = 0;
#s.t. C11 {i in {1..2}}: x[i] >= 0;
#s.t. C12 {i in {3..4}}: -.55 <= x[i] <= .55;
#s.t. C13 {i in {5..7}}: 196 <= x[i] <= 252;
#s.t. C14 {i in {8..9}}: -400 <= x[i] <= 800;


data;

param a := 50.176;
param b1 := .25;

let {j in 1..2} l[j] := 0;
let {j in 3..4} l[j] := -0.55;
let {j in 5..7} l[j] := 196;
let {j in 8..9} l[j] := -400;

let {j in 1..2} u[j] := Infinity;
let {j in 3..4} u[j] := 0.55;
let {j in 5..7} u[j] := 252;
let {j in 8..9} u[j] := 800;

let {j in 1..9} x[j] := 0;

#printf "optimal solution as starting point \n";
#let x[1] := 674.8881;
#let x[2] := 1134.170;
#let x[3] := 0.1335691;
#let x[4] := -0.3711526;
#let x[5] := 252;
#let x[6] := 252;
#let x[7] := 201.465;
#let x[8] := 426.661;
#let x[9] := 368.494;

display obj;

solve;

display x; 

display obj;

display obj - 5362.06928;
