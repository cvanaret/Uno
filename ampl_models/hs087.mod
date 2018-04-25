# hs87.mod	LOR2-RN-11-24
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# The problem given here is stated in the notation of Hock & Schittkowski,
# but is corrected to conform to the problem as stated in
#	D. M. Himmelblau, Applied Nonlinear Programming,
#	McGraw-Hill, 1972, pp. 413-414
# except that jump-discontinuities in the objective function are
# omitted by stating it as the sum of two piecewise-linear terms.
# Errors in the problem statement by Hock & Schittkowski are noted in
# comments below.  (dmg, 19970820)


# Nonlinear electrical network

# Ref.: A.R.Colville. A Comparative Study on Nonlinear Programming
# Codes. IBM Scientific Center Report 320-2949, no.6, 1968.

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 106.

# Number of variables:  11 (6 before presolve and linearization of pl terms)
# Number of constraints:  24 (16 before presolve and linearization of pl terms)
# Objective convex piece-wise linear
# Nonlinear constraints

param a;
param b;
param c;
param d0;
param d := cos(d0);
param e0;
param e := sin(e0);
param lim1 > 0;
param lim2 > 0;
param lim3 >= lim2;
param rate1 >= 0;
param rate2 >= 0;
param rate3 >= 0;
param rate4 >= 0;
param rate5 >= 0;

var x1 >= 0, <= 400;
var x2 >= 0, <= 1000;
var x3 >= 340, <= 420;
var x4 >= 340, <= 420;
var x5 >= -1000, <= 1000;	# Hock & Schittkowski say <= 10000
var x6 >= 0, <= 0.5236;

minimize obj:
    <<lim1; rate1, rate2>> x1 + <<lim2, lim3; rate3, rate4, rate5>> x2;

s.t. e1: x1 = 300 - x3*x4*cos(b - x6)/a + c*x3^2*d/a;
s.t. e2: x2 = -x3*x4*cos(b + x6)/a + c*x4^2*d/a;
s.t. e3: x5 = -x3*x4*sin(b + x6)/a + c*x4^2*e/a;
s.t. e4: 200 - x3*x4*sin(b - x6)/a + c*x3^2*e/a = 0;

data;

param a := 131.078;
param b := 1.48477;	# Hock & Schittkowski say 1.48577
param c := 0.90798;
param d0 := 1.47588;
param e0 := 1.47588;
param lim1 := 300;
param lim2 := 100;
param lim3 := 200;
param rate1 := 30;
param rate2 := 31;
param rate3 := 28;
param rate4 := 29;
param rate5 := 30;

let x1 := 390;
let x2 := 1000;
let x3 := 419.5;
let x4 := 340.5;
let x5 := 198.175;
let x6 := 0.5;

# Hock & Schittkowski have an incorrectly placed decimal point in
# the x5 component of the solution they give, which should be
# 21.307... rather than 213.07...

# The formulation stated here has initial objective value
# obj = 41490 and optimal value
# obj = 8827.5977 at

#printf "optimal solution as starting point \n";
# let x1 := 107.8119
# let x2 := 196.3187
# let x3 := 373.8307
# let x4 := 420
# let x5 := 21.30716
# let x6 := 0.153292

display obj;

solve;

display x1, x2, x3, x4, x5;

display obj;

display obj - 8827.5977;
