# hs114.mod	QOR2-MY-10-31
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Alkylation process

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 123.

# Number of variables:  10
# Number of constraints:  31
# Objective separable
# Objective nonconvex
# Nonlinear constraints

param n := 10;
set I := 1..n;

param lb{I};
param ub{i in I} > lb[i];
param x0{I};

param a := .99;
param b := .9;

var x{i in I} >= lb[i] <= ub[i] := x0[i];

minimize obj: 5.04*x[1] + .035*x[2] + 10*x[3] + 3.36*x[5] - .063*x[4]*x[7];

var G1 = 35.82 - .222*x[10] - b*x[9];
var G2 = -133 + 3*x[7] - a*x[10];
var G5 = 1.12*x[1] + .13167*x[1]*x[8] - .00667*x[1]*x[8]^2 - a*x[4];
var G6 = 57.425 + 1.098*x[8] - .038*x[8]^2 + .325*x[6] - a*x[7];
s.t. g1:  G1 >= 0;
s.t. g2:  G2 >= 0;
s.t. g3:  -G1 + x[9]*(1/b - b) >= 0;
s.t. g4:  -G2 + (1/a - a)*x[10] >= 0;
s.t. g5:  G5 >= 0;
s.t. g6:  G6 >= 0;
s.t. g7:  -G5 + (1/a - a)*x[4] >= 0;
s.t. g8:  -G6 + (1/a - a)*x[7] >= 0;
s.t. g9:  1.22*x[4] - x[1] - x[5] = 0;
s.t. g10: 98000*x[3]/(x[4]*x[9] + 1000*x[3]) - x[6] = 0;
s.t. g11: (x[2] + x[5])/x[1] - x[8] = 0;

data;
param :	lb	ub	x0 :=
1	.00001	2000	1745
2	.00001	16000	12000
3	.00001	120	110
4	.00001	5000	3048
5	.00001	2000	1974
6	85	93	89.2
7	90	95	92.8
8	3	12	8
9	1.2	4	3.6
10	145	162	145
;

display obj;

solve;

display x;

display obj;

display obj + 1768.80696;
