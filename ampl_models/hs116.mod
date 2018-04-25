# hs116.mod	LQR2-MN-13-41
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# 3-stage membrane separation

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 124.

# Number of variables: 13
# Number of constraints: 41
# Objective linear
# Nonlinear constraints


param N > 0 integer, := 13;
set I := 1 .. N;

var x {i in I} >= 0;

param a > 0;
param b > 0;
param c > 0;
param d > 0;
param e > 0;
param f > 0;


minimize obj:
     x[11] + x[12] + x[13];

s.t. c1: x[3] - x[2] >= 0;
s.t. c2: x[2] - x[1] >= 0;
s.t. c3: 1 - a * x[7] + a * x[8] >= 0;
s.t. c4: x[11] + x[12] + x[13] >= 50;
s.t. c5: x[13] - b * x[10] + c * x[3] * x[10] >= 0;
s.t. c6: x[5] - d * x[2] - e * x[2] * x[5] + f * x[2]^2 >= 0;
s.t. c7: x[6] - d * x[3] - e * x[3] * x[6] + f * x[3]^2 >= 0;
s.t. c8: x[4] - d * x[1] - e * x[1] * x[4] + f * x[1]^2 >= 0;
s.t. c9: x[12] - b * x[9] + c * x[2] * x[9] >= 0;
s.t. c10: x[11] - b * x[8] + c * x[1] * x[8] >= 0;
s.t. c11: x[5] * x[7] - x[1] * x[8] - x[4] * x[7] + x[4] * x[8] >= 0;
s.t. c12: 1 - a * (x[2] * x[9] + x[5] * x[8] - x[1] * x[8] - x[6] * x[9]) -
          x[5] - x[6] >= 0;
s.t. c13: x[2] * x[9] - x[3] * x[10] - x[6] * x[9] - 500 * x[2] +
          500 * x[6] + x[2] * x[10] >= 0;
s.t. c14: x[2] - 0.9 - a * (x[2] * x[10] - x[3] * x[10]) >= 0;
s.t. c15: x[11] + x[12] + x[13] <= 250;

s.t. b1: 0.1 <= x[1] <= 1;
s.t. b2: 0.1 <= x[2] <= 1;
s.t. b3: 0.1 <= x[3] <= 1;
s.t. b4: 0.0001 <= x[4] <= 0.1;
s.t. b5: 0.1 <= x[5] <= 0.9;
s.t. b6: 0.1 <= x[6] <= 0.9;
s.t. b7: 0.1 <= x[7] <= 1000;
s.t. b8: 0.1 <= x[8] <= 1000;
s.t. b9: 500 <= x[9] <= 1000;
s.t. b10: 0.1 <= x[10] <= 500;
s.t. b11: 1 <= x[11] <= 150;
s.t. b12: 0.0001 <= x[12] <= 150;
s.t. b13: 0.0001 <= x[13] <= 150;


data;
param a := 0.002;
param b := 1.262626;
param c := 1.231059;
param d := 0.03475;
param e := 0.975;
param f := 0.00975;
var x :=
    1 0.5  2 0.8  3 0.9  4 0.1  5 0.14  6 0.5  7 489  8 80  9 650
   10 450  11 150  12 150  13 150;

display obj;

solve;

display x;

display obj;

display obj - 97.588409;
