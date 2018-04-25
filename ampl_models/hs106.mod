# hs106.mod	LQR2-MN-8-22
# Original AMPL coding by Elena Bobrovnikova (summer 1996 at Bell Labs).

# Heat exchanger design

# Ref.: W. Hock and K. Schittkowski, Test Examples for Nonlinear Programming
# Codes.  Lecture Notes in Economics and Mathematical Systems, v. 187,
# Springer-Verlag, New York, 1981, p. 115.

# Number of variables: 8
# Number of constraints:  22
# Objective linear
# Nonlinear constraints

param N integer, := 8;
set I := 1..N;

param a >= 0;
param b >= 0;
param c >= 0;
param d >= 0;
param e >= 0;
param f >= 0;
param g >= 0;
param h >= 0;

var x{I};

minimize obj:
     x[1] + x[2] + x[3];

s.t. c1: 1 - a * (x[4] + x[6]) >= 0;
s.t. c2: 1 - a * (x[5] + x[7] - x[4]) >= 0;
s.t. c3: 1 - b * (x[8] - x[5]) >= 0;
s.t. c4: x[1] * x[6] - c * x[4] - d * x[1] + e >= 0;
s.t. c5: x[2] * x[7] - f * x[5] - x[2] * x[4] + f * x[4] >= 0;
s.t. c6: x[3] * x[8] - g - x[3] * x[5] + h * x[5] >= 0;
s.t. c7: 100 <= x[1] <= 10000;
s.t. c8 {i in {2,3}}: 1000 <= x[i] <= 10000;
s.t. c9 {i in 4..8}: 10 <= x[i] <= 1000;


data;
param a := 0.0025;
param b := 0.01;
param c := 833.3325;
param d := 100;
param e := 83333.33;
param f := 1250;
param g := 1250000;
param h := 2500;

var x :=
    1  5000   2  5000   3  5000   4  200   5  350   6  150   7  225   8  425;

#printf "optimal solution as starting point \n";
#var x :=
#    1  579.3167
#    2  1359.943
#    3  5110.071
#    4  182.0174
#    5  295.5985
#    6  217.9799
#    7  286.4162
#    8  395.5979
#   ;

display obj;

solve;

display x;

display obj;

display obj - 7049.330923;
