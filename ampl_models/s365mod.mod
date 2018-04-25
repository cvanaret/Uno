# AMPL Model by Hande Y. Benson
#
# Copyright (C) 2001 Princeton University
# All Rights Reserved
#
# Permission to use, copy, modify, and distribute this software and
# its documentation for any purpose and without fee is hereby
# granted, provided that the above copyright notice appear in all
# copies and that the copyright notice and this
# permission notice appear in all supporting documentation.                     

#   Source: a modification (Ph. Toint) of problem 365 in
#   K. Schittkowski,
#   " More Test Problems for Nonlinear Programming Codes",
#   Springer Verlag, Berlin, 1987.

#   SIF input: Ph. Toint, March 1991.

#   classification QOR2-AY-7-5

param N := 7;

var x{1..N};
var P = sqrt(x[2]^2) + x[3]^2;
var Q = sqrt(x[3]^2) + (x[2]-x[1])^2;

minimize f:
x[1]*x[3];

subject to cons1:
(x[4]-x[6])^2 + (x[5]-x[7])^2 - 4 >= 0;

subject to cons2:
(x[3]*x[4] - x[2]*x[5])/P - 1 >= 0;

subject to cons3:
(x[3]*x[6] - x[2]*x[7])/P - 1 >= 0;

subject to cons4:
(x[1]*x[3] + (x[2]-x[1])*x[5] - x[3]*x[4])/Q - 1 >= 0;

subject to cons5:
(x[1]*x[3] + (x[2]-x[1])*x[7] - x[3]*x[6])/Q - 1 >= 0;

subject to cons6:
x[1] >= 0.5;

subject to cons7:
x[3] >= 0.5;

subject to cons8:
x[5] >= 1;

subject to cons9:
x[7] >= 1;

data;
var x:=
1	3
2	0.01
3	2
4	-1.5
5	1.5
6	5
7	0;

solve;
display f;
display x;
