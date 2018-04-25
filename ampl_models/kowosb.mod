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

#   Source:  Problem 15 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-MN-4-0

param N:=4;
param M:=11;

param x_init{1..N};
var x{i in 1..N} := x_init[i];

param y{1..M};
param u{1..M};

minimize f:
	sum {i in 1..M} ( y[i]-x[1]*(u[i]^2+u[i]*x[2])/(u[i]^2+u[i]*x[3]+x[4]) )^2;

data;
param x_init := 1 0.25 2 0.39 3 0.415 4 0.39;
param:
	y	u:=
1	0.1957	4.0000
2	0.1947	2.0000
3	0.1735	1.0000
4	0.1600	0.5000
5	0.0844	0.2500
6	0.0627	0.1670
7	0.0456	0.1250
8	0.0342	0.1000
9	0.0323	0.0833
10	0.0235	0.0714
11	0.0246	0.0625;

solve;
display f;
display x;
