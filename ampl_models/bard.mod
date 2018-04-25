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

#   Source: Problem 3 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#16.
#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-3-0

param N := 3;
param M := 15;

param y{1..M};
param u{i in 1..M} := i;
param v{i in 1..M} := 16-i;
param w{i in 1..M} := min(u[i],v[i]);

var x{1..N} := 1;

minimize f:
	sum {i in 1..M} ( y[i]-(x[1]+u[i]/(v[i]*x[2]+w[i]*x[3])) )^2;

data;
param y:=
1	0.14
2	0.18
3	0.22
4	0.25
5	0.29
6	0.32
7	0.35
8	0.39
9	0.37
10	0.58
11	0.73
12	0.96
13	1.34
14	2.10
15	4.39;

solve;
display f;
display x;
