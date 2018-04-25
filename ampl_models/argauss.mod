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

#   Source: Problem 9 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#28
#   SIF input: Ph. Toint, Dec 1989.

#   classification NOR2-AN-3-15

param xinit{1..3};
param rhs{1..15};
var x{i in 1..3} := xinit[i];

minimize f: 0;
subject to cons{i in 1..15}:
	 x[1]*exp(-0.5*x[2]*(0.5*(8-i)-x[3])^2) = rhs[i];

data;
param rhs:=
1	0.0009
2	0.0044
3	0.0175
4	0.0540
5	0.1295
6	0.2420
7	0.3521
8	0.3989
9	0.3521
10	0.2420
11	0.1295
12	0.0540
13	0.0175
14	0.0044
15	0.0009;

param xinit:=
1	0.4
2	1.0
3	0.0;

solve;
display f;
display x;

