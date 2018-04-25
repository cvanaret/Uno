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

#   Source:  Problem 17 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See alos Buckley#32 (p. 77).

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-MN-5-0

param N:=5;
param M:=33;

param t{i in 1..M} := 10*(i-1);
param y{1..M};

param x_init{1..N};
var x{i in 1..N} := x_init[i];

minimize f:
	sum {i in 1..M} ( y[i]-x[1]-x[2]*exp(-t[i]*x[4])-x[3]*exp(-t[i]*x[5]) )^2;

data;
param x_init := 1 0.5 2 1.5 3 -1 4 0.01 5 0.02;
param y:=
1	0.844
2	0.908
3	0.932
4	0.936
5	0.925
6	0.908
7	0.881
8	0.850
9	0.818
10	0.784
11	0.751
12	0.718
13	0.685
14	0.658
15	0.628
16	0.603
17	0.580
18	0.558
19	0.538
20	0.522
21	0.506
22	0.490
23	0.478
24	0.467
25	0.457
26	0.448
27	0.438
28	0.431
29	0.424
30	0.420
31	0.414
32	0.411
33	0.406;

option loqo_options "verbose=2 timing=1 convex";
solve;
display f;
display x;

