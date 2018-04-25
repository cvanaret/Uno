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

#   Source:  Problem 19 in
#   J.J. More', B.S. Garbow and K.E. Hillstrom,
#   "Testing Unconstrained Optimization Software",
#   ACM Transactions on Mathematical Software, vol. 7(1), pp. 17-41, 1981.

#   See also Buckley#32 (p.78).

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-MN-11-0

param N:=11;
param M := 65;

param x_init{1..N};
var x{i in 1..N} := x_init[i];

param t{i in 1..M} := (i-1)/10;
param y{1..M};

minimize f:
	sum {i in 1..M} ( y[i]-x[1]*exp(-t[i]*x[5])-x[2]*exp(-(t[i]-x[9])^2*x[6])-x[3]*exp(-(t[i]-x[10])^2*x[7]) -
	x[4]*exp(-(t[i]-x[11])^2*x[8]) )^2;

data;
param x_init:=
1	1.3
2	0.65
3	0.65
4	0.7
5	0.6
6	3
7	5
8	7
9	2
10	4.5
11	5.5;
param y:=
1	1.366
2	1.191
3	1.112
4	1.013
5	0.991
6	0.885
7	0.831
8	0.847
9	0.786
10	0.725
11	0.746
12	0.679
13	0.608
14	0.655
15	0.616
16	0.606
17	0.602
18	0.626
19	0.651	
20	0.724
21	0.649
22	0.649
23	0.694
24	0.644	
25	0.624
26	0.661
27	0.612
28	0.558
29	0.533	
30	0.495
31	0.500
32	0.423
33	0.395
34	0.375
35	0.372
36	0.391
37	0.396
38	0.405
39	0.428
40	0.429
41	0.523
42	0.562
43	0.607
44	0.653
45	0.672
46	0.708
47	0.633	
48	0.668
49	0.645
50	0.632
51	0.591
52	0.559
53	0.597
54	0.625
55	0.739
56	0.710
57	0.729
58	0.720
59	0.636
60	0.581
61	0.428
62	0.292
63	0.162
64	0.098
65	0.054;

solve;
display f;
display x;
