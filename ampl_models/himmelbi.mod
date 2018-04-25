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

#   Source: problem 23 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.

#   SIF input: Ph. Toint, March 1991.
#              minor correction by Ph. Shott, Jan 1995.

#   classification OLR2-MN-100-12

param N:=100;
param NT:=20;
param a{1..5, 1..NT};
param b{1..NT};
param u{1..NT};
param c{1..5};
param nw := sum {i in 1..5} c[i];

var x{1..5,1..NT} <= nw;

minimize f:
	sum {j in 1..NT} u[j]*(-1.0+prod {i in 1..5} a[i,j]^x[i,j]);

subject to cons1:
	sum {i in 1..5} x[i,1] -b[1]>= 0;
subject to cons2:
	sum {i in 1..5} x[i,6] -b[6]>= 0;
subject to cons3:
	sum {i in 1..5} x[i,10] -b[10]>= 0;
subject to cons4:
	sum {i in 1..5} x[i,14] -b[14]>= 0;
subject to cons5:
	sum {i in 1..5} x[i,15] -b[15]>= 0;
subject to cons6:
	sum {i in 1..5} x[i,16] -b[16]>= 0;
subject to cons7:
	sum {i in 1..5} x[i,20] -b[20]>= 0;
subject to cons8{i in 1..5}:
	sum {j in 1..NT} x[i,j] - c[i]<= 0;

data;
param a(tr):
	1	2	3	4	5 := 
1	1	0.84	0.96	1	0.92
2	0.95	0.83	0.95	1	0.94
3	1	0.85	0.95	1	0.92
4	1	0.84	0.96	1	0.95
5	1	0.85	0.96	1	0.95
6	0.85	0.81	0.90	1	0.98
7	0.90	0.81	0.92	1	0.98
8	0.85	0.82	0.91	1	1
9	0.8	0.8	0.92	1	1
10	1	0.86	0.95	0.96	0.90
11	1	1	0.99	0.91	0.95
12	1	0.98	0.98	0.92	0.96
13	1	1	0.99	0.91	0.91
14	1	0.88	0.98	0.92	0.98
15	1	0.87	0.97	0.98	0.99
16	1	0.88	0.98	0.93	0.99
17	1	0.85	0.95	1	1
18	0.95	0.84	0.92	1	1
19	1	0.85	0.93	1	1
20	1	0.85	0.92	1	1;

param:	b	u:=
1	30	60
2	0	50
3	0	50
4	0	75
5	0	40
6	100	60
7	0	35
8	0	30
9	0	25
10	40	150
11	0	30
12	0	45
13	0	125
14	50	200
15	70	200
16	35	130
17	0	100
18	0	100
19	0	100
20	10	150;

param c:=
1	200
2	100
3	300
4	150
5	250;

solve; display f; display x;
