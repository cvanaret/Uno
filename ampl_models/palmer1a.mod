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

#   Source:
#   M. Palmer, Edinburgh, private communication.

#   SIF input: Nick Gould, 1990.

#   classification SBR2-RN-6-0

set M := {1..35};

param X{M};
param Y{M};

var A0 := 1.0;
var A2 := 1.0;
var A4 := 1.0;
var A6 := 1.0;
var B >= 0, := 1.0;
var C >= 0, := 1.0;

minimize palmer:
	sum {m in M} ( Y[m] - A0 - A2*(X[m]^2) - A4*(X[m]^4) - A6*(X[m]^6) - ( B / ( C + (X[m]^2) ) ) )^2;

data;  
param X:=
1 -1.788963
2 -1.745329
3 -1.658063
4 -1.570796
5 -1.483530
6 -1.396263
7 -1.308997
8 -1.218612
9 -1.134464
10 -1.047198
11 -0.872665
12 -0.698132
13 -0.523599
14 -0.349066
15 -0.174533
16 0.0000000
17 1.788963
18 1.745329
19 1.658063
20 1.570796
21 1.483530
22 1.396263
23 1.308997
24 1.218612
25 1.134464
26 1.047198
27 0.872665
28 0.698132
29 0.523599
30 0.349066
31 0.174533
32 -1.8762289
33 -1.8325957
34 1.8762289
35 1.8325957;

param Y:=
1 78.596218
2 65.77963
3 43.96947
4 27.038816
5 14.6126
6 6.2614
7 1.538330
8 0.000000
9 1.188045
10 4.6841
11 16.9321
12 33.6988
13 52.3664
14 70.1630
15 83.4221
16 88.3995
17 78.596218
18 65.77963
19 43.96947
20 27.038816
21 14.6126
22 6.2614
23 1.538330
24 0.000000
25 1.188045
26 4.6841
27 16.9321
28 33.6988
29 52.3664
30 70.1630
31 83.4221
32 108.18086
33 92.733676
34 108.18086
35 92.733676;


solve;

display A0,A2,A4,A6,B,C;
