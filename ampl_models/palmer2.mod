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

#   classification SBR2-RN-4-0

set M := {1..23};

param X{M};
param Y{M};

var A := 1.0;
var B >= 0, := 1.0;
var C >= 0, := 1.0;
var D >= 0, := 1.0;

minimize palmer:
	sum {m in M} ( Y[m] - 
(A*(X[m]^2) + ( B / ( C + (X[m]^2)/D ) ) ) )^2;

data;  
param X:=
1                   -1.745329
2                   -1.570796
3                   -1.396263
4                   -1.221730
5                   -1.047198
6                   -0.937187
7                   -0.872665
8                   -0.698132
9                   -0.523599
10                  -0.349066
11                  -0.174533
12                  0.0
13                  0.174533
14                  0.349066
15                  0.523599
16                  0.698132
17                  0.872665
18                  0.937187
19                  1.047198
20                  1.221730
21                  1.396263
22                  1.570796
23                  1.745329;

param Y:=
1                  72.676767
2                  40.149455
3                  18.8548
4                  6.4762
5                  0.8596
6                  0.00000
7                  0.2730
8                  3.2043
9                  8.1080
10                 13.4291
11                 17.714
12                 19.4529
13                 17.7149
14                 13.4291
15                 8.1080
16                 3.2053
17                 0.2730
18                 0.00000
19                 0.8596
20                 6.4762
21                 18.8548
22                 40.149455
23                 72.676767;

solve;

display A,B,C,D;
