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

#   classification QUR2-RN-8-0

set M := {1..23};

param X{M};
param Y{M};

var A0 := 1.0;
var A2 := 1.0;
var A4 := 1.0;
var A6 := 1.0;
var A8 := 1.0;
var A10 := 1.0;
var A12 := 1.0;
var A14 := 1.0;

minimize palmer:
	sum {m in M} (Y[m] -
(A0 + A2*X[m]^2 + A4*X[m]^4 + A6*X[m]^6 + A8*X[m]^8 + A10*X[m]^10 +
A12*X[m]^12 + A14*X[m]^14) )^2;

data;
param X:= 
1                   -1.658063
2                   -1.570796
3                   -1.396263
4                   -1.221730
5                   -1.047198
6                   -0.872665
7                   -0.741119
8                   -0.698132
9                   -0.523599
10                  -0.349066
11                  -0.174533
12                  0.0
13                  0.174533
14                  0.349066
15                  0.523599
16                  0.698132
17                  0.741119
18                  0.872665
19                  1.047198
20                  1.221730
21                  1.396263
22                  1.570796
23                  1.658063;

param Y:=
1                  67.27625
2                  52.8537
3                  30.2718
4                  14.9888
5                  5.5675
6                  0.92603
7                  0.0
8                  0.085108
9                  1.867422
10                 5.014768
11                 8.263520
12                 9.8046208
13                 8.263520
14                 5.014768
15                 1.867422
16                 0.085108
17                 0.0
18                 0.92603
19                 5.5675
20                 14.9888
21                 30.2718
22                 52.8537
23                 67.27625;

solve;

display A0,A2,A4,A6,A8,A10,A12,A14;
