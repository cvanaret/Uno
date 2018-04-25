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

#   SIF input: Nick Gould, 1992.

#   classification SBR2-RN-6-0

set M := {12..24};

param X{M};
param Y{M};

var A0 := 1.0;
var A2 := 1.0;
var A4 := 1.0;
var A6 := 1.0;
var B >= 0.00001, := 1.0;
var C >= 0.00001, := 1.0;

minimize palmer:
	sum {m in M} (Y[m] -
(A0 + A2*X[m]^2 + A4*X[m]^4 + A6*X[m]^6 + B/(C+X[m]^2)) )^2;

data;
param X:= 
12                 0.000000
13                 1.570796
14                 1.396263
15                 1.221730
16                 1.047198
17                 0.872665
18                 0.785398
19                 0.732789
20                 0.698132
21                 0.610865
22                 0.523599
23                 0.349066
24                 0.174533
;
param Y:=
12                 10.678659
13                 75.414511
14                 41.513459
15                 20.104735
16                  7.432436
17                  1.298082
18                  0.171300
19                  0.000000
20                  0.068203
21                  0.774499
22                  2.070002
23                  5.574556
24                  9.026378
;

solve;

display A0,A2,A4,A6,B,C;
