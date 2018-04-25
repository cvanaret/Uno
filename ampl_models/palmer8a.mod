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

set M := {12..23};

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
13                 0.174533
14                 0.314159
15                 0.436332
16                 0.514504
17                 0.610865
18                 0.785398
19                 0.959931
20                 1.134464
21                 1.308997
22                 1.483530
23                 1.570796
;
param Y:=
12                  4.757534
13                  3.121416
14                  1.207606
15                  0.131916
16                  0.000000
17                  0.258514
18                  3.380161
19                 10.762813
20                 23.745996
21                 44.471864
22                 76.541947
23                 97.874528
;

solve;

display A0,A2,A4,A6,B,C;
