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

#   classification SBR2-RN-8-0

set M := {12..23};

param X{M};
param Y{M};

param b := X[13];
param a := -b;
param diff := 2*b;

param t{i in M, j in 0..14} := if (j==0) then 1.0 else
			       if (j==1) then (2*X[i]-a-b)/diff else
			       2*t[i,j-1]*(2*X[i]-a-b)/diff - t[i,j-2];
var A0 := 1.0;
var A2 := 1.0;
var A4 := 1.0;
var A6 := 1.0;
var A8 := 1.0;
var A10 := 1.0;
var B >= 0, := 1.0;
var C >= 1e-5, := 1.0;

minimize palmer:
	sum {m in M} (Y[m] -
(A0*t[m,0] + A2*t[m,2] + A4*t[m,4] + A6*t[m,6] + A8*t[m,8] + A10*t[m,10] + B/(C+X[m]^2)) )^2;

data;
param X:= 
12                  0.000000
13                  1.570796
14                  1.396263
15                  1.308997
16                  1.221730
17                  1.125835
18                  1.047198
19                  0.872665
20                  0.698132
21                  0.523599
22                  0.349066
23                  0.174533
;
param Y:=
12                 83.57418
13                 81.007654
14                 18.983286
15                  8.051067
16                  2.044762
17                  0.000000
18                  1.170451
19                 10.479881
20                 25.785001
21                 44.126844
22                 62.822177
23                 77.719674
;

display palmer;
solve;

display A0,A2,A4,A6,B,C;
