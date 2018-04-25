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
#   M.  Palmer, Edinburgh, private communication.

#   SIF input: Nick Gould, 1992.

#   classification SBR2-RN-8-0

set M := {12..24};

param X{M};
param Y{M};

var A0 := 1.0;
var A2 := 1.0;
var A4 := 1.0;
var A6 := 1.0;
var A8 := 1.0;
var A10 := 1.0;
var K >= 0.0, := 1.0;
var L := 1.0;

minimize palmer:
	sum {m in M} (Y[m] -
(A0 + A2*X[m]^2 + A4*X[m]^4 + A6*X[m]^6 + A8*X[m]^8 + A10*X[m]^10 + L*exp(-K*X[m]^2)) )^2;

data;
param X:=
12                 0.000000
13                 0.139626
14                 0.261799
15                 0.436332
16                 0.565245
17                 0.512942
18                 0.610865
19                 0.785398
20                 0.959931
21                 1.134464
22                 1.308997
23                 1.483530
24                 1.658063
;
param Y:=
12                   4.419446
13                   3.564931
14                   2.139067
15                   0.404686
16                   0.000000
17                   0.035152
18                   0.146813
19                   2.718058
20                   9.474417
21                  26.132221
22                  41.451561
23                  72.283164
24                 117.630959
;

solve;

display A0,A2,A4,A6,A8,A10,K,L;
