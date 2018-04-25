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

#   Source: a modification (by Ph. Toint) of
#   L.C.W. Dixon, personnal communication, Jan 1991.

#   SIF input: Ph. Toint, Feb 1991.

#   classification SOR2-AN-V-V

param n := 100;

var x{i in 1..n} >= 1D-8, := if ( (i mod 2) == 1) then (i+1) else 1/(i+1);

minimize f:
	sum {i in 1..n-3} (100*(x[i+1]-x[i]^2)^2 + (x[i]-1.0)^2 + 90*(x[i+3]-x[i+2]^2)^2 + (x[i+2]-1.0)^2 + 10.1*(x[i+1]-1.0)^2 + 10.1*(x[i+3]-1.0)^2 + 19.8*((x[i+1]-1.0)*(x[i+3]-1.0)));

subject to cons1{i in 2..n by 2}:
	sum {j in 1..i} log(x[j]) = 0;

solve;
