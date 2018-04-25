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
#   A. Neumaier, private communication, 1991.

#   SIF input: Ph. L. Toint, May 1993.
#              minor correction by Ph. Shott, Jan 1995.

#   classification NOR2-AN-V-V

param N:=100;
param al{i in 1..N} := if (i < N) then (i+1)/N else 1;
param a{k in 1..N} := if (k=1) then sum {i in 1..N} al[i] else 
			sum {i in 1..N} exp(k*log(al[i]));
var x{i in 1..N} := (i-1)/N;

minimize f:0;
subject to cons1:
	(sum {i in 1..N} x[i] - a[1])^2 = 0;
subject to cons2{k in 2..N}:
	(sum {i in 1..N} x[i]^k - a[k])^2 = 0;
subject to cons3{i in 2..N}:
	x[i]-x[i-1] >= 0;

solve;
display f;
display x;
