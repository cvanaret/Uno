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
#   D. Shanno,
#   " On Variable Metric Methods for Sparse Hessians II: the New
#   Method",
#   MIS Tech report 27, University of Arizona (Tucson, UK), 1978.

#   See also Buckley #37 (p. 76) and Toint #15.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-V-0

param N:=10000;
var x{1..N} := -1;

minimize f:
	(x[1]-1)^2
	+ sum {i in 2..N} 100*(x[1]-x[i-1]^2)^2;

solve;
display f;
display x;
