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

#   Source: problem 5 (p. 89) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-2-0

param N:=2;
param xinit{1..2};
var x{i in 1..N} := xinit[i];

minimize f:
	(x[1]-1.0)^2+sum { i in 2..N} 100*(x[i]-x[i-1]^3)^2;

data;
param xinit:= 1 -1.2 2 1.0;

solve; display f; display x;
