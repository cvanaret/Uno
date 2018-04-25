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

#   Source:  problem 7 (p. 102) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification QUR2-AN-2-0

param xinit{1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	(-56*x[1]-256*x[2]+991+16*x[1]^2+16*x[2]^2-8*x[1]*x[2])/15;

data;
param xinit:= 1 3.0 2 8.0;

solve;
display f;
display x;
