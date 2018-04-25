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

#   Source:  problem 208 (p. 56) in
#   A.R. Buckley,
#   "Test functions for unconstrained minimization",
#   TR 1989CS-3, Mathematics, statistics and computing centre,
#   Dalhousie University, Halifax (CDN), 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification NQR2-AN-2-2

var x{1..2} := 0.5;

minimize f: 0;
subject to cons:
	(x[1]-0.1136*(x[1]+3.0*x[2])*(1-x[1])) = 0;
subject to cons2:
	(x[2]+7.5*(2.0*x[1]-x[2])*(1-x[2])) = 0;

solve; display f; display x;
