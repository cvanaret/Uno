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

#   Source:  The first problem given by
#   R. Fletcher,
#   "An optimal positive definite update for sparse Hessian matrices"
#   Numerical Analysis report NA/145, University of Dundee, 1992.

#   N.B. This formulation is incorrect. See FLETCBV2.SIF for
#        the correct version.

#   SIF input: Nick Gould, Oct 1992.

#   classification OUR2-AN-V-0

param n := 10000;
param kappa := 1.0;
param objscale := 1.0D+0;
param h := 1/(n+1);
param p := 1/objscale;

var x{i in 1..n} := i*h;

minimize f:
	0.5*p*(x[1])^2 +
	sum {i in 1..n-1} 0.5*p*(x[i]-x[i+1])^2 +
	0.5*p*(x[n])^2 +
	sum {i in 1..n} (p*(-1-2/h^2)*x[i]) +
	sum {i in 1..n} (-kappa*p*cos(x[i])/h^2);

solve;
display f;
display x;

