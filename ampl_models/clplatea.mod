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
#   J. Nocedal,
#   "Solving large nonlinear systems of equations arising in mechanics",
#   Proceedings of the Cocoyoc Numerical Analysis Conference, Mexico,
#   pp. 132-141, 1981.

#   SIF input: Ph. Toint, Dec 1989.

#   classification OXR2-MN-V-0

param p := 71;
param wght := -0.1;
param hp2 := 0.5*p^2;

var x{1..p, 1..p} := 0.0;

minimize f:
	sum {i in 2..p, j in 2..p} (
	0.5*(x[i,j]-x[i,j-1])^2+
	0.5*(x[i,j]-x[i-1,j])^2+
	hp2*(x[i,j]-x[i,j-1])^4+
	hp2*(x[i,j]-x[i-1,j])^4
	) + (wght*x[p,p]);

fix {j in 1..p} x[1,j] := 0.0;	

solve;
display f;
display x;
