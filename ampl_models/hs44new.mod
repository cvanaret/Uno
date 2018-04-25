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

#   Source: problem 44 in
#   W. Hock and K. Schittkowski,
#   "Test examples for nonlinear programming codes",
#   Lectures Notes in Economics and Mathematical Systems 187, Springer
#   Verlag, Heidelberg, 1981.

#   SIF input: Ph.L. Toint, October 1990.

#   classification QLR2-AN-4-6

param N:=4;
var x{1..N} >= 0, := 1.0;

minimize f:
	x[1]-x[2]-x[3]-x[1]*x[3]+x[1]*x[4]+x[2]*x[3]-x[2]*x[4];
subject to cons1:
	-x[1]-x[2]+8.0 >= 0;
subject to cons2:
	-4*x[1]-x[2]+12 >= 0;
subject to cons3:
	-3*x[1]-4*x[2]+12 >= 0;
subject to cons4:
	-2*x[3]-x[4]+8 >= 0;
subject to cons5:
	-x[3]-2*x[4]+8 >= 0;
subject to cons6:
	-x[3]-x[4]+5 >= 0;

solve; display f; display x;
