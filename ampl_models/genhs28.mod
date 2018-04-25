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

#   Source: a multi-dimensional extension of problem 28 in
#   W. Hock and K. Schittkowski,
#   "Test examples for nonlinear programming codes",
#   Lectures Notes in Economics and Mathematical Systems 187, Springer
#   Verlag, Heidelberg, 1981.

#   SIF input: Nick Gould, December 1991
#              minor correction by Ph. Shott, Jan 1995.

#   classification QLR2-AY-10-8

param N:=10;
var x{i in 1..N} := if (i=1) then -4.0 else 1.0;

minimize f:
	sum {i in 1..N-1} (x[i]+x[i+1])^2;

subject to cons{i in 1..N-2}:
	-1.0+x[i]+2*x[i+1]+3*x[i+2] = 0;

solve;
display f;
display x;
