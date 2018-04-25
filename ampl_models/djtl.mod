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

#   Source: modified version of problem 19 in
#   W. Hock and K. Schittkowski,
#   "Test examples for nonlinear programming codes",
#   Lectures Notes in Economics and Mathematical Systems 187, Springer
#   Verlag, Heidelberg, 1981.
#   that is meant to simulate the Lagrangian barrier objective function
#   for particular values of the shifts and multipliers

#   SIF input: A.R. Conn August 1993

#   classification OUR2-AN-2-0

#rvdb: this is stupid
var x1 := 15;
var x2 := -1;

minimize f:
	(x1-10)^3+(x2-20)^3
	+ if (-(x1-5)^2-(x2-5)^2+200+1 <= 0.0) then
	  1D10*(-(x1-5)^2-(x2-5)^2+200)^2 else
	  -log(-(x1-5)^2-(x2-5)^2+200+1)
	+ if ((x1-5)^2+(x2-5)^2-100+1 <= 0.0) then
	  1D10*((x1-5)^2+(x2-5)^2-100)^2 else
	  -log((x1-5)^2+(x2-5)^2-100+1)
	+ if ((x2-5)^2+(x1-6)^2+1 <= 0.0) then
	  1D10*((x2-5)^2+(x1-6)^2)^2 else
	  -log((x2-5)^2+(x1-6)^2+1)
	+ if (-(x2-5)^2-(x1-6)^2+82.81+1 <= 0.0) then
	  1D10*(-(x2-5)^2-(x1-6)^2+82.81)^2 else
	  -log(-(x2-5)^2-(x1-6)^2+82.81+1)
	+ if (100-x1+1 <= 0.0) then
	  1D10*(100-x1)^2 else
	  -log(100-x1+1)
	+ if (x1-13+1 <= 0.0) then
	  1D10*(x1-13)^2 else
	  -log(x1-13+1)
	+ if (100-x2+1 <= 0.0) then
	  1D10*(100-x2)^2 else
	  -log(100-x2+1)
	+ if (x2+1<= 0.0) then
	  1D10*(x2)^2 else
	  -log(x2+1)
;
solve;
display f;
display x1, x2;
