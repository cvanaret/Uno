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

#   Source: problem 35 in
#   W. Hock and K. Schittkowski,
#   "Test examples for nonlinear programming codes",
#   Lectures Notes in Economics and Mathematical Systems 187, Springer
#   Verlag, Heidelberg, 1981.

#   SIF input: A.R. Conn, April 1990

#   classification QLR2-AN-3-1

var x{1..3} := 0.5, >= 0.0;

minimize f:
	-8.0*x[1]-6.0*x[2]-4.0*x[3]+9.0+2.0*x[1]^2+2.0*x[2]^2+x[3]^2+2.0*x[1]*x[2]+2.0*x[1]*x[3];
subject to cons1:
	-x[1]-x[2]-2.0*x[3]+3.0 >= 0;
subject to cons2:
	x[2] = 0.5;

solve; display f; display x;
