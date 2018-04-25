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

#   Source: problem 21 in
#   W. Hock and K. Schittkowski,
#   "Test examples for nonlinear programming codes",
#   Lectures Notes in Economics and Mathematical Systems 187, Springer
#   Verlag, Heidelberg, 1981.

#   SIF input: A.R. Conn, April 1990

#   classification SLR2-AN-7-1

var x{i in 1..7} := if (i <= 2) then -1.0;

minimize f:
	-100+0.01*(x[1]^2+x[3]^2+x[5]^2+x[6]^2) + (x[2]^2+x[4]^2+x[7]^2);
subject to cons1:
	10*x[1]-x[2]-10 >= 0;
subject to cons2:
	2 <= x[1] <= 50;
subject to cons3:
	-50 <= x[2] <=50;
subject to cons4:
	x[3] <= 50;
subject to cons5:
	2 <= x[4];
subject to cons6:
	x[6] <= 0;
subject to cons7:
	0 <= x[7];

solve; display f; display x;
