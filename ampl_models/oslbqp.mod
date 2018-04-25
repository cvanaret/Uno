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

#   Source: Simple convex QP in OSL manual



#  Minimize   x1 + 2x5 - x8 + 1/2(x1##2 + x2##2 + x3##2 + x4##2
#                            + x5##2 + x6##2 + x7##2 + x8##2)
#   Subject to:
#   2.5 <= x1
#     0 <= x2 <= 4.1
#     0 <= x3
#     0 <= x4
#   0.5 <= x5 <= 4.0
#     0 <= x6
#     0 <= x7
#     0 <= x8 <= 4.3


#   SIF input: A.R. Conn, December 1992

#   classification QBR2-AN-8-0

var x{1..8} := 0.5;

minimize f:
	x[1]+2*x[5]-x[8]+0.5*sum{i in 1..8} x[i]^2;
subject to cons1:
	2.5 <= x[1];
subject to cons2:
	0 <= x[2] <= 4.1;
subject to cons3:
	0 <= x[3];
subject to cons4:
	0 <= x[4];
subject to cons5:
	0.5 <= x[5] <= 4.0;
subject to cons6:
	0 <= x[6];
subject to cons7:
	0 <= x[7];
subject to cons8:
	0 <= x[8] <= 4.3;

solve; display f; display x;
