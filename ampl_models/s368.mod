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
#   P. Wolfe,
#   "Explicit solution of an optimization problem",
#   Mathematical Programming 2, 258-260, 1972.

#   SIF input: Nick Gould, Oct 1992.

#   See also Schittkowski #368 (for N = 8)

#   classification OBR2-MN-V-0

param n := 100;

var x{1..n} <= 1.0, >= 0.0;

minimize f:
	sum {i in 1..n, j in 1..n} (-x[i]^2*x[j]^4) +
	sum {i in 1..n, j in 1..n} (x[i]^3*x[j]^3);

solve;
