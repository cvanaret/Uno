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

#   Source: problem 5 in
#   S. Nash,
#   "Newton-type minimization via the Lanczos process",
#   SIAM J. Num. Anal. 21, 1984, 770-788.

#   SIF input: Nick Gould, Oct 1992.
#              minor correction by Ph. Shott, Jan 1995.

#   classification SUR2-AN-V-0

param n := 500;

var x{1..n} := 1/(n+1);

minimize f:
	1.0 +
	sum {i in 2..n} 100*(x[i]-x[i-1]^2)^2 +
	sum {i in 2..n} (x[i]-1.0)^2;

solve;
display f;
display x;
