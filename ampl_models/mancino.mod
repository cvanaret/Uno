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
#   E. Spedicato,
#   "Computational experience with quasi-Newton algorithms for
#   minimization problems of moderate size",
#   Report N-175, CISE, Milano, 1975.

#   See also Buckley #51 (p. 72), Schittkowski #391 (for N = 30)

#   SIF input: Ph. Toint, Dec 1989.
#              correction by Ph. Shott, January, 1995.

#   classification SUR2-AN-V-0

param N := 100;

var x{i in 1..N} := -8.710996D-4*((i-50)^3 + 
		    sum {j in 1..N} sqrt(i/j)*((sin(log(sqrt(i/j))))^5+ 
		    (cos(log(sqrt(i/j))))^5));
var v{i in 1..N, j in 1..N} = sqrt (x[i]^2 +i/j);
var alpha{i in 1..N} = 1400*x[i] + (i-50)^3 + 
		       sum {j in 1..N} v[i,j]*((sin(log(v[i,j])))^5 + 
		       (cos(log(v[i,j])))^5);

minimize f:
sum {i in 1..N} alpha[i]^2;

solve;
display f;
display x;
