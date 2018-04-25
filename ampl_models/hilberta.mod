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
#   K. Schittkowski,
#   "More Test Examples for Nonlinear Programming Codes",
#   Springer Verlag, Heidelberg, 1987.

#   See also Buckley#19 (p. 59)

#   SIF input: Ph. Toint, Dec 1989.

#   classification QUR2-AN-V-0

param N := 10;

var x{1..N};
param A{i in 1..N, j in 1..N} := 1/(i+j-1);

minimize f:
sum {i in 1..N} x[i]*(sum {j in 1..N}A[i,j]*x[j]);

data;
var x:=
1	-4
2	-2;

solve;
display f;
display x;
