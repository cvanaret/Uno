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

#   Source: problem 57 in
#   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
#   "Performance of a multi-frontal scheme for partially separable
#   optimization"
#   Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.

#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AN-V-0

param N := 10000;
var x{i in 1..N} := if (i mod 2 = 0) then -1 else 1;

minimize f:
	sum {i in 1..N-2} (x[i]+x[i+1]+x[N])^4 + (x[1]-x[2])^2 + (x[N-1]+x[N])^2;

solve; display f; display x; 
