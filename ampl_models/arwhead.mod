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

#   Source: Problem 55 in
#   A.R. Conn, N.I.M. Gould, M. Lescrenier and Ph.L. Toint,
#   "Performance of a multifrontal scheme for partially separable
#   optimization",
#   Report 88/4, Dept of Mathematics, FUNDP (Namur, B), 1988.

#   SIF input: Ph. Toint, Dec 1989.

#   classification OUR2-AN-V-0

param N:=5000;
var x{1..N} := 1.0;

minimize f:
	sum {i in 1..N-1} (-4*x[i]+3.0) + sum {i in 1..N-1} (x[i]^2+x[N]^2)^2;

solve;
display f;
display x;
