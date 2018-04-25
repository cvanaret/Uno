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
#   Hans Mittelmann, private communication.

#   SIF input: N. Gould, Jan 1998

#   classification OOR2-AY-V-0

param n:=250;
var x{1..n} := 1.0;

minimize f:
	sum {i in 1..n-1, j in i+1..n} exp(x[j]-x[i]);
subject to cons1{i in 1..n-1}:
	x[i]+x[i+1] >= 0.0;
subject to cons2{i in 1..n-1}:
	exp(x[i])+exp(x[i+1]) <= 20.0;

solve;
display f;
display x;
