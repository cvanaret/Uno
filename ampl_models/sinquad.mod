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
#   N. Gould, private communication.

#   SIF input: N. Gould, Dec 1989.

#   classification OUR2-AY-V-0

param N:=10000;

var x{1..N} := 0.1;

minimize f:
	(x[1]-1)^4
	+ sum {i in 2..N-1} (sin(x[i]-x[N])-x[1]^2+x[i]^2)^2
	+ (x[N]^2-x[1]^2)^2
	;

solve;
display f;
display x;
