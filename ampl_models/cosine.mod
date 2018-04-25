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

#   SIF input: N. Gould, Jan 1996

#   classification OUR2-AN-V-0

param N:=10000;
var x{1..N} := 1.0;

minimize f:
	sum {i in 1..N-1} cos(-0.5*x[i+1]+x[i]^2);

solve; display f; display x;
