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
#   H. Zhang and X. Wang,
#   "Optimal sensor placement",
#   SIAM Review, vol. 35, p. 641, 1993.

#   SIF input: Nick Gould, June 1994

#   classification SUR2-AN-V-0

param N:=1000;
var theta{i in 1..N} := i/N;

minimize f:
	sum {i in 1..N, j in 1..N}
-(sin(theta[i])*sin(theta[j])*sin(theta[i]-theta[j]))^2;

solve;
display f;
display theta;
