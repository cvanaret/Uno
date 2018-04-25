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
#   M.J.D. Powell,
#   "A tolerant algorithm for linearly constrained optimization
#   calculations"'
#   Mathematical Programming 45(3), pp.561--562, 1989.

#   SDIF input: Ph. Toint and N. Gould, May 1990.

#   classification OLR2-AN-5-22

param R := 11;
param T{i in 1..R} := 5*(i-1)/(R-1);
param ET{i in 1..R} := exp(T[i]);
param pinit{0..2};

var P{i in 0..2} := pinit[i];
var Q{1..2} := 0.0;

minimize f:
	sum {i in 1..R} (
	(P[0]+P[1]*T[i]+P[2]*T[i]^2)/
	(ET[i]*(1+Q[1]*(T[i]-5)+Q[2]*(T[i]-5)^2))
	-1
	)^2;
subject to cons1{i in 1..R}:
	P[0]+P[1]*T[i]+P[2]*T[i]^2 - (T[i]-5)*ET[i]*Q[1] -
		(T[i]-5)^2*ET[i]*Q[2] -ET[i]>= 0;
subject to cons2{i in 1..R}:
	(T[i]-5)*Q[1] + (T[i]-5)^2*Q[2]+0.99999 >= 0;

data;
param pinit:= 0 1 1 1 2 6;

solve;
display f;
display P,Q;
