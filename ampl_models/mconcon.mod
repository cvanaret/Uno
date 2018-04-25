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

#   classification LOI2-MN-15-11

param N := 7;
param M := 4;

var P{1..N} := 965;
var Q{i in 1..M} := if (i>2) then -100 else 100;
var F{1..M} := 1000;

minimize f:
	sum {i in 1..N} (-P[i]);
subject to cons1_1:
	P[1]*abs(P[1])-P[2]*abs(P[2])-0.597053452*Q[1]*abs(Q[1])^0.8539 = 0;
subject to cons1_2:
	P[3]*abs(P[3])-P[4]*abs(P[4])-0.597053452*Q[2]*abs(Q[2])^0.8539 = 0;
subject to cons1_3:
	P[4]*abs(P[4])-P[5]*abs(P[5])-0.597053452*Q[3]*abs(Q[3])^0.8539 = 0;
subject to cons1_4:
	P[6]*abs(P[6])-P[7]*abs(P[7])-0.597053452*Q[4]*abs(Q[4])^0.8539 = 0;
subject to cons2:
	Q[1]-F[3] = 0;
subject to cons3:
	-Q[1]+F[1] = 0;
subject to cons4:
	Q[2]-F[1] = 0;
subject to cons5:
	-Q[2]+Q[3]+1000 = 0;
subject to cons6:
	-Q[3]-F[2] = 0;
subject to cons7:
	Q[4]+F[2] = 0;
subject to cons8:
	-Q[4]-F[4] = 0;
subject to cons9:
	P[3] <= 904.73;
subject to cons10:
	P[5] <= 904.73;
subject to cons11:
	P[1] <= 914.73;
subject to cons12:
	P[7] <= 914.73;
subject to cons13:
	F[4] <= 400;

solve;
display f;
display P,Q,F;
