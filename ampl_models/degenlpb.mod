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
#   T.C.T. Kotiah and D.I. Steinberg,
#   "Occurences of cycling and other phenomena arising in a class of
#   linear programming models",
#   Communications of the ACM, vol. 20, pp. 107-112, 1977.

#   SIF input: Ph. Toint, Aug 1990.

#   classification LLR2-AN-20-15

param N:=20;
param M:=15;
var x{1..N} := 1.0, <= 1.0, >= 0.0;

minimize f:
	-1*(0.01*x[2]+33.333*x[3]+100.0*x[4]+0.01*x[5]+33.343*x[6]+100.01*x[7]+33.333*x[8]+133.33*x[9]+100.0*x[10]);
subject to cons1:
	-0.70785+x[1]+2*x[2]+2*x[3]+2*x[4]+x[5]+2*x[6]+2*x[7]+x[8]+2*x[9]+x[10] = 0;
subject to cons2:
	0.326*x[1]-101*x[2]+200*x[5]+0.06*x[6]+0.02*x[7] = 0;
subject to cons3:
	0.0066667*x[1]-1.03*x[3]+200*x[6]+0.06*x[8]+0.02*x[9] = 0;
subject to cons4:
	0.00066667*x[1]-1.01*x[4]+200*x[7]+0.06*x[9]+0.02*x[10] = 0;
subject to cons5:
	0.978*x[2]-201*x[5]+100*x[11]+0.03*x[12]+0.01*x[13] = 0;
subject to cons6:
	0.01*x[2]+0.489*x[3]-101.03*x[6]+100*x[12]+0.03*x[14]+0.01*x[15] = 0;
subject to cons7:
	0.001*x[2]+0.489*x[4]-101.03*x[7]+100*x[13]+0.03*x[15]+0.01*x[16] = 0;
subject to cons8:
	0.001*x[3]+0.01*x[4]-1.04*x[9]+100*x[15]+0.03*x[18]+0.01*x[19] = 0;
subject to cons9:
	0.02*x[3]-1.06*x[8]+100*x[14]+0.03*x[17]+0.01*x[19] = 0;
subject to cons10:
	0.002*x[4]-1.02*x[10]+100*x[16]+0.03*x[19]+0.01*x[20] = 0;
subject to cons11:
	-2.5742D-6*x[11]+0.00252*x[13]-0.61975*x[16]+1.03*x[20] = 0;
subject to cons12:
	-0.00257*x[11] + 0.25221*x[12] - 6.2*x[14] + 1.09*x[17] = 0;
subject to cons13:
	0.00629*x[11]-0.20555*x[12]-4.1106*x[13]+101.04*x[15]+505.1*x[16]-256.72*x[19] = 0;
subject to cons14:
	0.00841*x[12]-0.08406*x[13]-0.20667*x[14]+20.658*x[16]+1.07*x[18]-10.5*x[19] = 0;
subject to cons15:
	-x[1]+300*x[2]+0.09*x[3]+0.03*x[4] = 0;

solve; display f; display x;
