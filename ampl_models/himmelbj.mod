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

#   Source: problem 6 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.

#   SIF input: Ph. Toint, March 1991.

#   classification OLR2-MY-45-14

param NSETS:=7;
param NS{1..7};
param NEQ := 6;
var X{k in 1..NSETS, j in 1..NS[k]} >= 1.0D-12, := 0.1;

minimize f:
 	X[1,2]*-7.69+ X[1,3]*-11.52+
 	X[1,4]*-36.60+ X[2,1]*-10.94+
 	X[2,8]*2.5966+ X[2,9]*-39.39+
 	X[2,10]*-21.35+ X[2,11]*-32.84+
 	X[2,12]*6.26+ X[3,1]*10.45+
 	X[3,3]*-0.5+ X[3,7]*2.2435+
 	X[3,9]*-39.39 +X[3,10]*-21.49+
 	X[3,11]*-32.84+ X[3,12]*6.12+
 	X[3,15]*-1.9028 +X[3,16]*-2.8889+
 	X[3,17]*-3.3622 +X[3,18]*-7.4854+
 	X[4,1]*-15.639 +X[4,3]*21.81+
 	X[5,1]*-16.79 +X[5,3]*18.9779+
 	X[6,2]*11.959+X[7,2]*12.899+
	sum {k in 1..NSETS, j in 1..NS[k]} X[k,j]*log(X[k,j]) -
	sum {j in 1..4} X[1,j]*log(sum {i in 1..4} X[1,i]) - 
	sum {j in 1..13} X[2,j]*log(sum {i in 1..13} X[2,i]) - 
	sum {j in 1..18} X[3,j]*log(sum {i in 1..18} X[3,i]) - 
	sum {j in 1..3} X[4,j]*log(sum {i in 1..3} X[4,i]) - 
	sum {j in 1..3} X[5,j]*log(sum {i in 1..3} X[5,i]) - 
	sum {j in 1..2} X[6,j]*log(sum {i in 1..2} X[6,i]) - 
	sum {j in 1..2} X[7,j]*log(sum {i in 1..2} X[7,i])  
	;

subject to cons1:
	X[1,1]+X[2,1]+X[3,1]+X[3,15]+X[3,16]*2+X[3,17]*3+X[3,18]*4 =
	0.652981;
subject to cons2:
	X[1,2]+X[2,2]+X[2,10]+X[2,11]+X[2,12]+X[3,2]+X[3,10]+X[3,11]+X[3,12] +
	X[6,2]+X[7,2] = 0.281941;
subject to cons3:
	X[1,3]+X[2,3]+X[3,3] = 3.705233;
subject to cons4:
	X[1,4]+X[2,4]+X[2,9]+X[2,11]-X[2,12]+X[3,4]+X[3,9]+X[3,11]-X[3,12]+X[4,1]
	-X[4,3]+X[5,1]-X[5,3]-X[6,2]-X[7,2]= 47.00022;
subject to cons5:
	X[1,4]+X[2,5]+X[2,9]+X[2,10]+X[2,11]+X[2,12]+X[3,5]+X[3,9]+
	X[3,10]+X[3,11]+X[3,12] = 47.02972;
subject to cons6:
	X[2,6]+X[3,6] = 0.08005;
subject to cons7:
	X[2,7]+X[3,7] = 0.08813;
subject to cons8:
	X[2,8]+X[3,8] = 0.04829;
subject to cons9:
	X[3,14]+X[3,15]+X[3,16]+X[3,17]+X[3,18] = 0.0022725;
subject to cons10:
	X[2,4]-X[2,5]-X[2,6]+X[2,7]+X[2,8]+X[2,10]-2*X[2,12]-X[2,13] -
	4*X[3,14]-3*X[3,15]-2*X[3,16]-X[3,17] = 0;
subject to cons11:
	-X[3,15]-2*X[3,16]-3*X[3,17]-4*X[3,18]+X[4,1]+X[4,2]+X[4,3] =0;
subject to cons12:
	X[5,1]+X[5,2]+X[5,3] = 0;
subject to cons13:
	-4*X[4,3]+X[6,1]+X[6,2] = 0;
subject to cons14:
	-4*X[5,3]+X[7,1]+X[7,2] = 0;
subject to cons15:
	X[2,13] = 0.0155;
subject to cons16:
	X[3,13] = 0.0211275;

data;
param NS:=
1	4
2	13
3	18
4	3
5	3
6	2
7	2;

option presolve_eps 0.000001;
solve;
display f;
display X;
