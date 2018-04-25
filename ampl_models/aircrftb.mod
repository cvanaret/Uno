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

#   Source: problem 9 in
#   J.J. More',"A collection of nonlinear model problems"
#   Proceedings of the AMS-SIAM Summer Seminar on the Computational
#   Solution of Nonlinear Systems of Equations, Colorado, 1988.
#   Argonne National Laboratory MCS-P60-0289, 1989.

#   SIF input: Ph. Toint, Dec 1989.

#   classification SXR2-RN-8-0

param m;
param n;

set N := 1..n;
set M := 1..m;

param A{M,N};
var X{N};
var AX{i in M} = sum{j in N} (A[i,j]*X[j]);

var P1=(-0.727*X[2]*X[3])+(8.39*X[3]*X[4])-(684.4*X[4]*X[5])+(63.5*X[4]*X[2]);
var P2 = (0.949*X[1]*X[3])+(0.173*X[1]*X[5]);
var P3 = (-0.716*X[1]*X[2])-(1.578*X[1]*X[4])+(1.132*X[4]*X[2]);
var P4 = -1*X[1]*X[5];
var P5 = X[1]*X[4];

minimize L2force:
	(AX[1]+P1)^2 +(AX[2]+P2)^2+(AX[3]+P3)^2+(AX[4]+P4)^2+(AX[5]+P5)^2;

subject to elevator:
	X[6] = -0.05;

subject to aileron:
	X[7] = 0.1;

subject to rudder:
	X[8] = 0.0;

data;
param m := 5;
param n := 8;

param A:
	1	2	3	4	5	6	7	8:=
1	-3.933	0.107	0.126	0	-9.99	0	-45.83	-7.64
2	0	-0.987	0	-22.95	0	-28.37	0	0
3	0.002	0	-0.235	0	5.67	0	-0.921	-6.51
4	0	1.0	0	-1.0	0	-0.168	0	0
5	0	0	-1.0	0	-0.196	0	-0.0071	0;

var X =
1 0
2 0
3 0
4 0
5 0
6 -0.05
7 0.1
8 0.0;


solve;

display X;
