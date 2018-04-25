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

#   Source: problem 32 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.

#   See Buckley#76 (p. 66)

#   SIF input: Ph. Toint, Dec 1989.

#   classification SUR2-AN-4-0

param a{1..7};
param b{1..7};
param xinit{1..4};
var x{i in 1..4} := xinit[i];

minimize f:
	10000*sum {i in 1..7} (-1 +
(x[1]^2+a[i]*x[2]^2+a[i]^2*x[3]^2)/(b[i]*(1+a[i]*x[4]^2)))^2;

data;
param:
	a		b:=
1	0.0		7.391
2	0.000428	11.18
3	0.001000	16.44
4	0.001610	16.20
5	0.002090	22.20
6	0.003480	24.02
7	0.005250	31.32;

param xinit:=
1	2.7
2	90.0
3	1500.0
4	10.0;

solve; display f; display x;
