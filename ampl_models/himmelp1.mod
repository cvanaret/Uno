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
#   B.N. Pshenichnyj
#   "The Linearization Method for Constrained Optimization",
#   Springer Verlag, SCM Series 22, Heidelberg, 1994

#   SIF input: Ph. Toint, December 1994.

#   classification OBR2-AN-2-0

param b{1..20};
param xinit{i in 1..2};
var x{i in 1..2} := xinit[i];

minimize f:
	-b[2]*x[1]-b[6]*x[2]-b[1]- ( b[3]*x[1]^2+b[4]*x[1]^3+b[5]*x[1]^4+ 
x[2]*(b[7]*x[1]+b[8]*x[1]^2+b[9]*x[1]^3+b[10]*x[1]^4) +
b[11]*x[2]^2+b[12]*x[2]^3+b[13]*x[2]^4 + b[14]/(1+x[2]) + (b[18]*x[1] +
b[15]*x[1]^2 + b[16]*x[1]^3)*x[2]^2 + (b[17]*x[1]^3+b[19]*x[1])*x[2]^3 +
b[20]*exp(0.0005*x[1]*x[2]) );

subject to cons1:
	0.0 <= x[1] <= 95.0;
subject to cons2:
	0.0 <= x[2] <= 75.0;

data;
param xinit:= 1 95 2 10;
param b:=
1	75.1963666677
2	-3.8112755343
3	0.1269366345
4	-0.0020567665
5	0.103450d-4
6	-6.8306567613
7	.0302344793
8	-0.0012813448
9	0.352599d-4
10	-0.2266d-6
11	0.2564581253
12	-.003460403
13	0.135139d-4
14	-28.1064434908
15	-0.52375d-5
16	-0.63d-8
17	0.7d-9
18	0.0003405462
19	-0.16638d-5
20	-2.8673112392;

solve;
display f;
display x;
