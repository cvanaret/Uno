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

#   Source: from problem 20 in
#   D.H. Himmelblau,
#   "Applied nonlinear programming",
#   McGraw-Hill, New-York, 1972.

#   SIF input: Ph. Toint, March 1991.

#   classification LOR2-MN-24-14

param F:=142.22471;
param B{1..24};
param C{1..12};
param D{1..12};
var x{1..24}:=0.04, >= 0.0;

minimize f:
	0.0693*x[1]+0.0577*x[2]+0.05*x[3]+0.2*x[4]+0.26*x[5]+0.55*x[6] +
0.06*x[7]+0.1*x[8]+0.12*x[9]+0.18*x[10]+0.1*x[11]+0.09*x[12]+0.0693*x[13]
+0.0577*x[14]+0.05*x[15]+0.2*x[16]+0.26*x[17]+0.55*x[18]+0.06*x[19] 
+0.1*x[20]+0.12*x[21]+0.18*x[22]+0.1*x[23]+0.09*x[24];

subject to cons1{i in 1..12}:
	sum {j in 1..12} x[i+12]*x[j]*40*B[i]/B[j] - 
	sum {j in 13..24} x[i]*x[j]*C[i]*B[i+12]/B[j] = 0 ;
subject to cons2:
	sum {i in 1..24} x[i]=1.0;
subject to cons3:
	sum {i in 1..12} (x[i]/D[i]+x[i+12]*F/B[i+12])=1.671;
data;
param B:=
1	44.094
2	58.12
3	58.12
4	137.4
5	120.9
6	170.9
7	62.501
8	84.94
9	133.425
10	82.507
11	46.07
12	60.097
13	44.094
14	58.12
15	58.12
16	137.4
17	120.9
18	170.9
19	62.501
20	84.94
21	133.425
22	82.507
23	46.07
24	60.097;

param:
	C	D:=
1	123.7	123.7
2	31.7	31.7
3	45.7	45.7
4	14.7	14.7
5	84.7	84.7
6	27.7	27.7
7	49.7	49.7
8	7.1	7.1
9	2.1	2.1
10	17.7	17.7
11	0.85	0.85
12	0.64	0.64;

solve;
display f;
display x;
