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

#   classification LQI2-RN-157-134

param Region{1..21};
param lopi{1..48};
param hipi{1..48};
param hisupply{1..7};
param loflow;
param hiflow;
param loprod{1..5};
param hiprod{1..5};

var Prod{i in 1..5}>= loprod[i], <= hiprod[i], := 0.0;
var Pi{i in 1..48}>= lopi[i], <= hipi[i], := 0.0;
var Flow{1..60}>=loflow, <=hiflow;
var Supply{i in 1..7}<= hisupply[i], := 0.0;
var Demand{1..37} := 0.0;

minimize cost:
Prod[1]+Prod[2]+Prod[3]+Prod[4]+Prod[5];
subject to node1:
Flow[1]-Flow[16]-Flow[56]=0;
subject to node2:
-1*Flow[1]+Demand[7]=0;
subject to node3:
Flow[2]+Flow[55]-Supply[4]=0;
subject to node4:
-1*Flow[2]+Demand[5]=0;
subject to node5:
Flow[3]+Demand[6]=0;
subject to node6:
-1*Flow[3]+Flow[4]+Demand[8]=0;
subject to node7:
-1*Flow[4]-Flow[35]+Flow[36]-Flow[47]+Flow[48]=0;
subject to node8:
Flow[5]+Flow[9]+Demand[21]-Supply[5]=0;
subject to node9:
-1*Flow[5]+Flow[6]+Demand[22]=0;
subject to node10:
-1*Flow[6]+Flow[7]-Flow[8]=0;
subject to node11:
-1*Flow[7]+Demand[9]=0;
subject to node12:
Flow[8]-Flow[10]+Flow[11]=0;
subject to node13:
-1*Flow[9]+Flow[10]+Demand[23]=0;
subject to node14:
-1*Flow[11]+Flow[12]+Demand[24]=0;
subject to node15:
-1*Flow[12]-Flow[15]+Flow[16]=0;
subject to node16:
Flow[13]-Flow[36]-Flow[48]-Flow[54]+Demand[10]=0;
subject to node17:
-1*Flow[13]+Flow[14]+Demand[25]=0;
subject to node18:
-1*Flow[14]+Flow[15]+Demand[26]=0;
subject to node19:
Flow[17]+Flow[27]+Flow[40]=0;
subject to node20:
-1*Flow[17]+Flow[18]+Demand[27]=0;
subject to node21:
-1*Flow[18]+Flow[19]+Demand[28]-Supply[1]=0;
subject to node22:
-1*Flow[19]+Demand[30]+Demand[31]=0;
subject to node23:
Flow[20]+Demand[36]+Demand[37]-Supply[6]=0;
subject to node24:
-1*Flow[20]+Flow[21]-Flow[33]-Flow[46]=0;
subject to node25:
-1*Flow[21]+Flow[22]+Demand[33]=0;
subject to node26:
-1*Flow[22]+Flow[23]+Demand[34]=0;
subject to node27:
-1*Flow[23]-Flow[24]+Flow[34]=0;
subject to node28:
Flow[24]-Flow[25]+Flow[26]+Flow[39]+Demand[19]=0;
subject to node29:
Flow[25]-Flow[52]+Flow[53]=0;
subject to node30:
-1*Flow[26]-Flow[30]+Flow[31]-Flow[39]-Flow[43]+Flow[44]+Demand[18]=0;
subject to node31:
-1*Flow[27]+Flow[28]-Flow[40]+Flow[41]+Demand[29]=0;
subject to node32:
-1*Flow[28]+Flow[29]-Flow[41]+Flow[42]=0;
subject to node33:
-1*Flow[29]+Flow[30]-Flow[42]+Flow[43]+Demand[32]=0;
subject to node34:
-1*Flow[31]+Flow[32]-Flow[34]-Flow[44]+Flow[45]=0;
subject to node35:
-1*Flow[32]+Flow[33]-Flow[45]+Flow[46]+Demand[35]=0;
subject to node36:
Flow[35]-Flow[38]+Flow[47]+Demand[11]=0;
subject to node37:
Flow[37]+Flow[49]-Supply[7]=0;
subject to node38:
-1*Flow[37]+Flow[38]+Demand[2]=0;
subject to node39:
-1*Flow[49]+Flow[50]+Demand[14]=0;
subject to node40:
-1*Flow[50]-Flow[51]+Flow[52]+Demand[15]=0;
subject to node41:
Flow[51]+Demand[17]-Supply[3]=0;
subject to node42:
-1*Flow[53]+Flow[54]+Demand[12]=0;
subject to node43:
-1*Flow[55]+Flow[56]+Demand[20]=0;
subject to node44:
Flow[57]+Flow[59]+Demand[1]-Supply[2]=0;
subject to node45:
-1*Flow[57]+Flow[58]+Demand[3]=0;
subject to node46:
-1*Flow[58]+Demand[4]=0;
subject to node47:
-1*Flow[59]+Flow[60]+Demand[13]=0;
subject to node48:
-1*Flow[60]+Demand[16]=0;
subject to arc1:
(Flow[1]*abs(Flow[1]))-(0.697133)*(Pi[1]-Pi[2])=0;
subject to arc2:
(Flow[2]*abs(Flow[2]))-(0.00265927)*(Pi[3]-Pi[4])=0;
subject to arc3:
(Flow[3]*abs(Flow[3]))-(0.172371)*(Pi[5]-Pi[6])=0;
subject to arc4:
(Flow[4]*abs(Flow[4]))-(0.258556)*(Pi[6]-Pi[7])=0;
subject to arc5:
(Flow[5]*abs(Flow[5]))-(0.190020)*(Pi[8]-Pi[9])=0;
subject to arc6:
(Flow[6]*abs(Flow[6]))-(4.18044)*(Pi[9]-Pi[10])=0;
subject to arc7:
(Flow[7]*abs(Flow[7]))-(0.321572)*(Pi[10]-Pi[11])=0;
subject to arc8:
(Flow[8]*abs(Flow[8]))+(0.0432263)*(Pi[10]-Pi[12])=0;
subject to arc9:
(Flow[9]*abs(Flow[9]))+(0.00611442)*(Pi[13]-Pi[8])=0;
subject to arc10:
(Flow[10]*abs(Flow[10]))+(0.0116174)*(Pi[12]-Pi[13])=0;
subject to arc11:
(Flow[11]*abs(Flow[11]))-(0.259358)*(Pi[12]-Pi[14])=0;
subject to arc12:
(Flow[12]*abs(Flow[12]))-(0.259358)*(Pi[14]-Pi[15])=0;
subject to arc13:
(Flow[13]*abs(Flow[13]))-(0.0285181)*(Pi[16]-Pi[17])=0;
subject to arc14:
(Flow[14]*abs(Flow[14]))-(0.0399253)*(Pi[17]-Pi[18])=0;
subject to arc15:
(Flow[15]*abs(Flow[15]))-(0.0285181)*(Pi[18]-Pi[15])=0;
subject to arc16:
(Flow[16]*abs(Flow[16]))+(0.0181479)*(Pi[1]-Pi[15])=0;
subject to arc17:
(Flow[17]*abs(Flow[17]))-(0.0864525)*(Pi[19]-Pi[20])=0;
subject to arc18:
(Flow[18]*abs(Flow[18]))-(0.0741022)*(Pi[20]-Pi[21])=0;
subject to arc19:
(Flow[19]*abs(Flow[19]))-(0.190020)*(Pi[21]-Pi[22])=0;
subject to arc20:
(Flow[20]*abs(Flow[20]))-(0.0235780)*(Pi[23]-Pi[24])=0;
subject to arc21:
(Flow[21]*abs(Flow[21]))-(4.18044)*(Pi[24]-Pi[25])=0;
subject to arc22:
(Flow[22]*abs(Flow[22]))-(0.522555)*(Pi[25]-Pi[26])=0;
subject to arc23:
(Flow[23]*abs(Flow[23]))-(1.39348)*(Pi[26]-Pi[27])=0;
subject to arc24:
(Flow[24]*abs(Flow[24]))+(0.696739)*(Pi[27]-Pi[28])=0;
subject to arc25:
(Flow[25]*abs(Flow[25]))+(0.418044)*(Pi[28]-Pi[29])=0;
subject to arc26:
(Flow[26]*abs(Flow[26]))-(0.000633373)*(Pi[28]-Pi[30])=0;
subject to arc27:
(Flow[27]*abs(Flow[27]))+(0.0332711)*(Pi[31]-Pi[19])=0;
subject to arc28:
(Flow[28]*abs(Flow[28]))-(0.0998133)*(Pi[31]-Pi[32])=0;
subject to arc29:
(Flow[29]*abs(Flow[29]))-(0.0332711)*(Pi[32]-Pi[33])=0;
subject to arc30:
(Flow[30]*abs(Flow[30]))-(0.0332711)*(Pi[33]-Pi[30])=0;
subject to arc31:
(Flow[31]*abs(Flow[31]))-(0.00954957)*(Pi[30]-Pi[34])=0;
subject to arc32:
(Flow[32]*abs(Flow[32]))-(0.0668470)*(Pi[34]-Pi[35])=0;
subject to arc33:
(Flow[33]*abs(Flow[33]))-(0.00607700)*(Pi[35]-Pi[24])=0;
subject to arc34:
(Flow[34]*abs(Flow[34]))-(0.410240)*(Pi[27]-Pi[34])=0;
subject to arc35:
(Flow[35]*abs(Flow[35]))-(0.0998133)*(Pi[36]-Pi[7])=0;
subject to arc36:
(Flow[36]*abs(Flow[36]))-(0.0499067)*(Pi[7]-Pi[16])=0;
subject to arc37:
(Flow[37]*abs(Flow[37]))-(0.348370)*(Pi[37]-Pi[38])=0;
subject to arc38:
(Flow[38]*abs(Flow[38]))-(0.696739)*(Pi[38]-Pi[36])=0;
subject to arc39:
(Flow[39]*abs(Flow[39]))-(0.0199627)*(Pi[28]-Pi[30])=0;
subject to arc40:
(Flow[40]*abs(Flow[40]))+(0.0193623)*(Pi[31]-Pi[19])=0;
subject to arc41:
(Flow[41]*abs(Flow[41]))-(0.0580870)*(Pi[31]-Pi[32])=0;
subject to arc42:
(Flow[42]*abs(Flow[42]))-(0.0193623)*(Pi[32]-Pi[33])=0;
subject to arc43:
(Flow[43]*abs(Flow[43]))-(0.0193623)*(Pi[33]-Pi[30])=0;
subject to arc44:
(Flow[44]*abs(Flow[44]))-(0.00954957)*(Pi[30]-Pi[34])=0;
subject to arc45:
(Flow[45]*abs(Flow[45]))-(0.0668470)*(Pi[34]-Pi[35])=0;
subject to arc46:
(Flow[46]*abs(Flow[46]))-(0.00607700)*(Pi[35]-Pi[24])=0;
subject to arc47:
(Flow[47]*abs(Flow[47]))-(2.09022)*(Pi[36]-Pi[7])=0;
subject to arc48:
(Flow[48]*abs(Flow[48]))-(0.696739)*(Pi[7]-Pi[16])=0;
subject to arc49:
(Flow[49]*abs(Flow[49]))-(0.0143027)*(Pi[37]-Pi[39])=0;
subject to arc50:
(Flow[50]*abs(Flow[50]))-(0.109654)*(Pi[39]-Pi[40])=0;
subject to arc51:
(Flow[51]*abs(Flow[51]))-(0.522555)*(Pi[41]-Pi[40])=0;
subject to arc52:
(Flow[52]*abs(Flow[52]))-(0.321572)*(Pi[40]-Pi[29])=0;
subject to arc53:
(Flow[53]*abs(Flow[53]))-(0.164096)*(Pi[29]-Pi[42])=0;
subject to arc54:
(Flow[54]*abs(Flow[54]))-(0.136747)*(Pi[42]-Pi[16])=0;
subject to arc55:
(Flow[55]*abs(Flow[55]))-(0.0399253)*(Pi[3]-Pi[43])=0;
subject to arc56:
(Flow[56]*abs(Flow[56]))+(0.0285181)*(Pi[1]-Pi[43])=0;
subject to arc57:
(Flow[57]*abs(Flow[57]))-(0.182329)*(Pi[44]-Pi[45])=0;
subject to arc58:
(Flow[58]*abs(Flow[58]))-(0.182329)*(Pi[45]-Pi[46])=0;
subject to arc59:
(Flow[59]*abs(Flow[59]))-(0.0332711)*(Pi[44]-Pi[47])=0;
subject to arc60:
(Flow[60]*abs(Flow[60]))-(0.00998133)*(Pi[47]-Pi[48])=0;
subject to region1:
-1*Demand[1]<=Region[1];
subject to region2:
-1*Demand[2]-Demand[3]-Demand[4]<=Region[2];
subject to region3:
-1*Demand[5]-Demand[6]<=Region[3];
subject to region4:
-1*Demand[7]-Demand[8]-Demand[9]-Demand[10]-Demand[11]<=Region[4];
subject to region5:
-1*Demand[12]-Demand[13]<=Region[5];
subject to region6:
-1*Demand[14]-Demand[15]-Demand[16]<=Region[6];
subject to region7:
-1*Demand[17]<=Region[7];
subject to region8:
-1*Demand[18]<=Region[8];
subject to region9:
-1*Demand[19]<=Region[9];
subject to region10:
-1*Demand[20]<=Region[10];
subject to region11:
-1*Demand[21]<=Region[11];
subject to region12:
-1*Demand[22]-Demand[23]<=Region[12];
subject to region13:
-1*Demand[24]-Demand[25]-Demand[26]<=Region[13];
subject to region14:
-1*Demand[27]-Demand[28]<=Region[14];
subject to region15:
-1*Demand[29]<=Region[15];
subject to region16:
-1*Demand[30]<=Region[16];
subject to region17:
-1*Demand[31]<=Region[17];
subject to region18:
-1*Demand[32]<=Region[18];
subject to region19:
-1*Demand[33]-Demand[34]-Demand[35]<=Region[19];
subject to region20:
-1*Demand[36]<=Region[20];
subject to region21:
-1*Demand[37]<=Region[21];
subject to prod1:
Supply[1]-Prod[1]<=0;
subject to prod2:
Supply[2]-Prod[2]<=0;
subject to prod3:
Supply[3]-Prod[3]<=0;
subject to prod4:
Supply[4]+Supply[5]+Supply[6]-Prod[4]<=0;
subject to prod5:
Supply[7]-Prod[5]<=0;

data;
param Region:=
1 	-8.90000E+00
2 	-3.30000E+00
3 	-5.40000E+00
4 	-1.45000E+01
5 	-2.20000E+00
6 	-5.80000E+00
7 	-3.50000E+00
8 	-1.70000E+00
9 	-1.50000E+00
10 	-8.00000E-01
11 	-1.30000E+00
12 	-2.30000E+00
13 	-2.20000E+00
14 	-1.10000E+00
15 	-7.00000E-01
16 	-1.50000E+00
17 	-4.30000E+00
18 	-1.50000E+00
19 	-4.30000E+00
20 	-1.10000E+00
21 	-5.00000E+00;

param loflow:=-361.167;
param hiflow:=361.167;

param hisupply:=
1 	1.00000E+30
2 	1.00000E+30
3 	1.00000E+30
4 	4.80000E+00
5 	2.45000E+01
6 	1.32000E+01
7 	1.00000E+30;

param loprod:=
1 	2.33333E+00
2 	4.92222E+00
3 	1.31333E+01
4 	9.55556E+00
5 	6.06667E+00;

param hiprod:=
1 	6.41667E+00
2 	1.35361E+01
3 	3.61167E+01
4 	2.62778E+01
5 	1.66833E+01;

param lopi:=
1 	0.00000E+00
2 	1.60000E+03
3 	0.00000E+00
4 	1.60000E+03
5 	1.60000E+03
6 	1.60000E+03
7 	0.00000E+00
8 	1.60000E+03
9 	1.60000E+03
10 	0.00000E+00
11 	1.60000E+03
12 	0.00000E+00
13 	1.60000E+03
14 	1.60000E+03
15 	0.00000E+00
16 	1.60000E+03
17 	1.60000E+03
18 	1.60000E+03
19 	0.00000E+00
20 	1.60000E+03
21 	1.60000E+03
22 	1.60000E+03
23 	1.60000E+03
24 	0.00000E+00
25 	1.60000E+03
26 	1.60000E+03
27 	0.00000E+00
28 	1.60000E+03
29 	0.00000E+00
30 	1.60000E+03
31 	1.60000E+03
32 	0.00000E+00
33 	1.60000E+03
34 	0.00000E+00
35 	1.60000E+03
36 	1.60000E+03
37 	0.00000E+00
38 	1.60000E+03
39 	1.60000E+03
40 	1.60000E+03
41 	1.60000E+03
42 	1.60000E+03
43 	1.60000E+03
44 	1.60000E+03
45 	1.60000E+03
46 	1.60000E+03
47 	1.60000E+03
48 	1.60000E+03;

var Flow:=
1 	2.28416E+00
2 	1.28924E+00
3 	-4.11076E+00
4 	-4.11076E+00
5 	8.34083E+00
6 	6.67356E+00
7 	5.15250E+00
8 	-1.52106E+00
9 	1.49619E+00
10 	8.63464E-01
11 	2.38452E+00
12 	1.84522E-01
13 	-2.80167E-01
14 	-2.80167E-01
15 	-2.80167E-01
16 	-9.56457E-02
17 	4.83333E-01
18 	0.00000E+00
19 	5.80000E+00
20 	3.62164E+00
21 	3.15920E+00
22 	0.00000E+00
23 	-1.14080E+00
24 	2.96863E+00
25 	6.56169E+00
26 	3.16455E-01
27 	-2.74176E-01
28 	-6.71257E-01
29 	-6.71257E-01
30 	-1.52215E+00
31 	-1.14514E+00
32 	-2.31222E-01
33 	-2.31221E-01
34 	1.82783E+00
35 	4.80828E-01
36 	-3.01828E-01
37 	1.06084E+01
38 	9.74452E+00
39 	1.77661E+00
40 	-2.09158E-01
41 	-5.12076E-01
42 	-5.12076E-01
43 	-1.16119E+00
44 	-1.14514E+00
45 	-2.31221E-01
46 	-2.31221E-01
47 	2.20035E+00
48 	-1.12775E+00
49 	1.91798E+00
50 	0.00000E+00
51 	1.15931E+01
52 	7.71110E+00
53 	1.14942E+00
54 	1.14942E+00
55 	3.17981E+00
56 	2.37981E+00
57 	2.43611E+00
58 	0.00000E+00
59 	2.20000E+00
60 	0.00000E+00;

param hipi:=
1 	6.40000E+03
2 	6.40000E+03
3 	6.40000E+03
4 	6.40000E+03
5 	6.40000E+03
6 	6.40000E+03
7 	6.40000E+03
8 	6.40000E+03
9 	6.40000E+03
10	 6.40000E+03
11	 6.40000E+03
12	 6.40000E+03
13	 6.40000E+03
14	 6.40000E+03
15	 6.40000E+03
16	 6.40000E+03
17	 6.40000E+03
18	 6.40000E+03
19	 6.40000E+03
20	 6.40000E+03
21	 4.48900E+03
22	 6.40000E+03
23	 6.40000E+03
24	 6.40000E+03
25 	6.40000E+03
26 	6.40000E+03
27 	6.40000E+03
28 	6.40000E+03
29 	6.40000E+03
30 	6.40000E+03
31 	6.40000E+03
32 	6.40000E+03
33 	6.40000E+03
34 	6.40000E+03
35 	6.40000E+03
36 	6.40000E+03
37 	4.48900E+03
38 	6.40000E+03
39 	6.40000E+03
40 	6.40000E+03
41 	6.40000E+03
42 	6.40000E+03
43 	6.40000E+03
44 	4.48900E+03
45 	6.40000E+03
46 	6.40000E+03
47 	6.40000E+03
48 	6.40000E+03;

solve;
